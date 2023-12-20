use anyhow::Result;
use log;

/// Node in the graph, can be at infinity
#[derive(Copy, Clone)]
pub enum Node {
    /// Node at infinity
    Infinity,

    /// Index of a finite node
    Value(usize),
}

impl Node {
    /// Checks equality between nodes
    pub fn equals(&self, node: &Node) -> bool {
        match (self, node) {
            (Node::Infinity, Node::Infinity) => true,
            (Node::Value(v1), Node::Value(v2)) => v1 == v2,
            (_, _) => false,
        }
    }

    /// Node to string
    pub fn to_string(&self) -> String {
        match self {
            Node::Infinity => "Node Infinity".to_string(),
            Node::Value(val) => format!("Node {}", val),
        }
    }

    /// Print node string
    pub fn print(&self) -> () {
        print!("{}", self.to_string());
    }

    /// Println node string
    pub fn println(&self) -> () {
        println!("{}", self.to_string());
    }
}

/// 2D Simplicial structure
pub struct SimplicialStructure2D {
    // i   : he1 \
    // i+1 : he2  | -> triangle
    // i+2 : he3 /
    // such that he2 = next(he1)
    // such that he3 = next(he2)
    // such that he1 = next(he3)
    // only first node is stored (last node is the next one, taking account of the %3)
    halfedge_first_node: Vec<Node>,
    halfedge_opposite: Vec<usize>,

    nb_triangles: usize,
}

#[derive(Copy, Clone)]
/// Halfedge iterator
pub struct IterHalfEdge<'a> {
    simplicial: &'a SimplicialStructure2D,
    ind_halfedge: usize,
}

#[derive(Copy, Clone)]
/// Triangle iterator
pub struct IterTriangle<'a> {
    simplicial: &'a SimplicialStructure2D,
    ind_triangle: usize,
}

impl SimplicialStructure2D {
    /// Simplicial structure initialisation
    pub fn new() -> SimplicialStructure2D {
        SimplicialStructure2D {
            halfedge_first_node: Vec::new(),
            halfedge_opposite: Vec::new(),
            nb_triangles: 0,
        }
    }

    /// Gets halfedge iterator from index
    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if ind_halfedge < self.halfedge_first_node.len() {
            Ok(IterHalfEdge {
                simplicial: self,
                ind_halfedge,
            })
        } else {
            Err(anyhow::Error::msg("Halfedge value not in simplicial"))
        }
    }

    /// Gets triangle iterator from index
    pub fn get_triangle(&self, ind_triangle: usize) -> Result<IterTriangle> {
        if ind_triangle < self.nb_triangles {
            Ok(IterTriangle {
                simplicial: self,
                ind_triangle,
            })
        } else {
            Err(anyhow::Error::msg("Triangle value not in simplicial"))
        }
    }

    /// Gets number of triangles
    pub fn get_nb_triangles(&self) -> usize {
        self.nb_triangles
    }

    fn insert_triangle(&mut self, nod1: Node, nod2: Node, nod3: Node) -> (usize, usize, usize) {
        let ind_first = self.halfedge_first_node.len();
        self.halfedge_first_node.push(nod1);
        self.halfedge_first_node.push(nod2);
        self.halfedge_first_node.push(nod3);
        self.nb_triangles = self.nb_triangles + 1;

        (ind_first, ind_first + 1, ind_first + 2)
    }

    fn replace_triangle(
        &mut self,
        ind_tri: usize,
        nod1: Node,
        nod2: Node,
        nod3: Node,
    ) -> (usize, usize, usize) {
        let ind_first = ind_tri * 3;
        self.halfedge_first_node[ind_first] = nod1;
        self.halfedge_first_node[ind_first + 1] = nod2;
        self.halfedge_first_node[ind_first + 2] = nod3;

        (ind_first, ind_first + 1, ind_first + 2)
    }

    /// Inserts a first triangle in the structure
    pub fn first_triangle(&mut self, nodes: [usize; 3]) -> Result<[IterTriangle; 4]> {
        if self.nb_triangles != 0 {
            return Err(anyhow::Error::msg("Already triangles in simplicial"));
        }
        let n0 = Node::Value(nodes[0]);
        let n1 = Node::Value(nodes[1]);
        let n2 = Node::Value(nodes[2]);
        let ninf = Node::Infinity;
        let first_tri = self.nb_triangles;
        let (h01, h12, h20) = self.insert_triangle(n0, n1, n2);
        let (hi2, h21, h1i) = self.insert_triangle(ninf, n2, n1);
        let (h2i, hi0, h02) = self.insert_triangle(n2, ninf, n0);
        let (h10, h0i, hi1) = self.insert_triangle(n1, n0, ninf);

        self.halfedge_opposite.push(h10);
        self.halfedge_opposite.push(h21);
        self.halfedge_opposite.push(h02);
        self.halfedge_opposite.push(h2i);
        self.halfedge_opposite.push(h12);
        self.halfedge_opposite.push(hi1);
        self.halfedge_opposite.push(hi2);
        self.halfedge_opposite.push(h0i);
        self.halfedge_opposite.push(h20);
        self.halfedge_opposite.push(h01);
        self.halfedge_opposite.push(hi0);
        self.halfedge_opposite.push(h1i);

        Ok([
            IterTriangle {
                simplicial: self,
                ind_triangle: first_tri,
            },
            IterTriangle {
                simplicial: self,
                ind_triangle: first_tri + 1,
            },
            IterTriangle {
                simplicial: self,
                ind_triangle: first_tri + 2,
            },
            IterTriangle {
                simplicial: self,
                ind_triangle: first_tri + 3,
            },
        ])
    }

    /// Inserts a new node in a triangle
    pub fn insert_node_within_triangle(
        &mut self,
        node: usize,
        ind_tri: usize,
    ) -> Result<[IterTriangle; 3]> {
        if ind_tri > self.nb_triangles {
            return Err(anyhow::Error::msg("Triangle index out of bounds"));
        }
        let h01 = ind_tri * 3;
        let h12 = ind_tri * 3 + 1;
        let h20 = ind_tri * 3 + 2;

        let n0 = self.halfedge_first_node[h01];
        let n1 = self.halfedge_first_node[h12];
        let n2 = self.halfedge_first_node[h20];
        let nn = Node::Value(node);

        let h10 = self.halfedge_opposite[h01];
        let h21 = self.halfedge_opposite[h12];
        let h02 = self.halfedge_opposite[h20];

        let (h01, h1n, hn0) = self.replace_triangle(ind_tri, n0, n1, nn);
        let (h12, h2n, hn1) = self.insert_triangle(n1, n2, nn);
        let (h20, h0n, hn2) = self.insert_triangle(n2, n0, nn);

        self.halfedge_opposite[h10] = h01;
        self.halfedge_opposite[h21] = h12;
        self.halfedge_opposite[h02] = h20;
        self.halfedge_opposite[h01] = h10;
        self.halfedge_opposite[h1n] = hn1;
        self.halfedge_opposite[hn0] = h0n;
        self.halfedge_opposite.push(h21);
        self.halfedge_opposite.push(hn2);
        self.halfedge_opposite.push(h1n);
        self.halfedge_opposite.push(h02);
        self.halfedge_opposite.push(hn0);
        self.halfedge_opposite.push(h2n);

        Ok([
            IterTriangle {
                simplicial: self,
                ind_triangle: ind_tri,
            },
            IterTriangle {
                simplicial: self,
                ind_triangle: self.nb_triangles - 2,
            },
            IterTriangle {
                simplicial: self,
                ind_triangle: self.nb_triangles - 1,
            },
        ])
    }

    /// Flips halfedge
    pub fn flip_halfedge(&mut self, ind_he: usize) -> () {
        let ind_he_opp = self.halfedge_opposite[ind_he];
        let ind_tri1 = ind_he / 3;
        let ind_tri2 = ind_he_opp / 3;

        let h01 = ind_tri1 * 3;
        let h12 = ind_tri1 * 3 + 1;
        let h20 = ind_tri1 * 3 + 2;

        let h01_opp = ind_tri2 * 3;
        let h12_opp = ind_tri2 * 3 + 1;
        let h20_opp = ind_tri2 * 3 + 2;

        let (hab, hbc) = if h01 == ind_he {
            (h12, h20)
        } else if h12 == ind_he {
            (h20, h01)
        } else {
            (h01, h12)
        };

        let (hcd, hda) = if h01_opp == ind_he_opp {
            (h12_opp, h20_opp)
        } else if h12_opp == ind_he_opp {
            (h20_opp, h01_opp)
        } else {
            (h01_opp, h12_opp)
        };

        let na = self.halfedge_first_node[hab];
        let nb = self.halfedge_first_node[hbc];
        let nc = self.halfedge_first_node[hcd];
        let nd = self.halfedge_first_node[hda];

        let hba = self.halfedge_opposite[hab];
        let hcb = self.halfedge_opposite[hbc];
        let hdc = self.halfedge_opposite[hcd];
        let had = self.halfedge_opposite[hda];

        let (hbc, hcd, hdb) = self.replace_triangle(ind_tri1, nb, nc, nd);
        let (hda, hab, hbd) = self.replace_triangle(ind_tri2, nd, na, nb);

        self.halfedge_opposite[hab] = hba;
        self.halfedge_opposite[hda] = had;
        self.halfedge_opposite[hbc] = hcb;
        self.halfedge_opposite[hcd] = hdc;

        self.halfedge_opposite[hbd] = hdb;
        self.halfedge_opposite[hdb] = hbd;

        self.halfedge_opposite[hba] = hab;
        self.halfedge_opposite[had] = hda;
        self.halfedge_opposite[hcb] = hbc;
        self.halfedge_opposite[hdc] = hcd;
    }

    /// Checks validity of simplicial graph (unit tests purposes)
    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        for ind_he in 0..self.halfedge_first_node.len() {
            let he = self.get_halfedge(ind_he)?;
            valid = valid && he.is_valid();
        }

        Ok(valid)
    }

    /// Println each triangle of the graph
    pub fn println(&self) -> () {
        for ind_tri in 0..self.nb_triangles {
            let tri = IterTriangle {
                simplicial: self,
                ind_triangle: ind_tri,
            };
            print!("  ");
            tri.println();
        }
    }
}

impl<'a> IterHalfEdge<'a> {
    /// Gets halfedge index
    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    /// First node
    pub fn first_node(&self) -> Node {
        self.simplicial.halfedge_first_node[self.ind_halfedge]
    }

    /// Last node
    pub fn last_node(&self) -> Node {
        let on_fac = self.ind_halfedge % 3;

        let ind_next = if on_fac == 2 {
            self.ind_halfedge - 2
        } else {
            self.ind_halfedge + 1
        };
        self.simplicial.halfedge_first_node[ind_next]
    }

    /// Next halfedge on same triangle
    pub fn next_halfedge(&self) -> IterHalfEdge<'a> {
        let on_fac = self.ind_halfedge % 3;

        let ind_next = if on_fac == 2 {
            self.ind_halfedge - 2
        } else {
            self.ind_halfedge + 1
        };

        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_next,
        }
    }

    /// Previous halfedge on same triangle
    pub fn prev_halfedge(&self) -> IterHalfEdge<'a> {
        let on_fac = self.ind_halfedge % 3;

        let ind_prev = if on_fac == 0 {
            self.ind_halfedge + 2
        } else {
            self.ind_halfedge - 1
        };
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_prev,
        }
    }

    /// Opposite halfedge: Same vertices in opposite order (on neighbor triangle)
    pub fn opposite_halfedge(&self) -> IterHalfEdge<'a> {
        let ind_opp = self.simplicial.halfedge_opposite[self.ind_halfedge];
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_opp,
        }
    }

    /// Triangle containing halfedge
    pub fn triangle(&self) -> IterTriangle<'a> {
        let ind_triangle = self.ind_halfedge / 3;
        IterTriangle {
            simplicial: self.simplicial,
            ind_triangle,
        }
    }

    /// Checks halfedge validity (unit test purposes)
    pub fn is_valid(&self) -> bool {
        let first_node = self.first_node();
        let last_node = self.last_node();

        let he_next = self.next_halfedge();
        let he_prev = self.prev_halfedge();
        let he_opp = self.opposite_halfedge();

        let mut valid = true;

        if !he_next.first_node().equals(&last_node) {
            log::error!("{}: Wrong next halfedge", self.to_string());
            valid = false;
        }
        if !he_prev.last_node().equals(&first_node) {
            log::error!("{}: Wrong previous halfedge", self.to_string());
            valid = false;
        }
        if !he_opp.first_node().equals(&last_node) || !he_opp.last_node().equals(&first_node) {
            log::error!("{}: Wrong opposite halfedge", self.to_string());
            valid = false;
        }

        valid
    }

    /// Halfedge to string
    pub fn to_string(&self) -> String {
        format!(
            "Edge {}: {} -> {}",
            self.ind(),
            self.first_node().to_string(),
            self.last_node().to_string()
        )
    }

    /// Print halfedge string
    pub fn print(&self) -> () {
        print!("{}", self.to_string());
    }

    /// Println halfedge string
    pub fn println(&self) -> () {
        println!("{}", self.to_string());
    }
}

impl<'a> IterTriangle<'a> {
    /// Gets triangle index
    pub fn ind(&self) -> usize {
        self.ind_triangle
    }

    /// Returns true if one of the nodes is infinity
    pub fn contains_infinity(&self) -> bool {
        let [he0, he1, he2] = self.halfedges();

        he0.first_node().equals(&Node::Infinity)
            || he1.first_node().equals(&Node::Infinity)
            || he2.first_node().equals(&Node::Infinity)
    }

    /// Surrounding halfedges (array of halfedge iterators)
    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        [
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: self.ind_triangle * 3,
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: self.ind_triangle * 3 + 1,
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: self.ind_triangle * 3 + 2,
            },
        ]
    }

    /// Nodes(array of nodes)
    pub fn nodes(&self) -> [Node; 3] {
        [
            self.simplicial.halfedge_first_node[self.ind_triangle * 3],
            self.simplicial.halfedge_first_node[self.ind_triangle * 3 + 1],
            self.simplicial.halfedge_first_node[self.ind_triangle * 3 + 2],
        ]
    }

    /// Triangle to string
    pub fn to_string(&self) -> String {
        let nod1 = self.simplicial.halfedge_first_node[self.ind_triangle * 3];
        let nod2 = self.simplicial.halfedge_first_node[self.ind_triangle * 3 + 1];
        let nod3 = self.simplicial.halfedge_first_node[self.ind_triangle * 3 + 2];
        format!(
            "Face {}: {} -> {} -> {}",
            self.ind(),
            nod1.to_string(),
            nod2.to_string(),
            nod3.to_string()
        )
    }

    /// Print triangle string
    pub fn print(&self) -> () {
        print!("{}", self.to_string());
    }

    /// Println triangle string
    pub fn println(&self) -> () {
        println!("{}", self.to_string());
    }
}
