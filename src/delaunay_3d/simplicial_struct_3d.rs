use std::vec;

use anyhow::Result;
use log;

/// For each triangle index within tetrahedron, associate list of vertices within tetrahedron
pub const TRIANGLE_SUBINDICES: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// For each triangle index, for each halfedge index, associate triangle and halfedge index within
/// tetrahedron
pub const NEIGHBOR_HALFEDGE: [[(usize, usize); 3]; 4] = [
    [(2, 1), (1, 1), (3, 1)],
    [(3, 2), (0, 1), (2, 0)],
    [(1, 2), (0, 0), (3, 0)],
    [(2, 2), (0, 2), (1, 0)],
];

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

/// 3D Simplicial structure
pub struct SimplicialStructure3D {
    // i   : nod0 \
    // i+1 : nod1  | -> tetrahedron
    // i+2 : nod2  |
    // i+3 : nod3 /
    // such that tri0 = (i+1, i+3, i+2)
    // such that tri1 = (i,   i+2, i+3)
    // such that tri2 = (i,   i+3, i+1)
    // such that tri3 = (i,   i+1, i+2)
    tet_nodes: Vec<Node>,
    halftriangle_opposite: Vec<usize>,

    nb_tetrahedra: usize,

    // structures to speed up tetrahedra insertion with Bowyer Watson algorithm
    should_rem_tet: Vec<bool>,
    should_keep_tet: Vec<bool>,
    tet_to_rem: Vec<usize>,
    tet_to_keep: Vec<usize>,
    tet_to_check: Vec<usize>,
}

#[derive(Copy, Clone)]
/// Halfedge iterator
pub struct IterHalfEdge<'a> {
    simplicial: &'a SimplicialStructure3D,
    ind_halftriangle: usize,
    ind_halfedge: usize,
}

#[derive(Copy, Clone)]
/// Half triangle iterator
pub struct IterHalfTriangle<'a> {
    simplicial: &'a SimplicialStructure3D,
    ind_halftriangle: usize,
}

#[derive(Copy, Clone)]
/// Triangle iterator
pub struct IterTetrahedron<'a> {
    simplicial: &'a SimplicialStructure3D,
    ind_tetrahedron: usize,
}

impl SimplicialStructure3D {
    /// Simplicial structure initialisation
    pub fn new() -> SimplicialStructure3D {
        SimplicialStructure3D {
            tet_nodes: Vec::new(),
            halftriangle_opposite: Vec::new(),
            nb_tetrahedra: 0,
            should_rem_tet: Vec::new(),
            should_keep_tet: Vec::new(),
            tet_to_rem: Vec::new(),
            tet_to_keep: Vec::new(),
            tet_to_check: Vec::new(),
        }
    }

    fn halftriangle(&self, ind_halftriangle: usize) -> IterHalfTriangle {
        IterHalfTriangle {
            simplicial: self,
            ind_halftriangle,
        }
    }

    /// Gets halfedge iterator from index
    pub fn get_halftriangle(&self, ind_halftriangle: usize) -> Result<IterHalfTriangle> {
        if ind_halftriangle < self.halftriangle_opposite.len() {
            Ok(self.halftriangle(ind_halftriangle))
        } else {
            Err(anyhow::Error::msg("Halftriangle value not in simplicial"))
        }
    }

    fn tetrahedron(&self, ind_tetrahedron: usize) -> IterTetrahedron {
        IterTetrahedron {
            simplicial: self,
            ind_tetrahedron,
        }
    }

    /// Gets tetrahedron iterator from index
    pub fn get_tetrahedron(&self, ind_tetrahedron: usize) -> Result<IterTetrahedron> {
        if ind_tetrahedron < self.nb_tetrahedra {
            Ok(self.tetrahedron(ind_tetrahedron))
        } else {
            Err(anyhow::Error::msg("Tetrahedron value not in simplicial"))
        }
    }

    /// Gets number of triangles
    pub fn get_nb_tetrahedra(&self) -> usize {
        self.nb_tetrahedra
    }

    /// Gets tetrahedra containing a specific node
    pub fn get_tetrahedra_containing(&self, node: &Node) -> Vec<IterTetrahedron> {
        let mut vec_tet = Vec::new();
        for i in 0..self.nb_tetrahedra {
            let first_nod = i << 2;
            for j in 0..4 {
                if self.tet_nodes[first_nod + j].equals(node) {
                    vec_tet.push(self.tetrahedron(i));
                    break;
                }
            }
        }

        vec_tet
    }

    /// Starts BW insertion, setting a first tetrahedron to remove
    pub fn bw_start(&mut self, ind_first_tetra: usize) -> Result<()> {
        if self.tet_to_check.len() != 0 || self.tet_to_keep.len() != 0 {
            return Err(anyhow::Error::msg(
                "Bowyer Watson algorithm already started",
            ));
        }
        self.bw_rem_tetra(ind_first_tetra);

        Ok(())
    }

    /// Gets next tetrahedron to check
    pub fn bw_tetra_to_check(&mut self) -> Option<usize> {
        loop {
            if let Some(ind_tetra) = self.tet_to_check.pop() {
                if self.should_rem_tet[ind_tetra] == false
                    && self.should_keep_tet[ind_tetra] == false
                {
                    return Some(ind_tetra);
                }
            } else {
                break;
            }
        }
        None
    }

    /// Sets tetrahedron to remove
    pub fn bw_rem_tetra(&mut self, ind_tetra: usize) -> () {
        let tri0 = ind_tetra << 2;
        let tri1 = tri0 + 1;
        let tri2 = tri0 + 2;
        let tri3 = tri0 + 3;
        let opp_tri0 = self.halftriangle_opposite[tri0];
        let opp_tri1 = self.halftriangle_opposite[tri1];
        let opp_tri2 = self.halftriangle_opposite[tri2];
        let opp_tri3 = self.halftriangle_opposite[tri3];
        self.tet_to_check.push(opp_tri0 >> 2);
        self.tet_to_check.push(opp_tri1 >> 2);
        self.tet_to_check.push(opp_tri2 >> 2);
        self.tet_to_check.push(opp_tri3 >> 2);
        self.should_rem_tet[ind_tetra] = true;
        self.tet_to_rem.push(ind_tetra);
    }

    /// Sets tetrahedron to keep
    pub fn bw_keep_tetra(&mut self, ind_tetra: usize) -> Result<()> {
        self.should_keep_tet[ind_tetra] = true;
        self.tet_to_keep.push(ind_tetra);
        Ok(())
    }

    /// BW insertion algorithm
    pub fn bw_insert_node(&mut self, nod: Node) -> Result<Vec<usize>> {
        if self.tet_to_check.len() != 0 {
            return Err(anyhow::Error::msg(
                "Cannot insert node if all tetrahedra are not checked",
            ));
        }

        // 1 - find boundary triangle
        let ind_tri_first = if let Some(&ind_tetra_keep) = self.tet_to_keep.last() {
            let tetra = self.tetrahedron(ind_tetra_keep);
            let tris = tetra.halftriangles();
            if tris[0].opposite().tetrahedron().should_rem() {
                tris[0].ind()
            } else if tris[1].opposite().tetrahedron().should_rem() {
                tris[1].ind()
            } else if tris[2].opposite().tetrahedron().should_rem() {
                tris[2].ind()
            } else if tris[3].opposite().tetrahedron().should_rem() {
                tris[3].ind()
            } else {
                return Err(anyhow::Error::msg("Isolated kept tetrahedron"));
            }
        } else {
            return Err(anyhow::Error::msg("No kept tetrahedron"));
        };

        // 2 - build boundary triangles graph
        let mut vec_tri = vec![ind_tri_first];
        let mut vec_nei: Vec<[Option<usize>; 3]> = vec![[None; 3]];
        let mut ind_cur = 0;
        loop {
            let cur_tri = IterHalfTriangle {
                simplicial: &self,
                ind_halftriangle: vec_tri[ind_cur],
            };
            let he = cur_tri.halfedges();
            for j in 0..3 {
                if vec_nei[ind_cur][j].is_none() {
                    let mut he_cur = he[j].opposite().neighbor().opposite();
                    let (ind_cur2, j2) = loop {
                        if !he_cur.triangle().tetrahedron().should_rem() {
                            let ind_tri2 = he_cur.triangle().ind();
                            let j2 = he_cur.triangle_subind();
                            let ind_cur2 = if let Some((i2, _)) =
                                vec_tri.iter().enumerate().find(|(_, &ind)| ind == ind_tri2)
                            {
                                i2
                            } else {
                                vec_tri.push(ind_tri2);
                                vec_nei.push([None; 3]);
                                vec_tri.len() - 1
                            };
                            break (ind_cur2, j2);
                        } else {
                            he_cur = he_cur.neighbor().opposite();
                        }
                    };
                    vec_nei[ind_cur][j] = Some(ind_cur2);
                    vec_nei[ind_cur2][j2] = Some(ind_cur);
                }
            }
            ind_cur = ind_cur + 1;
            if ind_cur >= vec_tri.len() {
                break;
            }
        }

        let mut added_tets = Vec::new();
        // 3 - create tetrahedra
        for i in 0..vec_tri.len() {
            let cur_tri = IterHalfTriangle {
                simplicial: &self,
                ind_halftriangle: vec_tri[i],
            };
            let [nod0, nod1, nod2] = cur_tri.nodes();
            if let Some(ind_add) = self.tet_to_rem.pop() {
                added_tets.push(ind_add);
                self.replace_tetrahedron(ind_add, nod0, nod2, nod1, nod);
            } else {
                added_tets.push(self.get_nb_tetrahedra());
                self.halftriangle_opposite.push(0);
                self.halftriangle_opposite.push(0);
                self.halftriangle_opposite.push(0);
                self.halftriangle_opposite.push(0);
                self.insert_tetrahedron(nod0, nod2, nod1, nod);
            };
        }

        // 4 - create links
        for i in 0..vec_tri.len() {
            let (tri0, tri1, tri2, tri3) = (
                added_tets[i] * 4,
                added_tets[i] * 4 + 1,
                added_tets[i] * 4 + 2,
                added_tets[i] * 4 + 3,
            );

            let ind_tri_nei = vec_tri[i];

            let ind_nei0 = vec_nei[i][1].unwrap();
            let ind_nei1 = vec_nei[i][0].unwrap();
            let ind_nei2 = vec_nei[i][2].unwrap();

            let ind_tet_nei0 = added_tets[ind_nei0];
            let ind_tet_nei1 = added_tets[ind_nei1];
            let ind_tet_nei2 = added_tets[ind_nei2];

            let ind_tri0_nei = if vec_nei[ind_nei0][0] == Some(i) {
                ind_tet_nei0 * 4 + 1
            } else if vec_nei[ind_nei0][1] == Some(i) {
                ind_tet_nei0 * 4 + 0
            } else {
                ind_tet_nei0 * 4 + 2
            };
            let ind_tri1_nei = if vec_nei[ind_nei1][0] == Some(i) {
                ind_tet_nei1 * 4 + 1
            } else if vec_nei[ind_nei1][1] == Some(i) {
                ind_tet_nei1 * 4 + 0
            } else {
                ind_tet_nei1 * 4 + 2
            };
            let ind_tri2_nei = if vec_nei[ind_nei2][0] == Some(i) {
                ind_tet_nei2 * 4 + 1
            } else if vec_nei[ind_nei2][1] == Some(i) {
                ind_tet_nei2 * 4 + 0
            } else {
                ind_tet_nei2 * 4 + 2
            };

            self.halftriangle_opposite[tri0] = ind_tri0_nei;
            self.halftriangle_opposite[tri1] = ind_tri1_nei;
            self.halftriangle_opposite[tri2] = ind_tri2_nei;
            self.halftriangle_opposite[tri3] = ind_tri_nei;
            self.halftriangle_opposite[ind_tri_nei] = tri3;
        }

        loop {
            if let Some(ind_tetra_keep) = self.tet_to_keep.pop() {
                self.should_keep_tet[ind_tetra_keep] = false;
            } else {
                break;
            }
        }

        Ok(added_tets)
    }

    /// Clean removed tetraedra
    pub fn clean_to_rem(&mut self) -> Result<()> {
        self.tet_to_rem.sort();
        loop {
            if let Some(ind_tetra_rem) = self.tet_to_rem.pop() {
                self.should_rem_tet[ind_tetra_rem] = false;
                self.mov_end_tetrahedron(ind_tetra_rem)?;
            } else {
                break;
            }
        }
        Ok(())
    }

    fn insert_tetrahedron(
        &mut self,
        nod1: Node,
        nod2: Node,
        nod3: Node,
        nod4: Node,
    ) -> (usize, usize, usize, usize) {
        let ind_first = self.tet_nodes.len();
        self.tet_nodes.push(nod1);
        self.tet_nodes.push(nod2);
        self.tet_nodes.push(nod3);
        self.tet_nodes.push(nod4);
        self.should_rem_tet.push(false);
        self.should_keep_tet.push(false);
        self.nb_tetrahedra = self.nb_tetrahedra + 1;

        (ind_first, ind_first + 1, ind_first + 2, ind_first + 3)
    }

    fn replace_tetrahedron(
        &mut self,
        ind_tetra: usize,
        nod1: Node,
        nod2: Node,
        nod3: Node,
        nod4: Node,
    ) -> (usize, usize, usize, usize) {
        let ind_first = ind_tetra * 4;
        self.tet_nodes[ind_first] = nod1;
        self.tet_nodes[ind_first + 1] = nod2;
        self.tet_nodes[ind_first + 2] = nod3;
        self.tet_nodes[ind_first + 3] = nod4;
        self.should_rem_tet[ind_tetra] = false;
        self.should_keep_tet[ind_tetra] = false;

        (ind_first, ind_first + 1, ind_first + 2, ind_first + 3)
    }

    fn mov_end_tetrahedron(&mut self, ind_tetra: usize) -> Result<()> {
        if ind_tetra != self.nb_tetrahedra - 1 {
            let ind_tri_opp1 = self.halftriangle_opposite[self.halftriangle_opposite.len() - 4];
            let ind_tri_opp2 = self.halftriangle_opposite[self.halftriangle_opposite.len() - 3];
            let ind_tri_opp3 = self.halftriangle_opposite[self.halftriangle_opposite.len() - 2];
            let ind_tri_opp4 = self.halftriangle_opposite[self.halftriangle_opposite.len() - 1];

            let [nod1, nod2, nod3, nod4] = self.tetrahedron(self.nb_tetrahedra - 1).nodes();

            let (ind_tri1, ind_tri2, ind_tri3, ind_tri4) =
                self.replace_tetrahedron(ind_tetra, nod1, nod2, nod3, nod4);

            self.halftriangle_opposite[ind_tri1] = ind_tri_opp1;
            self.halftriangle_opposite[ind_tri2] = ind_tri_opp2;
            self.halftriangle_opposite[ind_tri3] = ind_tri_opp3;
            self.halftriangle_opposite[ind_tri4] = ind_tri_opp4;

            self.halftriangle_opposite[ind_tri_opp1] = ind_tri1;
            self.halftriangle_opposite[ind_tri_opp2] = ind_tri2;
            self.halftriangle_opposite[ind_tri_opp3] = ind_tri3;
            self.halftriangle_opposite[ind_tri_opp4] = ind_tri4;
        }

        self.tet_nodes.pop();
        self.tet_nodes.pop();
        self.tet_nodes.pop();
        self.tet_nodes.pop();

        self.halftriangle_opposite.pop();
        self.halftriangle_opposite.pop();
        self.halftriangle_opposite.pop();
        self.halftriangle_opposite.pop();

        self.should_rem_tet.pop();
        self.should_keep_tet.pop();
        self.nb_tetrahedra = self.nb_tetrahedra - 1;

        Ok(())
    }

    /// Inserts a first tetrahedron in the structure
    pub fn first_tetrahedron(&mut self, nodes: [usize; 4]) -> Result<[IterTetrahedron; 4]> {
        if self.nb_tetrahedra != 0 {
            return Err(anyhow::Error::msg("Already tetrahedra in simplicial"));
        }
        let n0 = Node::Value(nodes[0]);
        let n1 = Node::Value(nodes[1]);
        let n2 = Node::Value(nodes[2]);
        let n3 = Node::Value(nodes[3]);
        let ni = Node::Infinity;
        let first_tetra = self.nb_tetrahedra;
        let (t132, t023, t031, t012) = self.insert_tetrahedron(n0, n1, n2, n3);
        let (t2i3, t13i, t1i2, t123) = self.insert_tetrahedron(n1, n2, n3, ni);
        let (t3i2, t02i, t0i3, t032) = self.insert_tetrahedron(n0, n3, n2, ni);
        let (t1i3, t03i, t0i1, t013) = self.insert_tetrahedron(n0, n1, n3, ni);
        let (t2i1, t01i, t0i2, t021) = self.insert_tetrahedron(n0, n2, n1, ni);

        self.halftriangle_opposite.push(t123); // t132
        self.halftriangle_opposite.push(t032); // t023
        self.halftriangle_opposite.push(t013); // t031
        self.halftriangle_opposite.push(t021); // t012
        self.halftriangle_opposite.push(t3i2); // t2i3
        self.halftriangle_opposite.push(t1i3); // t13i
        self.halftriangle_opposite.push(t2i1); // t1i2
        self.halftriangle_opposite.push(t132); // t123
        self.halftriangle_opposite.push(t2i3); // t3i2
        self.halftriangle_opposite.push(t0i2); // t02i
        self.halftriangle_opposite.push(t03i); // t0i3
        self.halftriangle_opposite.push(t023); // t032
        self.halftriangle_opposite.push(t13i); // t1i3
        self.halftriangle_opposite.push(t0i3); // t03i
        self.halftriangle_opposite.push(t01i); // t0i1
        self.halftriangle_opposite.push(t031); // t013
        self.halftriangle_opposite.push(t1i2); // t2i1
        self.halftriangle_opposite.push(t0i1); // t01i
        self.halftriangle_opposite.push(t02i); // t0i2
        self.halftriangle_opposite.push(t012); // t021

        Ok([
            IterTetrahedron {
                simplicial: self,
                ind_tetrahedron: first_tetra,
            },
            IterTetrahedron {
                simplicial: self,
                ind_tetrahedron: first_tetra + 1,
            },
            IterTetrahedron {
                simplicial: self,
                ind_tetrahedron: first_tetra + 2,
            },
            IterTetrahedron {
                simplicial: self,
                ind_tetrahedron: first_tetra + 3,
            },
        ])
    }

    /// Checks validity of simplicial graph (unit tests purposes)
    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        for ind_tetra in 0..self.nb_tetrahedra {
            let tetra = self.get_tetrahedron(ind_tetra)?;

            valid = valid && tetra.is_valid();
            for tri in tetra.halftriangles() {
                valid = valid && tri.is_valid();
                for he in tri.halfedges() {
                    valid = valid && he.is_valid();
                }
            }
        }

        Ok(valid)
    }

    /// Println each triangle of the graph
    pub fn println(&self) -> () {
        for ind_tetra in 0..self.nb_tetrahedra {
            let tetra = self.tetrahedron(ind_tetra);
            print!("  ");
            tetra.println();
        }
    }
}

impl<'a> IterHalfEdge<'a> {
    /// Gets subindex within triangle
    pub fn triangle_subind(&self) -> usize {
        self.ind_halfedge
    }

    /// First node
    pub fn first_node(&self) -> Node {
        let mod4 = self.ind_halftriangle % 4;
        let subdind = TRIANGLE_SUBINDICES[mod4];

        self.simplicial.tet_nodes[self.ind_halftriangle - mod4 + subdind[self.ind_halfedge]]
    }

    /// Last node
    pub fn last_node(&self) -> Node {
        let mod4 = self.ind_halftriangle % 4;
        let subdind = TRIANGLE_SUBINDICES[mod4];

        self.simplicial.tet_nodes
            [self.ind_halftriangle - mod4 + subdind[(self.ind_halfedge + 1) % 3]]
    }

    /// Next halfedge on same triangle
    pub fn next(&self) -> IterHalfEdge<'a> {
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halftriangle: self.ind_halftriangle,
            ind_halfedge: (self.ind_halfedge + 1) % 3,
        }
    }

    /// Previous halfedge on same triangle
    pub fn prev(&self) -> IterHalfEdge<'a> {
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halftriangle: self.ind_halftriangle,
            ind_halfedge: (self.ind_halfedge + 2) % 3,
        }
    }

    /// Opposite halfedge on opposite triangle
    pub fn opposite(&self) -> IterHalfEdge<'a> {
        let [he0, he1, he2] = self.triangle().opposite().halfedges();

        let last_node = self.last_node();
        if he0.first_node().equals(&last_node) {
            he0
        } else if he1.first_node().equals(&last_node) {
            he1
        } else {
            he2
        }
    }

    /// Opposite halfedge on neighbor triangle (same tetrahedron)
    pub fn neighbor(&self) -> IterHalfEdge<'a> {
        let mod_tri = self.ind_halftriangle % 4;

        let (neigh_tri, neigh_halfedge) = NEIGHBOR_HALFEDGE[mod_tri][self.ind_halfedge];

        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halftriangle: self.ind_halftriangle - mod_tri + neigh_tri,
            ind_halfedge: neigh_halfedge,
        }
    }

    /// Triangle containing halfedge
    pub fn triangle(&self) -> IterHalfTriangle<'a> {
        IterHalfTriangle {
            simplicial: self.simplicial,
            ind_halftriangle: self.ind_halftriangle,
        }
    }

    /// Checks halfedge validity (unit test purposes)
    pub fn is_valid(&self) -> bool {
        let first_node = self.first_node();
        let last_node = self.last_node();

        let he_next = self.next();
        let he_prev = self.prev();
        let he_opp = self.opposite();
        let he_nei = self.neighbor();

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
        if !he_nei.first_node().equals(&last_node) || !he_nei.last_node().equals(&first_node) {
            log::error!("{}: Wrong neighbor halfedge", self.to_string());
            valid = false;
        }

        valid
    }

    /// Halfedge to string
    pub fn to_string(&self) -> String {
        format!(
            "Edge: {} -> {}",
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

impl<'a> IterHalfTriangle<'a> {
    /// Gets half triangle index
    pub fn ind(&self) -> usize {
        self.ind_halftriangle
    }

    /// Returns true if one of the nodes is infinity
    pub fn contains_infinity(&self) -> bool {
        let [nod1, nod2, nod3] = self.nodes();

        nod1.equals(&Node::Infinity) || nod2.equals(&Node::Infinity) || nod3.equals(&Node::Infinity)
    }

    /// Tetrahedron containing halftriangle
    pub fn tetrahedron(&self) -> IterTetrahedron<'a> {
        IterTetrahedron {
            simplicial: self.simplicial,
            ind_tetrahedron: self.ind_halftriangle >> 2,
        }
    }

    /// Surrounding halfedges (array of halfedge iterators)
    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        [
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halftriangle: self.ind_halftriangle,
                ind_halfedge: 0,
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halftriangle: self.ind_halftriangle,
                ind_halfedge: 1,
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halftriangle: self.ind_halftriangle,
                ind_halfedge: 2,
            },
        ]
    }

    /// Nodes(array of nodes)
    pub fn nodes(&self) -> [Node; 3] {
        let mod4 = self.ind_halftriangle % 4;
        let subdind = TRIANGLE_SUBINDICES[mod4];
        [
            self.simplicial.tet_nodes[self.ind_halftriangle - mod4 + subdind[0]],
            self.simplicial.tet_nodes[self.ind_halftriangle - mod4 + subdind[1]],
            self.simplicial.tet_nodes[self.ind_halftriangle - mod4 + subdind[2]],
        ]
    }

    /// Opposite node on same tetrahedron
    pub fn opposite_node(&self) -> Node {
        self.simplicial.tet_nodes[self.ind()]
    }

    /// Opposite halftriangle on neighbor tetrahedron
    pub fn opposite(&self) -> IterHalfTriangle<'a> {
        IterHalfTriangle {
            simplicial: self.simplicial,
            ind_halftriangle: self.simplicial.halftriangle_opposite[self.ind_halftriangle],
        }
    }

    /// Checks halftriangle validity (unit test purposes)
    pub fn is_valid(&self) -> bool {
        let [nod0, nod1, nod2] = self.nodes();

        let [nod0o, nod1o, nod2o] = self.opposite().nodes();

        if nod0.equals(&nod0o) && nod1.equals(&nod2o) && nod2.equals(&nod1o) {
            ()
        } else if nod0.equals(&nod2o) && nod1.equals(&nod1o) && nod2.equals(&nod0o) {
            ()
        } else if nod0.equals(&nod1o) && nod1.equals(&nod0o) && nod2.equals(&nod2o) {
            ()
        } else {
            log::error!("{}: Wrong opposite halftriangle", self.to_string());
            log::error!("{}", self.opposite().to_string());
            return false;
        };

        true
    }

    /// Triangle to string
    pub fn to_string(&self) -> String {
        let [nod1, nod2, nod3] = self.nodes();
        format!(
            "Triangle {}: {} -> {} -> {}",
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

impl<'a> IterTetrahedron<'a> {
    /// Gets tetrahedron index
    pub fn ind(&self) -> usize {
        self.ind_tetrahedron
    }

    fn should_rem(&self) -> bool {
        self.simplicial.should_rem_tet[self.ind_tetrahedron]
    }

    fn bw_to_keep(&self) -> bool {
        self.simplicial.should_keep_tet[self.ind_tetrahedron]
    }

    /// Returns true if one of the nodes is infinity
    pub fn contains_infinity(&self) -> bool {
        let ind_first = self.ind_tetrahedron << 2;
        self.simplicial.tet_nodes[ind_first + 0].equals(&Node::Infinity)
            || self.simplicial.tet_nodes[ind_first + 1].equals(&Node::Infinity)
            || self.simplicial.tet_nodes[ind_first + 2].equals(&Node::Infinity)
            || self.simplicial.tet_nodes[ind_first + 3].equals(&Node::Infinity)
    }

    /// Surrounding halftriangles (array of halftriangle iterators)
    pub fn halftriangles(&self) -> [IterHalfTriangle<'a>; 4] {
        let ind_first = self.ind_tetrahedron << 2;
        [
            IterHalfTriangle {
                simplicial: self.simplicial,
                ind_halftriangle: ind_first,
            },
            IterHalfTriangle {
                simplicial: self.simplicial,
                ind_halftriangle: ind_first + 1,
            },
            IterHalfTriangle {
                simplicial: self.simplicial,
                ind_halftriangle: ind_first + 2,
            },
            IterHalfTriangle {
                simplicial: self.simplicial,
                ind_halftriangle: ind_first + 3,
            },
        ]
    }

    /// Nodes(array of nodes)
    pub fn nodes(&self) -> [Node; 4] {
        let ind_first = self.ind_tetrahedron << 2;
        [
            self.simplicial.tet_nodes[ind_first],
            self.simplicial.tet_nodes[ind_first + 1],
            self.simplicial.tet_nodes[ind_first + 2],
            self.simplicial.tet_nodes[ind_first + 3],
        ]
    }

    /// Checks validity of tetrahedron (for unit test purposes)
    pub fn is_valid(&self) -> bool {
        if self.should_rem() || self.bw_to_keep() {
            log::error!("{}: non cleaned tetrahedron", self.to_string());
            false
        } else {
            let [n0, n1, n2, n3] = self.nodes();

            if n0.equals(&n1)
                || n0.equals(&n2)
                || n0.equals(&n3)
                || n1.equals(&n2)
                || n1.equals(&n3)
                || n2.equals(&n3)
            {
                log::error!("{}: Wrong set of nodes", self.to_string());
                false
            } else {
                true
            }
        }
    }

    /// Tetrahedron to string
    pub fn to_string(&self) -> String {
        let [nod1, nod2, nod3, nod4] = self.nodes();
        format!(
            "Tetrahedron {}: {} -> {} -> {} -> {}",
            self.ind(),
            nod1.to_string(),
            nod2.to_string(),
            nod3.to_string(),
            nod4.to_string()
        )
    }

    /// Print tetrahedron string
    pub fn print(&self) -> () {
        print!("{}", self.to_string());
    }

    /// Println tetrahedron string
    pub fn println(&self) -> () {
        println!("{}", self.to_string());
    }
}
