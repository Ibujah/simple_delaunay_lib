use anyhow::Result;
use std::collections::HashSet;

#[derive(Copy, Clone)]
pub enum Node {
    Infinity,
    Value(usize),
}

impl Node {
    pub fn equals(&self, node: &Node) -> bool {
        match (self, node) {
            (Node::Infinity, Node::Infinity) => true,
            (Node::Value(v1), Node::Value(v2)) => v1 == v2,
            (_, _) => false,
        }
    }
}

/// 2D Simplicial structure
pub struct SimplicialStructure2D {
    halfedge_nodes: Vec<[Node; 2]>,
    halfedge_opposite: Vec<usize>,
    halfedge_next: Vec<usize>,
    halfedge_prev: Vec<usize>,
    halfedge_triangle: Vec<usize>,

    triangle_halfedges: Vec<[usize; 3]>,
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
    pub fn new() -> SimplicialStructure2D {
        SimplicialStructure2D {
            halfedge_nodes: Vec::new(),
            halfedge_opposite: Vec::new(),
            halfedge_next: Vec::new(),
            halfedge_prev: Vec::new(),
            halfedge_triangle: Vec::new(),

            triangle_halfedges: Vec::new(),
        }
    }

    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if ind_halfedge < self.halfedge_nodes.len() {
            Ok(IterHalfEdge {
                simplicial: self,
                ind_halfedge,
            })
        } else {
            Err(anyhow::Error::msg("Halfedge value not in simplicial"))
        }
    }

    pub fn get_triangle(&self, ind_triangle: usize) -> Result<IterTriangle> {
        if ind_triangle < self.triangle_halfedges.len() {
            Ok(IterTriangle {
                simplicial: self,
                ind_triangle,
            })
        } else {
            Err(anyhow::Error::msg("Triangle value not in simplicial"))
        }
    }

    fn insert_edge(&mut self, ind_points: [Node; 2]) -> [usize; 2] {
        let he1 = self.insert_halfedge(ind_points);
        let he2 = self.insert_halfedge([ind_points[1], ind_points[0]]);

        self.halfedge_opposite.insert(he1, he2);
        self.halfedge_opposite.insert(he2, he1);

        [he1, he2]
    }

    fn insert_triangle(&mut self, ind_halfedges: [usize; 3]) -> usize {
        self.triangle_halfedges
            .insert(self.ind_new_triangle, ind_halfedges);

        self.halfedge_triangle
            .insert(ind_halfedges[0], self.ind_new_triangle);
        self.halfedge_triangle
            .insert(ind_halfedges[1], self.ind_new_triangle);
        self.halfedge_triangle
            .insert(ind_halfedges[2], self.ind_new_triangle);

        self.halfedge_next
            .insert(ind_halfedges[0], ind_halfedges[1]);
        self.halfedge_next
            .insert(ind_halfedges[1], ind_halfedges[2]);
        self.halfedge_next
            .insert(ind_halfedges[2], ind_halfedges[0]);

        self.halfedge_prev
            .insert(ind_halfedges[0], ind_halfedges[2]);
        self.halfedge_prev
            .insert(ind_halfedges[1], ind_halfedges[0]);
        self.halfedge_prev
            .insert(ind_halfedges[2], ind_halfedges[1]);

        self.ind_new_triangle = self.ind_new_triangle + 1;
        self.ind_new_triangle - 1
    }

    fn remove_halfedge(&mut self, ind_halfedge: usize) -> () {
        self.halfedge_opposite.remove(&ind_halfedge);
        let nodes = self.halfedge_nodes.remove(&ind_halfedge).unwrap();

        if let Node::Value(ind_nod) = nodes[0] {
            self.node_halfedges
                .get_mut(&ind_nod)
                .unwrap()
                .remove(&ind_halfedge);
        } else {
            self.inf_halfedges.remove(&ind_halfedge);
        }
    }

    fn remove_triangle(&mut self, ind_triangle: usize) -> () {
        let he = self.triangle_halfedges.remove(&ind_triangle).unwrap();

        for ind_he in he {
            self.halfedge_triangle.remove(&ind_he);
            self.halfedge_next.remove(&ind_he);
            self.halfedge_prev.remove(&ind_he);
            let &ind_he_opp = self.halfedge_opposite.get(&ind_he).unwrap();
            if !self.halfedge_triangle.contains_key(&ind_he_opp) {
                self.remove_halfedge(ind_he);
                self.remove_halfedge(ind_he_opp);
            }
        }
    }

    pub fn first_triangle(&mut self, nodes: [usize; 3]) -> Result<[usize; 4]> {
        self.insert_node(nodes[0])?;
        self.insert_node(nodes[1])?;
        self.insert_node(nodes[2])?;

        let [hi0, h0i] = self.insert_edge([Node::Infinity, Node::Value(nodes[0])]);
        let [hi1, h1i] = self.insert_edge([Node::Infinity, Node::Value(nodes[1])]);
        let [hi2, h2i] = self.insert_edge([Node::Infinity, Node::Value(nodes[2])]);
        let [h01, h10] = self.insert_edge([Node::Value(nodes[0]), Node::Value(nodes[1])]);
        let [h12, h21] = self.insert_edge([Node::Value(nodes[1]), Node::Value(nodes[2])]);
        let [h20, h02] = self.insert_edge([Node::Value(nodes[2]), Node::Value(nodes[0])]);

        let tri012 = self.insert_triangle([h01, h12, h20]);
        let tri21i = self.insert_triangle([h21, h1i, hi2]);
        let tri02i = self.insert_triangle([h02, h2i, hi0]);
        let tri10i = self.insert_triangle([h10, h0i, hi1]);

        Ok([tri012, tri21i, tri02i, tri10i])
    }

    pub fn insert_node_within_path(
        &mut self,
        node: usize,
        path: &Vec<usize>,
        triangles_to_rem: &HashSet<usize>,
    ) -> Result<Vec<usize>> {
        for &ind_triangle in triangles_to_rem.iter() {
            self.remove_triangle(ind_triangle);
        }
        self.insert_node(node)?;
        let ind_node = Node::Value(node);
        let mut vec_tri = Vec::new();
        for &ind_halfedge in path.iter() {
            let &nodes = self.halfedge_nodes.get(&ind_halfedge).unwrap();
            let [ind_he1, _] = self.insert_edge([nodes[1], ind_node]);
            let ind_he2 = self.insert_halfedge([ind_node, nodes[0]]);
            let ind_tri = self.insert_triangle([ind_he1, ind_he2, ind_halfedge]);
            vec_tri.push(ind_tri);
        }
        Ok(vec_tri)
    }

    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        let node_inf = self.get_node_inf();
        valid = valid && node_inf.is_valid();

        for (&ind_node, _) in self.node_halfedges.iter() {
            let node = self.get_node(ind_node)?;
            valid = valid && node.is_valid();
        }

        for (&ind_he, _) in self.halfedge_nodes.iter() {
            let he = self.get_halfedge(ind_he)?;
            valid = valid && he.is_valid();
        }

        Ok(valid)
    }

    pub fn get_triangle_indices(&self) -> Vec<usize> {
        self.triangle_halfedges
            .iter()
            .map(|(&ind, _)| ind)
            .collect()
    }
}

impl<'a> IterNode<'a> {
    /// Gets node index
    pub fn ind(&self) -> Node {
        self.ind_node
    }

    /// Gets list of halfedges starting at this node
    pub fn halfedges(&self) -> Vec<IterHalfEdge<'a>> {
        if let Node::Value(ind_nod) = self.ind_node {
            self.simplicial.node_halfedges[&ind_nod]
                .iter()
                .fold(Vec::new(), |mut v, &x| {
                    v.push(IterHalfEdge {
                        simplicial: self.simplicial,
                        ind_halfedge: x,
                    });
                    v
                })
        } else {
            self.simplicial
                .inf_halfedges
                .iter()
                .fold(Vec::new(), |mut v, &x| {
                    v.push(IterHalfEdge {
                        simplicial: self.simplicial,
                        ind_halfedge: x,
                    });
                    v
                })
        }
    }

    pub fn is_valid(&self) -> bool {
        for he in self.halfedges() {
            if !he.first_node().ind().equals(&self.ind()) {
                match self.ind() {
                    Node::Infinity => println!("Non coherent halfedge on infinity node"),
                    Node::Value(val) => println!("Non coherent halfedge on node {}", val),
                }
                return false;
            }
        }
        true
    }

    pub fn print(&self) -> () {
        match self.ind() {
            Node::Infinity => print!("Node Infinity"),
            Node::Value(val) => print!("Node {}", val),
        }
    }

    pub fn println(&self) -> () {
        self.print();
        println!("");
    }
}

impl<'a> IterHalfEdge<'a> {
    /// Gets halfedge index
    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    /// Gets both nodes iterator
    pub fn nodes(&self) -> [IterNode; 2] {
        let ind_nodes = self.simplicial.halfedge_nodes[&self.ind_halfedge];
        [
            IterNode {
                simplicial: self.simplicial,
                ind_node: ind_nodes[0],
            },
            IterNode {
                simplicial: self.simplicial,
                ind_node: ind_nodes[1],
            },
        ]
    }

    /// First node iterator
    pub fn first_node(&self) -> IterNode<'a> {
        IterNode {
            simplicial: self.simplicial,
            ind_node: self.simplicial.halfedge_nodes[&self.ind_halfedge][0],
        }
    }

    /// Last node iterator
    pub fn last_node(&self) -> IterNode<'a> {
        IterNode {
            simplicial: self.simplicial,
            ind_node: self.simplicial.halfedge_nodes[&self.ind_halfedge][1],
        }
    }

    /// Next halfedge on same face
    pub fn next_halfedge(&self) -> IterHalfEdge<'a> {
        let &ind_next = self
            .simplicial
            .halfedge_next
            .get(&self.ind_halfedge)
            .unwrap();
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_next,
        }
    }

    /// Previous halfedge on same face
    pub fn prev_halfedge(&self) -> IterHalfEdge<'a> {
        let &ind_prev = self
            .simplicial
            .halfedge_prev
            .get(&self.ind_halfedge)
            .unwrap();
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_prev,
        }
    }

    /// Opposite halfedge: Same vertices in opposite order (on neighbor face)
    pub fn opposite_halfedge(&self) -> IterHalfEdge<'a> {
        let &ind_opp = self
            .simplicial
            .halfedge_opposite
            .get(&self.ind_halfedge)
            .unwrap();
        IterHalfEdge {
            simplicial: self.simplicial,
            ind_halfedge: ind_opp,
        }
    }

    /// Face containing halfedge
    pub fn face(&self) -> IterTriangle<'a> {
        let &ind_triangle = self
            .simplicial
            .halfedge_triangle
            .get(&self.ind_halfedge)
            .unwrap();
        IterTriangle {
            simplicial: self.simplicial,
            ind_triangle,
        }
    }

    pub fn is_valid(&self) -> bool {
        let mut in_first_node = false;
        let first_node = self.first_node();
        let last_node = self.last_node();
        for he in first_node.halfedges().iter() {
            if he.ind() == self.ind() {
                in_first_node = true;
                break;
            }
        }

        if !in_first_node {
            self.print();
            println!(": Wrong first node");
        }

        let mut valid = true;

        if self.simplicial.halfedge_prev.contains_key(&self.ind()) {
            if !self
                .prev_halfedge()
                .last_node()
                .ind()
                .equals(&first_node.ind())
            {
                self.print();
                println!(": Wrong previous halfedge");
                valid = false;
            }
        } else {
            self.print();
            println!(": No previous halfedge");
            valid = false;
        }

        if self.simplicial.halfedge_next.contains_key(&self.ind()) {
            if !self
                .next_halfedge()
                .first_node()
                .ind()
                .equals(&last_node.ind())
            {
                self.print();
                println!(": Wrong next halfedge");
                valid = false;
            }
        } else {
            self.print();
            println!(": No next halfedge");
            valid = false;
        }

        if self.simplicial.halfedge_opposite.contains_key(&self.ind()) {
            let he_opp = self.opposite_halfedge();
            if !he_opp.first_node().ind().equals(&last_node.ind())
                || !he_opp.last_node().ind().equals(&first_node.ind())
            {
                self.print();
                println!(": Wrong opposite halfedge");
                valid = false;
            }
        } else {
            self.print();
            println!(": No opposite halfedge");
            valid = false;
        }

        if self.simplicial.halfedge_triangle.contains_key(&self.ind()) {
            let tri = self.face();
            let [he1, he2, he3] = tri.halfedges();
            if he1.ind() != self.ind() && he2.ind() != self.ind() && he3.ind() != self.ind() {
                self.print();
                println!(": Wrong triangle");
                valid = false;
            }
        } else {
            self.print();
            println!(": No triangle");
            valid = false;
        }

        valid && in_first_node
    }

    pub fn print(&self) -> () {
        print!("Edge {}: ", self.ind());
        self.first_node().print();
        print!(" -> ");
        self.last_node().print();
    }

    pub fn println(&self) -> () {
        self.print();
        println!("");
    }
}

impl<'a> IterTriangle<'a> {
    /// Gets triangle index
    pub fn ind(&self) -> usize {
        self.ind_triangle
    }

    /// Surrounding halfedges (array of halfedge iterators)
    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        let &triangle = self
            .simplicial
            .triangle_halfedges
            .get(&self.ind_triangle)
            .unwrap();

        [
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: triangle[0],
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: triangle[1],
            },
            IterHalfEdge {
                simplicial: self.simplicial,
                ind_halfedge: triangle[2],
            },
        ]
    }

    pub fn print(&self) -> () {
        print!("Face {}: ", self.ind());
        let [he0, he1, he2] = self.halfedges();
        he0.first_node().print();
        print!(" -> ");
        he1.first_node().print();
        print!(" -> ");
        he2.first_node().print();
    }

    pub fn println(&self) -> () {
        self.print();
        println!("");
    }
}
