use anyhow::Result;
use robust::{self, Coord};
use std::time::Instant;

use log;

use super::geometry_operations_2d::{build_hilbert_curve, is_convex};
use super::simplicial_struct_2d::{self, Node, SimplicialStructure2D};

pub enum ExtendedTriangle {
    Triangle([[f64; 2]; 3]),
    Segment([[f64; 2]; 2]),
}

/// 2D Delaunay structure
pub struct DelaunayStructure2D {
    simpl_struct: simplicial_struct_2d::SimplicialStructure2D,
    vertex_coordinates: Vec<[f64; 2]>,
    indices_to_insert: Vec<usize>,
    inserted_indices: Vec<usize>,
}

impl DelaunayStructure2D {
    pub fn new() -> DelaunayStructure2D {
        DelaunayStructure2D {
            simpl_struct: simplicial_struct_2d::SimplicialStructure2D::new(),
            vertex_coordinates: Vec::new(),
            indices_to_insert: Vec::new(),
            inserted_indices: Vec::new(),
        }
    }

    pub fn add_vertex_to_insert(&mut self, to_insert: [f64; 2]) -> () {
        self.indices_to_insert.push(self.vertex_coordinates.len());
        self.vertex_coordinates.push(to_insert);
    }

    pub fn add_vertices_to_insert(&mut self, to_insert: &Vec<[f64; 2]>) -> () {
        for &vert in to_insert.iter() {
            self.indices_to_insert.push(self.vertex_coordinates.len());
            self.vertex_coordinates.push(vert);
        }
    }

    pub fn get_simplicial(&self) -> &SimplicialStructure2D {
        &self.simpl_struct
    }

    pub fn get_vertices(&self) -> &Vec<[f64; 2]> {
        &self.vertex_coordinates
    }

    pub fn get_extended_triangle(&self, ind_triangle: usize) -> Result<ExtendedTriangle> {
        let [node1, node2, node3] = self.get_simplicial().get_triangle(ind_triangle)?.nodes();

        let ext_tri = match (node1, node2, node3) {
            (
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt2 = self.get_vertices()[ind_v2];
                let pt3 = self.get_vertices()[ind_v3];
                ExtendedTriangle::Segment([pt2, pt3])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt3 = self.get_vertices()[ind_v3];
                ExtendedTriangle::Segment([pt3, pt1])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Infinity,
            ) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt2 = self.get_vertices()[ind_v2];
                ExtendedTriangle::Segment([pt1, pt2])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt2 = self.get_vertices()[ind_v2];
                let pt3 = self.get_vertices()[ind_v3];
                ExtendedTriangle::Triangle([pt1, pt2, pt3])
            }
            (_, _, _) => {
                return Err(anyhow::Error::msg("Case should not happen"));
            }
        };

        Ok(ext_tri)
    }

    fn is_vertex_strict_in_circle(&self, ind_vert: usize, ind_tri: usize) -> Result<bool> {
        let vert = self.get_vertices()[ind_vert];
        let ext_tri = self.get_extended_triangle(ind_tri)?;

        let sign = match ext_tri {
            ExtendedTriangle::Triangle(tri) => robust::incircle(
                Coord {
                    x: tri[0][0],
                    y: tri[0][1],
                },
                Coord {
                    x: tri[1][0],
                    y: tri[1][1],
                },
                Coord {
                    x: tri[2][0],
                    y: tri[2][1],
                },
                Coord {
                    x: vert[0],
                    y: vert[1],
                },
            ),
            ExtendedTriangle::Segment(lin) => robust::orient2d(
                Coord {
                    x: lin[0][0],
                    y: lin[0][1],
                },
                Coord {
                    x: lin[1][0],
                    y: lin[1][1],
                },
                Coord {
                    x: vert[0],
                    y: vert[1],
                },
            ),
        };
        Ok(sign > 0.)
    }

    fn is_triangle_flat(&self, ind_tri: usize) -> Result<bool> {
        let ext_tri = self.get_extended_triangle(ind_tri)?;

        let flat = if let ExtendedTriangle::Triangle(tri) = ext_tri {
            let sign = robust::orient2d(
                Coord {
                    x: tri[0][0],
                    y: tri[0][1],
                },
                Coord {
                    x: tri[1][0],
                    y: tri[1][1],
                },
                Coord {
                    x: tri[2][0],
                    y: tri[2][1],
                },
            );
            sign == 0.
        } else {
            false
        };
        Ok(flat)
    }

    fn choose_he<'a>(
        &self,
        vec_edg: &Vec<simplicial_struct_2d::IterHalfEdge<'a>>,
        vert: &[f64; 2],
    ) -> Option<simplicial_struct_2d::IterHalfEdge<'a>> {
        for &he in vec_edg {
            let ind1 = he.first_node();
            let ind2 = he.last_node();
            if let (Node::Value(v1), Node::Value(v2)) = (ind1, ind2) {
                let pt1 = self.get_vertices()[v1];
                let pt2 = self.get_vertices()[v2];
                let sign = robust::orient2d(
                    Coord {
                        x: pt1[0],
                        y: pt1[1],
                    },
                    Coord {
                        x: pt2[0],
                        y: pt2[1],
                    },
                    Coord {
                        x: vert[0],
                        y: vert[1],
                    },
                );
                if he.face().contains_infinity() {
                    if sign <= 0. {
                        return Some(he);
                    }
                } else if sign < 0. {
                    return Some(he);
                }
            }
        }
        None
    }

    fn walk_by_visibility(&self, ind_vert: usize, ind_starting_triangle: usize) -> Result<usize> {
        let vert = self.get_vertices()[ind_vert];
        let mut ind_tri_cur = ind_starting_triangle;
        let start_tri = self.get_simplicial().get_triangle(ind_tri_cur)?;
        let mut vec_edg: Vec<simplicial_struct_2d::IterHalfEdge> =
            start_tri.halfedges().iter().map(|&he| he).collect();
        let mut side = false;
        loop {
            if let Some(he) = self.choose_he(&vec_edg, &vert) {
                let he_opp = he.opposite_halfedge();
                ind_tri_cur = he_opp.face().ind();
                vec_edg.clear();
                if side {
                    vec_edg.push(he_opp.next_halfedge());
                    vec_edg.push(he_opp.prev_halfedge());
                } else {
                    vec_edg.push(he_opp.prev_halfedge());
                    vec_edg.push(he_opp.next_halfedge());
                }
                side = !side;
            } else {
                return Ok(ind_tri_cur);
            }
        }
    }

    fn should_flip_halfedge(&self, ind_he: usize) -> Result<bool> {
        let he = self.get_simplicial().get_halfedge(ind_he)?;
        let ind_tri_abd = he.face().ind();
        let node_a = he.prev_halfedge().first_node();
        let node_b = he.first_node();

        let ind_tri_bcd = he.opposite_halfedge().face().ind();
        let node_c = he.opposite_halfedge().prev_halfedge().first_node();
        let node_d = he.opposite_halfedge().first_node();

        match (node_a, node_b, node_c, node_d) {
            (Node::Value(ind_node_a), Node::Value(_), Node::Value(ind_node_c), Node::Value(_)) => {
                Ok(self.is_vertex_strict_in_circle(ind_node_c, ind_tri_abd)?
                    || self.is_vertex_strict_in_circle(ind_node_a, ind_tri_bcd)?)
            }
            (Node::Infinity, Node::Value(_), Node::Value(ind_node_c), Node::Value(_)) => {
                Ok(self.is_vertex_strict_in_circle(ind_node_c, ind_tri_abd)?
                    || self.is_triangle_flat(ind_tri_bcd)?)
            }
            (
                Node::Value(ind_node_a),
                Node::Infinity,
                Node::Value(ind_node_c),
                Node::Value(ind_node_d),
            ) => {
                let pt_a = self.get_vertices()[ind_node_a];
                let pt_c = self.get_vertices()[ind_node_c];
                let pt_d = self.get_vertices()[ind_node_d];
                Ok(is_convex(pt_c, pt_d, pt_a) == 1)
            }
            (Node::Value(ind_node_a), Node::Value(_), Node::Infinity, Node::Value(_)) => {
                Ok(self.is_triangle_flat(ind_tri_abd)?
                    || self.is_vertex_strict_in_circle(ind_node_a, ind_tri_bcd)?)
            }
            (
                Node::Value(ind_node_a),
                Node::Value(ind_node_b),
                Node::Value(ind_node_c),
                Node::Infinity,
            ) => {
                let pt_a = self.get_vertices()[ind_node_a];
                let pt_b = self.get_vertices()[ind_node_b];
                let pt_c = self.get_vertices()[ind_node_c];
                Ok(is_convex(pt_a, pt_b, pt_c) == 1)
            }
            (_, _, _, _) => Err(anyhow::Error::msg("Multiple infinity linked together")),
        }
    }

    pub fn update_delaunay(&mut self) -> Result<()> {
        if self.get_vertices().len() < 3 {
            return Err(anyhow::Error::msg(
                "Needs at least 3 vertices to compute Delaunay",
            ));
        }

        let now = Instant::now();
        self.indices_to_insert = build_hilbert_curve(self.get_vertices(), &self.indices_to_insert);
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        log::info!("Hilbert curve computed in {}ms", nano as f32 / 1e6);

        let now = Instant::now();
        // first triangle insertion
        if self.get_vertices().len() == self.indices_to_insert.len() {
            let ind1 = self.indices_to_insert.pop().unwrap();
            let ind2 = self.indices_to_insert.pop().unwrap();
            let pt1 = self.get_vertices()[ind1];
            let pt2 = self.get_vertices()[ind2];

            let mut aligned = Vec::new();

            loop {
                if let Some(ind3) = self.indices_to_insert.pop() {
                    let pt3 = self.get_vertices()[ind3];

                    let sign = robust::orient2d(
                        Coord {
                            x: pt1[0],
                            y: pt1[1],
                        },
                        Coord {
                            x: pt2[0],
                            y: pt2[1],
                        },
                        Coord {
                            x: pt3[0],
                            y: pt3[1],
                        },
                    );

                    if sign > 0. {
                        self.inserted_indices.push(ind1);
                        self.inserted_indices.push(ind2);
                        self.inserted_indices.push(ind3);
                        self.simpl_struct.first_triangle([ind1, ind2, ind3])?
                    } else if sign < 0. {
                        self.inserted_indices.push(ind1);
                        self.inserted_indices.push(ind2);
                        self.inserted_indices.push(ind3);
                        self.simpl_struct.first_triangle([ind1, ind3, ind2])?
                    } else {
                        aligned.push(ind3);
                        continue;
                    };
                } else {
                    return Err(anyhow::Error::msg(
                        "Could not find three non aligned points",
                    ));
                }

                break;
            }

            self.indices_to_insert.append(&mut aligned);
        }
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        log::info!("First triangle computed in {}ms", nano as f32 / 1e6);

        let mut walk_ms = 0;
        let mut insert_ms = 0;
        let mut flip_ms = 0;
        let mut step = 0;
        loop {
            if let Some(ind_v) = self.indices_to_insert.pop() {
                let now = Instant::now();
                let ind_triangle =
                    self.walk_by_visibility(ind_v, self.simpl_struct.get_nb_triangles() - 1)?;

                let duration = now.elapsed();
                let milli = duration.as_nanos();
                walk_ms = walk_ms + milli;

                step = step + 1;

                let now = Instant::now();
                let mut he_to_evaluate = Vec::new();
                let [he1, he2, he3] = self.simpl_struct.get_triangle(ind_triangle)?.halfedges();
                he_to_evaluate.push(he1.opposite_halfedge().ind());
                he_to_evaluate.push(he2.opposite_halfedge().ind());
                he_to_evaluate.push(he3.opposite_halfedge().ind());
                let _ = self
                    .simpl_struct
                    .insert_node_within_triangle(ind_v, ind_triangle)?;

                let duration = now.elapsed();
                let milli = duration.as_nanos();
                insert_ms = insert_ms + milli;

                let now = Instant::now();
                while let Some(ind_he) = he_to_evaluate.pop() {
                    if self.should_flip_halfedge(ind_he)? {
                        let he = self.get_simplicial().get_halfedge(ind_he)?;
                        let ind_he_add1 = he.prev_halfedge().opposite_halfedge().ind();
                        let ind_he_add2 = he.next_halfedge().opposite_halfedge().ind();
                        let ind_he_add3 = he
                            .opposite_halfedge()
                            .prev_halfedge()
                            .opposite_halfedge()
                            .ind();
                        let ind_he_add4 = he
                            .opposite_halfedge()
                            .next_halfedge()
                            .opposite_halfedge()
                            .ind();
                        self.simpl_struct.flip_halfedge(ind_he);
                        he_to_evaluate.push(ind_he_add1);
                        he_to_evaluate.push(ind_he_add2);
                        he_to_evaluate.push(ind_he_add3);
                        he_to_evaluate.push(ind_he_add4);
                    }
                }
                self.inserted_indices.push(ind_v);

                let duration = now.elapsed();
                let milli = duration.as_nanos();
                flip_ms = flip_ms + milli;
            } else {
                break;
            }
        }
        log::info!("Walks computed in {}ms", walk_ms as f32 / 1e6);
        log::info!("Insertions computed in {}ms", insert_ms as f32 / 1e6);
        log::info!("Flips computed in {}ms", flip_ms as f32 / 1e6);

        Ok(())
    }

    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        for ind_tri in 0..self.get_simplicial().get_nb_triangles() {
            if self.is_triangle_flat(ind_tri)? {
                log::error!("Flat triangle: ");
                self.get_simplicial().get_triangle(ind_tri)?.println();
                valid = false;
            }
            for &ind_vert in self.inserted_indices.iter() {
                let in_circle = self.is_vertex_strict_in_circle(ind_vert, ind_tri)?;
                if in_circle {
                    log::error!("Non Delaunay triangle: ");
                    self.get_simplicial().get_triangle(ind_tri)?.println();
                    valid = false;
                }
            }
        }

        Ok(valid)
    }
}
