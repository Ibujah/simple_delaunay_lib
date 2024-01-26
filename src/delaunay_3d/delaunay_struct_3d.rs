use anyhow::Result;
use robust::{insphere, orient3d, Coord3D};
use std::time::Instant;

use super::geometry_operations_3d::build_hilbert_curve_3d;
use super::simplicial_struct_3d::{IterHalfTriangle, Node, SimplicialStructure3D};

/// Extended tetrahedron, including point at infinity
pub enum ExtendedTetrahedron {
    /// Regular tetrahedron
    Tetrahedron([[f64; 3]; 4]),
    /// Tetrahedron with a point at infinity
    Triangle([[f64; 3]; 3]),
}

/// 3D Delaunay structure
pub struct DelaunayStructure3D {
    simpl_struct: SimplicialStructure3D,
    vertex_coordinates: Vec<[f64; 3]>,
    walk_ns: u128,
    insert_ns: u128,
}

impl DelaunayStructure3D {
    /// Delaunay structure initialisation
    pub fn new() -> DelaunayStructure3D {
        DelaunayStructure3D {
            simpl_struct: SimplicialStructure3D::new(),
            vertex_coordinates: Vec::new(),
            walk_ns: 0,
            insert_ns: 0,
        }
    }

    /// Gets simplicial structure
    pub fn get_simplicial(&self) -> &SimplicialStructure3D {
        &self.simpl_struct
    }

    /// Gets graph vertices
    pub fn get_vertices(&self) -> &Vec<[f64; 3]> {
        &self.vertex_coordinates
    }

    /// Gets extended tetrahedron from index
    pub fn get_extended_tetrahedron(&self, ind_tetrahedron: usize) -> Result<ExtendedTetrahedron> {
        let [node1, node2, node3, node4] = self
            .get_simplicial()
            .get_tetrahedron(ind_tetrahedron)?
            .nodes();

        let ext_tri = match (node1, node2, node3, node4) {
            (Node::Infinity, Node::Value(ind_v2), Node::Value(ind_v3), Node::Value(ind_v4)) => {
                let pt2 = self.get_vertices()[ind_v2];
                let pt3 = self.get_vertices()[ind_v3];
                let pt4 = self.get_vertices()[ind_v4];
                ExtendedTetrahedron::Triangle([pt2, pt4, pt3])
            }
            (Node::Value(ind_v1), Node::Infinity, Node::Value(ind_v3), Node::Value(ind_v4)) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt3 = self.get_vertices()[ind_v3];
                let pt4 = self.get_vertices()[ind_v4];
                ExtendedTetrahedron::Triangle([pt1, pt3, pt4])
            }
            (Node::Value(ind_v1), Node::Value(ind_v2), Node::Infinity, Node::Value(ind_v4)) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt2 = self.get_vertices()[ind_v2];
                let pt4 = self.get_vertices()[ind_v4];
                ExtendedTetrahedron::Triangle([pt1, pt4, pt2])
            }
            (Node::Value(ind_v1), Node::Value(ind_v2), Node::Value(ind_v3), Node::Infinity) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt2 = self.get_vertices()[ind_v2];
                let pt3 = self.get_vertices()[ind_v3];
                ExtendedTetrahedron::Triangle([pt1, pt2, pt3])
            }
            (
                Node::Value(ind_v1),
                Node::Value(ind_v2),
                Node::Value(ind_v3),
                Node::Value(ind_v4),
            ) => {
                let pt1 = self.get_vertices()[ind_v1];
                let pt2 = self.get_vertices()[ind_v2];
                let pt3 = self.get_vertices()[ind_v3];
                let pt4 = self.get_vertices()[ind_v4];
                ExtendedTetrahedron::Tetrahedron([pt1, pt2, pt3, pt4])
            }
            (_, _, _, _) => {
                return Err(anyhow::Error::msg("Case should not happen"));
            }
        };

        Ok(ext_tri)
    }

    fn is_vertex_in_sphere(&self, ind_vert: usize, ind_tetra: usize) -> Result<bool> {
        let vert = self.get_vertices()[ind_vert];
        let ext_tri = self.get_extended_tetrahedron(ind_tetra)?;

        let sign = match ext_tri {
            ExtendedTetrahedron::Tetrahedron(tri) => insphere(
                Coord3D {
                    x: tri[0][0],
                    y: tri[0][1],
                    z: tri[0][2],
                },
                Coord3D {
                    x: tri[1][0],
                    y: tri[1][1],
                    z: tri[1][2],
                },
                Coord3D {
                    x: tri[2][0],
                    y: tri[2][1],
                    z: tri[2][2],
                },
                Coord3D {
                    x: tri[3][0],
                    y: tri[3][1],
                    z: tri[3][2],
                },
                Coord3D {
                    x: vert[0],
                    y: vert[1],
                    z: vert[2],
                },
            ),
            ExtendedTetrahedron::Triangle(lin) => orient3d(
                Coord3D {
                    x: lin[0][0],
                    y: lin[0][1],
                    z: lin[0][2],
                },
                Coord3D {
                    x: lin[1][0],
                    y: lin[1][1],
                    z: lin[1][2],
                },
                Coord3D {
                    x: lin[2][0],
                    y: lin[2][1],
                    z: lin[2][2],
                },
                Coord3D {
                    x: vert[0],
                    y: vert[1],
                    z: vert[2],
                },
            ),
        };
        Ok(sign >= 0.)
    }

    fn is_vertex_strict_in_sphere(&self, ind_vert: usize, ind_tetra: usize) -> Result<bool> {
        let vert = self.get_vertices()[ind_vert];
        let ext_tri = self.get_extended_tetrahedron(ind_tetra)?;

        let sign = match ext_tri {
            ExtendedTetrahedron::Tetrahedron(tri) => insphere(
                Coord3D {
                    x: tri[0][0],
                    y: tri[0][1],
                    z: tri[0][2],
                },
                Coord3D {
                    x: tri[1][0],
                    y: tri[1][1],
                    z: tri[1][2],
                },
                Coord3D {
                    x: tri[2][0],
                    y: tri[2][1],
                    z: tri[2][2],
                },
                Coord3D {
                    x: tri[3][0],
                    y: tri[3][1],
                    z: tri[3][2],
                },
                Coord3D {
                    x: vert[0],
                    y: vert[1],
                    z: vert[2],
                },
            ),
            ExtendedTetrahedron::Triangle(lin) => orient3d(
                Coord3D {
                    x: lin[0][0],
                    y: lin[0][1],
                    z: lin[0][2],
                },
                Coord3D {
                    x: lin[1][0],
                    y: lin[1][1],
                    z: lin[1][2],
                },
                Coord3D {
                    x: lin[2][0],
                    y: lin[2][1],
                    z: lin[2][2],
                },
                Coord3D {
                    x: vert[0],
                    y: vert[1],
                    z: vert[2],
                },
            ),
        };
        Ok(sign > 0.)
    }

    fn is_tetrahedron_flat(&self, ind_tri: usize) -> Result<bool> {
        let ext_tri = self.get_extended_tetrahedron(ind_tri)?;

        let flat = if let ExtendedTetrahedron::Tetrahedron(tri) = ext_tri {
            let sign = orient3d(
                Coord3D {
                    x: tri[0][0],
                    y: tri[0][1],
                    z: tri[0][2],
                },
                Coord3D {
                    x: tri[1][0],
                    y: tri[1][1],
                    z: tri[1][2],
                },
                Coord3D {
                    x: tri[2][0],
                    y: tri[2][1],
                    z: tri[2][2],
                },
                Coord3D {
                    x: tri[3][0],
                    y: tri[3][1],
                    z: tri[3][2],
                },
            );
            sign == 0.
        } else {
            false
        };
        Ok(flat)
    }

    fn choose_tri<'a>(
        &self,
        vec_tri: &Vec<IterHalfTriangle<'a>>,
        vert: &[f64; 3],
    ) -> Option<IterHalfTriangle<'a>> {
        for &tri in vec_tri {
            let [nod1, nod2, nod3] = tri.nodes();
            if let (Node::Value(v1), Node::Value(v2), Node::Value(v3)) = (nod1, nod2, nod3) {
                let pt1 = self.get_vertices()[v1];
                let pt2 = self.get_vertices()[v2];
                let pt3 = self.get_vertices()[v3];
                let sign = orient3d(
                    Coord3D {
                        x: pt1[0],
                        y: pt1[1],
                        z: pt1[2],
                    },
                    Coord3D {
                        x: pt2[0],
                        y: pt2[1],
                        z: pt2[2],
                    },
                    Coord3D {
                        x: pt3[0],
                        y: pt3[1],
                        z: pt3[2],
                    },
                    Coord3D {
                        x: vert[0],
                        y: vert[1],
                        z: vert[2],
                    },
                );
                if tri.tetrahedron().contains_infinity() {
                    if sign <= 0. {
                        return Some(tri);
                    }
                } else if sign < 0. {
                    return Some(tri);
                }
            }
        }
        None
    }

    fn walk_check_all(&self, ind_vert: usize) -> Result<usize> {
        for ind_tetra_cur in 0..self.get_simplicial().get_nb_tetrahedra() {
            if self.is_tetrahedron_flat(ind_tetra_cur)? {
                continue;
            }
            if self.is_vertex_in_sphere(ind_vert, ind_tetra_cur)? {
                return Ok(ind_tetra_cur);
            }
        }
        Err(anyhow::Error::msg("Could not find sphere containing point"))
    }

    fn walk_by_visibility(
        &self,
        ind_vert: usize,
        ind_starting_tetrahedron: usize,
    ) -> Result<usize> {
        let vert = self.get_vertices()[ind_vert];
        let mut ind_tetra_cur = ind_starting_tetrahedron;
        let start_tetra = self.get_simplicial().get_tetrahedron(ind_tetra_cur)?;
        let mut vec_tri: Vec<IterHalfTriangle> =
            start_tetra.halftriangles().iter().map(|&tri| tri).collect();
        let mut side = 0;
        let mut nb_visited = 0;
        let th_visited = self.get_simplicial().get_nb_tetrahedra() >> 2;
        loop {
            if nb_visited > th_visited {
                break Err(anyhow::Error::msg("Could not find sphere containing point"));
            }
            if let Some(tri) = self.choose_tri(&vec_tri, &vert) {
                nb_visited = nb_visited + 1;
                let tri_opp = tri.opposite();
                ind_tetra_cur = tri_opp.tetrahedron().ind();
                vec_tri.clear();
                let hes = tri_opp.halfedges();
                vec_tri.push(hes[(0 + side) % 3].neighbor().triangle());
                vec_tri.push(hes[(1 + side) % 3].neighbor().triangle());
                vec_tri.push(hes[(2 + side) % 3].neighbor().triangle());
                side = (side + 1) % 3;
            } else {
                if self.is_vertex_in_sphere(ind_vert, ind_tetra_cur)? {
                    break Ok(ind_tetra_cur);
                } else {
                    break Err(anyhow::Error::msg("Could not find sphere containing point"));
                }
            }
        }
    }

    /// Inserts point using Bowyer Watson method
    fn insert_bw(&mut self, ind_vert: usize, ind_first_tetra: usize) -> Result<Vec<usize>> {
        self.simpl_struct.bw_start(ind_first_tetra)?;

        loop {
            if let Some(ind_tetra) = self.simpl_struct.bw_tetra_to_check() {
                if self.is_vertex_in_sphere(ind_vert, ind_tetra)? {
                    self.simpl_struct.bw_rem_tetra(ind_tetra);
                } else {
                    self.simpl_struct.bw_keep_tetra(ind_tetra)?;
                }
            } else {
                break;
            }
        }

        let nod = Node::Value(ind_vert);
        self.simpl_struct.bw_insert_node(nod)
    }

    fn insert_vertex_helper(&mut self, ind_vertex: usize, near_to: usize) -> Result<usize> {
        let now = Instant::now();
        let ind_tetrahedron = if let Ok(ind) = self.walk_by_visibility(ind_vertex, near_to) {
            ind
        } else {
            self.simpl_struct.clean_to_rem()?;
            self.walk_check_all(ind_vertex)?
        };
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        self.walk_ns = self.walk_ns + nano;

        let now = Instant::now();
        let added_tetra = self.insert_bw(ind_vertex, ind_tetrahedron)?;
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        self.insert_ns = self.insert_ns + nano;

        Ok(added_tetra[0])
    }

    fn insert_first_tetrahedron(&mut self, indices_to_insert: &mut Vec<usize>) -> Result<()> {
        let now = Instant::now();
        // first tetrahedron insertion
        if self.get_vertices().len() == indices_to_insert.len() {
            let ind1 = indices_to_insert.pop().unwrap();
            let ind2 = indices_to_insert.pop().unwrap();
            let pt1 = self.get_vertices()[ind1];
            let pt2 = self.get_vertices()[ind2];

            let mut aligned = Vec::new();
            let vec12 = [pt2[0] - pt1[0], pt2[1] - pt1[1], pt2[2] - pt1[2]];

            let i3 = indices_to_insert
                .iter()
                .rev()
                .enumerate()
                .map(|(e, &ind)| (e, self.get_vertices()[ind]))
                .map(|(e, pt)| (e, [pt[0] - pt1[0], pt[1] - pt1[1], pt[2] - pt1[2]]))
                .map(|(e, vec)| (e, vec[0] * vec12[0] + vec[1] * vec12[1] + vec[2] * vec12[2]))
                .map(|(e, scal)| if scal < 0.0 { (e, -scal) } else { (e, scal) })
                .max_by(|(_, val1), (_, val2)| val1.partial_cmp(val2).unwrap())
                .map(|(e, _)| e)
                .unwrap();

            let ind3 = indices_to_insert.remove(i3);
            let pt3 = self.get_vertices()[ind3];

            loop {
                if let Some(ind4) = indices_to_insert.pop() {
                    let pt4 = self.get_vertices()[ind4];

                    let sign = robust::orient3d(
                        Coord3D {
                            x: pt1[0],
                            y: pt1[1],
                            z: pt1[2],
                        },
                        Coord3D {
                            x: pt2[0],
                            y: pt2[1],
                            z: pt2[2],
                        },
                        Coord3D {
                            x: pt3[0],
                            y: pt3[1],
                            z: pt3[2],
                        },
                        Coord3D {
                            x: pt4[0],
                            y: pt4[1],
                            z: pt4[2],
                        },
                    );

                    if sign > 0. {
                        self.simpl_struct
                            .first_tetrahedron([ind1, ind2, ind3, ind4])?
                    } else if sign < 0. {
                        self.simpl_struct
                            .first_tetrahedron([ind1, ind3, ind2, ind4])?
                    } else {
                        aligned.push(ind4);
                        continue;
                    };
                } else {
                    return Err(anyhow::Error::msg("Could not find four non aligned points"));
                }

                break;
            }
            indices_to_insert.append(&mut aligned);
        }
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        log::info!("First tetrahedron computed in {}ms", nano as f32 / 1e6);
        Ok(())
    }

    /// insert a single vertex in the structure
    pub fn insert_vertex(&mut self, vertex: [f64; 3], near_to: Option<usize>) -> Result<()> {
        if self.simpl_struct.get_nb_tetrahedra() == 0 {
            return Err(anyhow::Error::msg(
                "Needs at least 1 tetrahedron to insert a single point",
            ));
        }
        let indices_to_insert = self.vertex_coordinates.len();
        self.vertex_coordinates.push(vertex);
        self.insert_vertex_helper(
            indices_to_insert,
            near_to.unwrap_or(self.simpl_struct.get_nb_tetrahedra() - 1),
        )?;
        self.simpl_struct.clean_to_rem()?;
        log::info!("Walks computed in {}ms", self.walk_ns as f32 / 1e6);
        log::info!("Insertions computed in {}ms", self.insert_ns as f32 / 1e6);
        Ok(())
    }

    /// Updates delaunay graph, including newly inserted vertices
    pub fn insert_vertices(
        &mut self,
        to_insert: &Vec<[f64; 3]>,
        reorder_points: bool,
    ) -> Result<()> {
        let mut indices_to_insert = Vec::new();
        for &vert in to_insert.iter() {
            indices_to_insert.push(self.vertex_coordinates.len());
            self.vertex_coordinates.push(vert);
        }

        if self.get_vertices().len() < 4 {
            return Err(anyhow::Error::msg(
                "Needs at least 4 vertices to compute Delaunay",
            ));
        }

        if reorder_points {
            let now = Instant::now();
            indices_to_insert = build_hilbert_curve_3d(self.get_vertices(), &indices_to_insert);
            let duration = now.elapsed();
            let nano = duration.as_nanos();
            log::info!("Hilbert curve computed in {}ms", nano as f32 / 1e6);
        }

        if self.simpl_struct.get_nb_tetrahedra() == 0 {
            self.insert_first_tetrahedron(&mut indices_to_insert)?;
        }

        let mut last_added = self.simpl_struct.get_nb_tetrahedra() - 1;
        loop {
            if let Some(ind_vertex) = indices_to_insert.pop() {
                last_added = self.insert_vertex_helper(ind_vertex, last_added)?;
            } else {
                break;
            }
        }
        self.simpl_struct.clean_to_rem()?;
        log::info!("Walks computed in {}ms", self.walk_ns as f32 / 1e6);
        log::info!("Insertions computed in {}ms", self.insert_ns as f32 / 1e6);

        Ok(())
    }

    /// Checks Delaunay graph validity (unit tests purpose)
    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        if !self.get_simplicial().is_valid()? {
            return Ok(false);
        }

        for ind_tetra in 0..self.get_simplicial().get_nb_tetrahedra() {
            if self.is_tetrahedron_flat(ind_tetra)? {
                log::warn!(
                    "Flat tetrahedron: {}",
                    self.get_simplicial()
                        .get_tetrahedron(ind_tetra)?
                        .to_string()
                );
                continue;
            }
            for ind_vert in 0..self.vertex_coordinates.len() {
                let in_sphere = self.is_vertex_strict_in_sphere(ind_vert, ind_tetra)?;
                if in_sphere {
                    log::error!(
                        "Non Delaunay tetrahedron: {}",
                        self.get_simplicial()
                            .get_tetrahedron(ind_tetra)?
                            .to_string()
                    );

                    valid = false;
                    break;
                }
            }
        }

        Ok(valid)
    }
}
