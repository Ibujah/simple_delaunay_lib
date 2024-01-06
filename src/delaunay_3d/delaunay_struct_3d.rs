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
    indices_to_insert: Vec<usize>,
    inserted_indices: Vec<usize>,
}

impl DelaunayStructure3D {
    /// Delaunay structure initialisation
    pub fn new() -> DelaunayStructure3D {
        DelaunayStructure3D {
            simpl_struct: SimplicialStructure3D::new(),
            vertex_coordinates: Vec::new(),
            indices_to_insert: Vec::new(),
            inserted_indices: Vec::new(),
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

    /// Adds a vertex to insert within the structure (need to call update_delaunay())
    pub fn add_vertex_to_insert(&mut self, to_insert: [f64; 3]) -> () {
        self.indices_to_insert.push(self.vertex_coordinates.len());
        self.vertex_coordinates.push(to_insert);
    }

    /// Adds a set of vertices to insert within the structure (need to call update_delaunay())
    pub fn add_vertices_to_insert(&mut self, to_insert: &Vec<[f64; 3]>) -> () {
        for &vert in to_insert.iter() {
            self.indices_to_insert.push(self.vertex_coordinates.len());
            self.vertex_coordinates.push(vert);
        }
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
        loop {
            if let Some(tri) = self.choose_tri(&vec_tri, &vert) {
                let tri_opp = tri.opposite();
                ind_tetra_cur = tri_opp.tetrahedron().ind();
                vec_tri.clear();
                let hes = tri_opp.halfedges();
                vec_tri.push(hes[(0 + side) % 3].neighbor().triangle());
                vec_tri.push(hes[(1 + side) % 3].neighbor().triangle());
                vec_tri.push(hes[(2 + side) % 3].neighbor().triangle());
                side = (side + 1) % 3;
            } else {
                return Ok(ind_tetra_cur);
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

    /// Updates delaunay graph, including newly inserted vertices
    pub fn update_delaunay(&mut self) -> Result<()> {
        if self.get_vertices().len() < 3 {
            return Err(anyhow::Error::msg(
                "Needs at least 3 vertices to compute Delaunay",
            ));
        }

        let now = Instant::now();
        self.indices_to_insert =
            build_hilbert_curve_3d(self.get_vertices(), &self.indices_to_insert);
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
            let vec12 = [pt2[0] - pt1[0], pt2[1] - pt1[1], pt2[2] - pt1[2]];

            let i3 = self
                .indices_to_insert
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

            let ind3 = self.indices_to_insert.remove(i3);
            let pt3 = self.get_vertices()[ind3];

            loop {
                if let Some(ind4) = self.indices_to_insert.pop() {
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
                        self.inserted_indices.push(ind1);
                        self.inserted_indices.push(ind2);
                        self.inserted_indices.push(ind3);
                        self.inserted_indices.push(ind4);
                        self.simpl_struct
                            .first_tetrahedron([ind1, ind2, ind3, ind4])?
                    } else if sign < 0. {
                        self.inserted_indices.push(ind1);
                        self.inserted_indices.push(ind2);
                        self.inserted_indices.push(ind3);
                        self.inserted_indices.push(ind4);
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
            self.indices_to_insert.append(&mut aligned);
        }
        self.is_valid()?;
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        log::info!("First tetrahedron computed in {}ms", nano as f32 / 1e6);

        let mut walk_ns = 0;
        let mut insert_ns = 0;
        let mut step = 0;
        let mut last_added = self.simpl_struct.get_nb_tetrahedra() - 1;
        loop {
            if let Some(ind_v) = self.indices_to_insert.pop() {
                let now = Instant::now();
                let ind_tetrahedron = self.walk_by_visibility(ind_v, last_added)?;

                let duration = now.elapsed();
                let nano = duration.as_nanos();
                walk_ns = walk_ns + nano;

                step = step + 1;

                let now = Instant::now();
                let added_tetra = self.insert_bw(ind_v, ind_tetrahedron)?;
                let duration = now.elapsed();
                let nano = duration.as_nanos();
                insert_ns = insert_ns + nano;
                last_added = added_tetra[0];

                self.inserted_indices.push(ind_v);
            } else {
                break;
            }
        }
        self.simpl_struct.clean_to_rem()?;
        log::info!("Walks computed in {}ms", walk_ns as f32 / 1e6);
        log::info!("Insertions computed in {}ms", insert_ns as f32 / 1e6);

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
            for &ind_vert in self.inserted_indices.iter() {
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
