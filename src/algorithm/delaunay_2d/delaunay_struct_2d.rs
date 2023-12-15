use anyhow::Result;
use nalgebra::base::*;
use robust::{self, Coord};
use std::collections::{HashMap, HashSet};
use std::time::Instant;

use svg::node::element;
use svg::node::element::path::Data;
use svg::Document;

use super::geometry_2d::{Circle, ExtendedCircle, Line};
use super::geometry_operations_2d::{
    build_hilbert_curve, circle_center_and_radius, line_normal_and_factor,
};
use super::simplicial_struct_2d::{self, Node, SimplicialStructure2D};

pub enum ExtendedTriangle {
    Triangle([[f64; 2]; 3]),
    Segment([[f64; 2]; 2]),
}

pub fn draw_triangle(
    document: Document,
    pts: &[Vector2<f64>; 3],
    color: &str,
    stroke_width: f64,
) -> Document {
    let [pt0, pt1, pt2] = pts;
    let data = Data::new()
        .move_to((pt0[0], pt0[1]))
        .line_by((pt1[0] - pt0[0], pt1[1] - pt0[1]))
        .line_by((pt2[0] - pt1[0], pt2[1] - pt1[1]))
        .close();

    let path = element::Path::new()
        .set("fill", "none")
        .set("stroke", color)
        .set("stroke-width", stroke_width)
        .set("d", data);

    document.add(path)
}

pub fn draw_circle(document: Document, ctr: &Vector2<f32>, rad: f32) -> Document {
    let circle = element::Circle::new()
        .set("cx", ctr[0])
        .set("cy", ctr[1])
        .set("r", rad)
        .set("stroke", "green")
        .set("stroke-width", 0.1)
        .set("fill", "none");

    document.add(circle)
}

/// 2D Delaunay structure
pub struct DelaunayStructure2D {
    simpl_struct: simplicial_struct_2d::SimplicialStructure2D,
    vertex_coordinates: Vec<[f64; 2]>,
    indices_to_insert: Vec<usize>,
    extended_triangles: HashMap<usize, ExtendedTriangle>,
}

impl DelaunayStructure2D {
    pub fn new() -> DelaunayStructure2D {
        DelaunayStructure2D {
            simpl_struct: simplicial_struct_2d::SimplicialStructure2D::new(),
            vertex_coordinates: Vec::new(),
            indices_to_insert: Vec::new(),
            extended_triangles: HashMap::new(),
        }
    }

    pub fn add_vertex_to_insert(&mut self, to_insert: Vector2<f64>) -> () {
        self.indices_to_insert.push(self.vertex_coordinates.len());
        self.vertex_coordinates.push([to_insert.x, to_insert.y]);
    }

    pub fn add_vertices_to_insert(&mut self, to_insert: &Vec<Vector2<f64>>) -> () {
        for &vert in to_insert.iter() {
            self.indices_to_insert.push(self.vertex_coordinates.len());
            self.vertex_coordinates.push([vert.x, vert.y]);
        }
    }

    pub fn get_simplicial(&self) -> &SimplicialStructure2D {
        &self.simpl_struct
    }

    pub fn get_vertices(&self) -> Vec<Vector2<f64>> {
        self.vertex_coordinates
            .iter()
            .map(|[x, y]| Vector2::new(*x, *y))
            .collect()
    }

    pub fn get_extended_triangles(&self) -> &HashMap<usize, ExtendedTriangle> {
        &self.extended_triangles
    }

    pub fn get_extended_circle(&self, ind: usize) -> Result<ExtendedCircle> {
        let ext_tri = self
            .extended_triangles
            .get(&ind)
            .ok_or(anyhow::Error::msg("Index not in triangles"))?;

        let res = match ext_tri {
            ExtendedTriangle::Triangle(tri) => {
                let pt1 = Vector2::new(tri[0][0] as f32, tri[0][1] as f32);
                let pt2 = Vector2::new(tri[1][0] as f32, tri[1][1] as f32);
                let pt3 = Vector2::new(tri[2][0] as f32, tri[2][1] as f32);
                let (ctr, rad) = circle_center_and_radius(&pt1, &pt2, &pt3)
                    .ok_or(anyhow::Error::msg("Could not compute circle"))?;

                ExtendedCircle::Circle(Circle::new(ctr, rad))
            }
            ExtendedTriangle::Segment(lin) => {
                let pt1 = Vector2::new(lin[0][0] as f32, lin[0][1] as f32);
                let pt2 = Vector2::new(lin[1][0] as f32, lin[1][1] as f32);

                let (nor, fac) = line_normal_and_factor(&pt1, &pt2);

                ExtendedCircle::Line(Line::new(nor, fac))
            }
        };

        Ok(res)
    }

    fn add_extended_triangle(&mut self, ind_triangle: usize) -> Result<()> {
        let triangle = self.simpl_struct.get_triangle(ind_triangle)?;
        let [he1, he2, he3] = triangle.halfedges();

        let node1 = he1.first_node();
        let node2 = he2.first_node();
        let node3 = he3.first_node();

        match (node1.ind(), node2.ind(), node3.ind()) {
            (
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt2 = self.vertex_coordinates[ind_v2];
                let pt3 = self.vertex_coordinates[ind_v3];
                if self
                    .extended_triangles
                    .insert(ind_triangle, ExtendedTriangle::Segment([pt2, pt3]))
                    .is_some()
                {
                    return Err(anyhow::Error::msg("Triangle index already exists"));
                }
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt3 = self.vertex_coordinates[ind_v3];
                if self
                    .extended_triangles
                    .insert(ind_triangle, ExtendedTriangle::Segment([pt3, pt1]))
                    .is_some()
                {
                    return Err(anyhow::Error::msg("Triangle index already exists"));
                }
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Infinity,
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt2 = self.vertex_coordinates[ind_v2];
                if self
                    .extended_triangles
                    .insert(ind_triangle, ExtendedTriangle::Segment([pt1, pt2]))
                    .is_some()
                {
                    return Err(anyhow::Error::msg("Triangle index already exists"));
                }
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt2 = self.vertex_coordinates[ind_v2];
                let pt3 = self.vertex_coordinates[ind_v3];
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
                // let sign = geometry_exact_computation_2d::ccw([pt1, pt2, pt3]);

                if sign > 0. {
                    if self
                        .extended_triangles
                        .insert(ind_triangle, ExtendedTriangle::Triangle([pt1, pt2, pt3]))
                        .is_some()
                    {
                        return Err(anyhow::Error::msg("Triangle index already exists"));
                    }
                } else if sign < 0. {
                    self.simpl_struct.get_triangle(ind_triangle)?.println();
                    return Err(anyhow::Error::msg("Misoriented triangle"));
                } else {
                    return Err(anyhow::Error::msg("Flat triangle"));
                }
            }
            (_, _, _) => {
                return Err(anyhow::Error::msg("Case should not happen"));
            }
        }

        Ok(())
    }

    fn rem_extended_triangle(&mut self, ind_triangle: usize) -> () {
        self.extended_triangles.remove(&ind_triangle);
    }

    fn is_vertex_in_circle(&self, ind_vert: usize, ind_tri: usize) -> bool {
        let vert = self.vertex_coordinates[ind_vert];
        let ext_tri = self.extended_triangles.get(&ind_tri).unwrap();

        let in_circle = match ext_tri {
            ExtendedTriangle::Triangle(tri) => {
                //let sign = geometry_exact_computation_2d::incircle(*tri, vert);
                let sign = robust::incircle(
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
                );
                sign >= 0.
            }
            ExtendedTriangle::Segment(lin) => {
                let sign = robust::orient2d(
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
                );
                // let sign = geometry_exact_computation_2d::ccw([lin[0], lin[1], vert]);

                if sign == 0. {
                    let [[x1, y1], [x2, y2]] = lin;
                    let [x, y] = vert;
                    let scal1 = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
                    let scal2 = (x1 - x2) * (x - x2) + (y1 - y2) * (y - y2);

                    scal1 > 0.0 && scal2 > 0.0
                } else {
                    sign > 0.
                }
            }
        };
        in_circle
    }

    fn choose_he<'a>(
        &self,
        vec_edg: &Vec<simplicial_struct_2d::IterHalfEdge<'a>>,
        vert: &[f64; 2],
    ) -> Option<simplicial_struct_2d::IterHalfEdge<'a>> {
        for &he in vec_edg {
            let ind1 = he.first_node().ind();
            let ind2 = he.last_node().ind();
            if let (Node::Value(v1), Node::Value(v2)) = (ind1, ind2) {
                let pt1 = self.vertex_coordinates[v1];
                let pt2 = self.vertex_coordinates[v2];
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
                if sign < 0. {
                    return Some(he);
                }

                // let sign = geometry_exact_computation_2d::ccw([pt1, pt2, *vert]);
                // if sign < 0 {
                //     return Some(he);
                // }
            }
        }
        None
    }

    fn walk_by_visibility(&self, ind_vert: usize, ind_starting_triangle: usize) -> Result<usize> {
        let vert = self.vertex_coordinates[ind_vert];
        let mut ind_tri_cur = ind_starting_triangle;
        let start_tri = self.simpl_struct.get_triangle(ind_tri_cur)?;
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

    fn search_path_surrounding_v0(
        &self,
        ind_vert: usize,
        ind_triangle: usize,
    ) -> Result<(Vec<usize>, HashSet<usize>, u128)> {
        if !self.is_vertex_in_circle(ind_vert, ind_triangle) {
            return Err(anyhow::Error::msg("Vertex not in circle"));
        }
        let [he1, he2, he3] = self.simpl_struct.get_triangle(ind_triangle)?.halfedges();
        let mut tri_to_rem = HashSet::new();
        tri_to_rem.insert(ind_triangle);
        let mut path = Vec::new();

        let mut he_to_analyse = vec![he1, he2, he3];

        let mut dur_in_circle = 0;
        loop {
            if let Some(he) = he_to_analyse.pop() {
                let he_opp = he.opposite_halfedge();
                let tri_opp = he_opp.face();
                let now = Instant::now();
                let in_circle = self.is_vertex_in_circle(ind_vert, tri_opp.ind());
                let duration = now.elapsed();
                let nano = duration.as_nanos();
                dur_in_circle = dur_in_circle + nano;
                if in_circle {
                    if !tri_to_rem.contains(&tri_opp.ind()) {
                        tri_to_rem.insert(tri_opp.ind());
                        he_to_analyse.push(he_opp.next_halfedge());
                        he_to_analyse.push(he_opp.prev_halfedge());
                    }
                } else {
                    // let ind1 = he.first_node().ind();
                    // let ind2 = he.last_node().ind();
                    // if let (Node::Value(v1), Node::Value(v2)) = (ind1, ind2) {
                    //     let vert1 = self.vertex_coordinates[v1];
                    //     let vert2 = self.vertex_coordinates[v2];
                    //     let vert = self.vertex_coordinates[ind_vert];
                    //     let sign = geometry_exact_computation_2d::ccw(&[vert1, vert2, vert]);

                    //     if sign < 0 {
                    //         return Err(anyhow::Error::msg("Misoriented triangle"));
                    //     } else if sign == 0 {
                    //         return Err(anyhow::Error::msg("Flat triangle"));
                    //     }
                    // }
                    path.push(he.ind());
                }
            } else {
                break;
            }
        }

        Ok((path, tri_to_rem, dur_in_circle))
    }

    fn search_path_surrounding(
        &self,
        ind_vert: usize,
        ind_triangle: usize,
    ) -> Result<(Vec<usize>, HashSet<usize>, u128)> {
        if !self.is_vertex_in_circle(ind_vert, ind_triangle) {
            return Err(anyhow::Error::msg("Vertex not in circle"));
        }
        let [he1, _, _] = self.simpl_struct.get_triangle(ind_triangle)?.halfedges();
        let mut tri_to_rem = HashSet::new();
        let mut path = Vec::new();

        let mut he_cur = he1;
        let first_node = he1.first_node().ind();

        let mut dur_in_circle = 0;
        loop {
            let tri_cur = he_cur.face().ind();
            tri_to_rem.insert(tri_cur);
            let he_opp = he_cur.opposite_halfedge();
            let tri_opp = he_opp.face();

            let now = Instant::now();
            let in_circle = if tri_to_rem.contains(&tri_opp.ind()) {
                true
            } else {
                self.is_vertex_in_circle(ind_vert, tri_opp.ind())
            };
            let duration = now.elapsed();
            let nano = duration.as_nanos();
            dur_in_circle = dur_in_circle + nano;
            if in_circle {
                he_cur = he_opp.next_halfedge();
            } else {
                path.push(he_cur.ind());
                if he_cur.last_node().ind().equals(&first_node) {
                    break;
                }
                he_cur = he_cur.next_halfedge();
            }
        }

        Ok((path, tri_to_rem, dur_in_circle))
    }

    pub fn update_delaunay(&mut self) -> Result<()> {
        if self.vertex_coordinates.len() < 3 {
            return Err(anyhow::Error::msg(
                "Needs at least 3 vertices to compute Delaunay",
            ));
        }

        let now = Instant::now();
        self.indices_to_insert = build_hilbert_curve(&self.get_vertices(), &self.indices_to_insert);
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        println!("Hilbert curve computed in {}ms", nano as f32 / 1e6);

        let now = Instant::now();
        // first triangle insertion
        let mut last_added_triangle =
            if self.vertex_coordinates.len() == self.indices_to_insert.len() {
                let ind1 = self.indices_to_insert.pop().unwrap();
                let ind2 = self.indices_to_insert.pop().unwrap();

                let pt1 = self.vertex_coordinates[ind1];
                let pt2 = self.vertex_coordinates[ind2];

                let mut aligned = Vec::new();

                let ind_tri1 = loop {
                    let ind3 = self.indices_to_insert.pop().unwrap();
                    let pt3 = self.vertex_coordinates[ind3];

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
                    // let sign = geometry_exact_computation_2d::ccw([pt1, pt2, pt3]);

                    let [ind_tri1, ind_tri2, ind_tri3, ind_tri4] = if sign > 0. {
                        self.simpl_struct.first_triangle([ind1, ind2, ind3])?
                    } else if sign < 0. {
                        self.simpl_struct.first_triangle([ind1, ind3, ind2])?
                    } else {
                        aligned.push(ind3);
                        continue;
                    };
                    self.add_extended_triangle(ind_tri1)?;
                    self.add_extended_triangle(ind_tri2)?;
                    self.add_extended_triangle(ind_tri3)?;
                    self.add_extended_triangle(ind_tri4)?;
                    break ind_tri1;
                };

                self.indices_to_insert.append(&mut aligned);

                if !self.simpl_struct.is_valid()? {
                    return Err(anyhow::Error::msg(
                        "Simplicial structure not valid anymore (first step)",
                    ));
                }

                ind_tri1
            } else {
                *self.extended_triangles.iter().next().unwrap().0
            };
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        println!("First triangle computed in {}ms", nano as f32 / 1e6);

        let mut walk_ms = 0;
        let mut surr_path_ms = 0;
        let mut in_circle_ms = 0;
        let mut insert_ms = 0;
        let mut update_ms = 0;
        let mut step = 0;
        loop {
            if let Some(ind_v) = self.indices_to_insert.pop() {
                let now = Instant::now();
                // let ind_triangle = self.search_circle_containing(ind_v)?;
                let ind_triangle = self.walk_by_visibility(ind_v, last_added_triangle)?;
                let duration = now.elapsed();
                let milli = duration.as_nanos();
                walk_ms = walk_ms + milli;

                let now = Instant::now();
                // println!("v1");
                let (path, tri_to_rem, dur_in_circle) =
                    self.search_path_surrounding_v0(ind_v, ind_triangle)?;
                // println!("v2");
                // let (path, tri_to_rem, dur_in_circle) =
                //     self.search_path_surrounding(ind_v, ind_triangle)?;
                let duration = now.elapsed();
                let milli = duration.as_nanos();
                surr_path_ms = surr_path_ms + milli;
                in_circle_ms = in_circle_ms + dur_in_circle;

                step = step + 1;

                let now = Instant::now();
                let tri_to_add =
                    self.simpl_struct
                        .insert_node_within_path(ind_v, &path, &tri_to_rem)?;
                let duration = now.elapsed();
                let milli = duration.as_nanos();
                insert_ms = insert_ms + milli;

                let now = Instant::now();
                //update triangles
                for &ind_tri in tri_to_rem.iter() {
                    self.rem_extended_triangle(ind_tri);
                }
                for &ind_tri in tri_to_add.iter() {
                    self.add_extended_triangle(ind_tri)?;
                }
                let duration = now.elapsed();
                let milli = duration.as_nanos();
                update_ms = update_ms + milli;
                for tri_added in tri_to_add.iter() {
                    let tri = self.extended_triangles.get(tri_added).unwrap();
                    match *tri {
                        ExtendedTriangle::Triangle(_) => last_added_triangle = *tri_added,
                        ExtendedTriangle::Segment(_) => (),
                    }
                }
            } else {
                break;
            }
        }
        println!("Walks computed in {}ms", walk_ms as f32 / 1e6);
        println!(
            "Surrounding paths computed in {}ms ({}ms for incircle computation)",
            surr_path_ms as f32 / 1e6,
            in_circle_ms as f32 / 1e6
        );
        println!("Insertions computed in {}ms", insert_ms as f32 / 1e6);
        println!("Circle updates computed in {}ms", update_ms as f32 / 1e6);

        // if !self.simpl_struct.is_valid()? {
        //     return Err(anyhow::Error::msg("Simplicial structure not valid anymore"));
        // }

        Ok(())
    }

    pub fn draw_svg(&self, name: String, highlight: &HashSet<usize>) -> Result<()> {
        let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));

        let rect = element::Rectangle::new()
            .set("x", -50)
            .set("y", -50)
            .set("width", 1100)
            .set("height", 1100)
            .set("fill", "white");
        document = document.add(rect);

        for ind_triangle in self.get_simplicial().get_triangle_indices() {
            let tri = self.get_simplicial().get_triangle(ind_triangle)?;

            let [h1, h2, h3] = tri.halfedges();
            let ind_pt1 = h1.first_node().ind();
            let ind_pt2 = h2.first_node().ind();
            let ind_pt3 = h3.first_node().ind();

            if let (Node::Value(val1), Node::Value(val2), Node::Value(val3)) =
                (ind_pt1, ind_pt2, ind_pt3)
            {
                let pt1 = self.get_vertices()[val1];
                let pt2 = self.get_vertices()[val2];
                let pt3 = self.get_vertices()[val3];
                if highlight.contains(&ind_triangle) {
                    document = draw_triangle(
                        document,
                        &[pt1 * 1000., pt2 * 1000., pt3 * 1000.],
                        "red",
                        2.0,
                    );
                } else {
                    document = draw_triangle(
                        document,
                        &[pt1 * 1000., pt2 * 1000., pt3 * 1000.],
                        "black",
                        1.0,
                    );
                }
            }
        }

        for ind_triangle in self.get_simplicial().get_triangle_indices() {
            let extended_cir = self.get_extended_circle(ind_triangle)?;
            if let ExtendedCircle::Circle(circle) = extended_cir {
                document = draw_circle(document, &(circle.center * 1000.), circle.radius * 1000.);
            }
        }

        svg::save(name, &document).unwrap();
        Ok(())
    }
}
