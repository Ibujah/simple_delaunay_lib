use anyhow::Result;
use nalgebra::base::*;
use robust::{self, Coord};
use std::collections::HashSet;
use std::time::Instant;

use svg::node::element;
use svg::node::element::path::Data;
use svg::Document;

use super::geometry_2d::{Circle, ExtendedCircle, Line};
use super::geometry_operations_2d::{
    build_hilbert_curve, circle_center_and_radius, is_convex, line_normal_and_factor,
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

    pub fn get_extended_triangle(&self, ind_triangle: usize) -> Result<ExtendedTriangle> {
        let triangle = self.simpl_struct.get_triangle(ind_triangle)?;
        let [he1, he2, he3] = triangle.halfedges();

        let node1 = he1.first_node();
        let node2 = he2.first_node();
        let node3 = he3.first_node();

        let ext_tri = match (node1, node2, node3) {
            (
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt2 = self.vertex_coordinates[ind_v2];
                let pt3 = self.vertex_coordinates[ind_v3];
                ExtendedTriangle::Segment([pt2, pt3])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Infinity,
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt3 = self.vertex_coordinates[ind_v3];
                ExtendedTriangle::Segment([pt3, pt1])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Infinity,
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt2 = self.vertex_coordinates[ind_v2];
                ExtendedTriangle::Segment([pt1, pt2])
            }
            (
                simplicial_struct_2d::Node::Value(ind_v1),
                simplicial_struct_2d::Node::Value(ind_v2),
                simplicial_struct_2d::Node::Value(ind_v3),
            ) => {
                let pt1 = self.vertex_coordinates[ind_v1];
                let pt2 = self.vertex_coordinates[ind_v2];
                let pt3 = self.vertex_coordinates[ind_v3];
                ExtendedTriangle::Triangle([pt1, pt2, pt3])
            }
            (_, _, _) => {
                return Err(anyhow::Error::msg("Case should not happen"));
            }
        };

        Ok(ext_tri)
    }

    pub fn get_extended_circle(&self, ind: usize) -> Result<ExtendedCircle> {
        let ext_tri = self.get_extended_triangle(ind)?;

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

    fn is_vertex_in_circle(&self, ind_vert: usize, ind_tri: usize) -> Result<i8> {
        let vert = self.vertex_coordinates[ind_vert];
        let ext_tri = self.get_extended_triangle(ind_tri)?;

        let in_circle = match ext_tri {
            ExtendedTriangle::Triangle(tri) => {
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
                if sign > 0. {
                    1
                } else if sign < 0. {
                    -1
                } else {
                    0
                }
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

                if sign > 0. {
                    1
                } else if sign < 0. {
                    -1
                } else {
                    let [[x1, y1], [x2, y2]] = lin;
                    let [x, y] = vert;
                    let scal1 = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
                    let scal2 = (x1 - x2) * (x - x2) + (y1 - y2) * (y - y2);

                    if scal1 > 0.0 && scal2 > 0.0 {
                        0
                    } else {
                        -1
                    }
                }
            }
        };
        Ok(in_circle)
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

    fn should_flip_halfedge(&self, ind_he: usize) -> Result<bool> {
        let he = self.simpl_struct.get_halfedge(ind_he)?;
        let ind_tri_abd = he.face().ind();
        let node_a = he.prev_halfedge().first_node();
        let node_b = he.first_node();

        let ind_tri_bcd = he.opposite_halfedge().face().ind();
        let node_c = he.opposite_halfedge().prev_halfedge().first_node();
        let node_d = he.opposite_halfedge().first_node();

        match (node_a, node_b, node_c, node_d) {
            (Node::Value(ind_node_a), Node::Value(_), Node::Value(ind_node_c), Node::Value(_)) => {
                Ok(self.is_vertex_in_circle(ind_node_c, ind_tri_abd)? == 1
                    || self.is_vertex_in_circle(ind_node_a, ind_tri_bcd)? == 1)
            }
            (Node::Infinity, Node::Value(_), Node::Value(ind_node_c), Node::Value(_)) => {
                Ok(self.is_vertex_in_circle(ind_node_c, ind_tri_abd)? == 1
                    || self.is_triangle_flat(ind_tri_bcd)?)
            }
            (
                Node::Value(ind_node_a),
                Node::Infinity,
                Node::Value(ind_node_c),
                Node::Value(ind_node_d),
            ) => {
                let pt_a = self.vertex_coordinates[ind_node_a];
                let pt_c = self.vertex_coordinates[ind_node_c];
                let pt_d = self.vertex_coordinates[ind_node_d];
                Ok(is_convex(pt_c, pt_d, pt_a) == 1)
            }
            (Node::Value(ind_node_a), Node::Value(_), Node::Infinity, Node::Value(_)) => {
                Ok(self.is_triangle_flat(ind_tri_abd)?
                    || self.is_vertex_in_circle(ind_node_a, ind_tri_bcd)? == 1)
            }
            (
                Node::Value(ind_node_a),
                Node::Value(ind_node_b),
                Node::Value(ind_node_c),
                Node::Infinity,
            ) => {
                let pt_a = self.vertex_coordinates[ind_node_a];
                let pt_b = self.vertex_coordinates[ind_node_b];
                let pt_c = self.vertex_coordinates[ind_node_c];
                Ok(is_convex(pt_a, pt_b, pt_c) == 1)
            }
            (_, _, _, _) => Err(anyhow::Error::msg("Multiple infinity linked together")),
        }
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
        if self.vertex_coordinates.len() == self.indices_to_insert.len() {
            let ind1 = self.indices_to_insert.pop().unwrap();
            let ind2 = self.indices_to_insert.pop().unwrap();
            let pt1 = self.vertex_coordinates[ind1];
            let pt2 = self.vertex_coordinates[ind2];

            let mut aligned = Vec::new();

            loop {
                if let Some(ind3) = self.indices_to_insert.pop() {
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
                        "Simplicial structure not valid anymore (first step)",
                    ));
                }

                break;
            }

            self.indices_to_insert.append(&mut aligned);

            if !self.simpl_struct.is_valid()? {
                return Err(anyhow::Error::msg(
                    "Simplicial structure not valid anymore (first step)",
                ));
            }
        }
        let duration = now.elapsed();
        let nano = duration.as_nanos();
        println!("First triangle computed in {}ms", nano as f32 / 1e6);

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
                        let he = self.simpl_struct.get_halfedge(ind_he)?;
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
        println!("Walks computed in {}ms", walk_ms as f32 / 1e6);
        println!("Insertions computed in {}ms", insert_ms as f32 / 1e6);
        println!("Flips computed in {}ms", flip_ms as f32 / 1e6);

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

        for ind_triangle in 0..self.get_simplicial().get_nb_triangles() {
            let tri = self.get_simplicial().get_triangle(ind_triangle)?;

            let [h1, h2, h3] = tri.halfedges();
            let ind_pt1 = h1.first_node();
            let ind_pt2 = h2.first_node();
            let ind_pt3 = h3.first_node();

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

        for ind_triangle in 0..self.get_simplicial().get_nb_triangles() {
            if let Ok(extended_cir) = self.get_extended_circle(ind_triangle) {
                if let ExtendedCircle::Circle(circle) = extended_cir {
                    document =
                        draw_circle(document, &(circle.center * 1000.), circle.radius * 1000.);
                }
            }
        }

        svg::save(name, &document).unwrap();
        Ok(())
    }

    pub fn is_valid(&self) -> Result<bool> {
        let mut valid = true;

        for ind_tri in 0..self.simpl_struct.get_nb_triangles() {
            if self.is_triangle_flat(ind_tri)? {
                print!("not valid: ");
                self.simpl_struct.get_triangle(ind_tri)?.println();
                valid = false;
            }
            for &ind_vert in self.inserted_indices.iter() {
                let in_circle = self.is_vertex_in_circle(ind_vert, ind_tri)? == 1;
                if in_circle {
                    print!("not valid: ");
                    self.simpl_struct.get_triangle(ind_tri)?.println();
                    valid = false;
                }
            }
        }

        Ok(valid)
    }
}
