use anyhow::Result;
use env_logger;
use log;
use nalgebra::base::*;
use rand::Rng;
use std::collections::HashSet;
use std::time::Instant;

use svg::node::element;
use svg::node::element::path::Data;
use svg::Document;

use delaunay_lib::delaunay::delaunay_2d::delaunay_struct_2d::{
    DelaunayStructure2D, ExtendedTriangle,
};
use delaunay_lib::delaunay::delaunay_2d::geometry_operations_2d::{
    build_hilbert_curve, circle_center_and_radius, line_normal_and_factor,
};
use delaunay_lib::delaunay::delaunay_2d::simplicial_struct_2d::Node;

#[derive(Copy, Clone)]
pub struct Circle {
    pub center: Vector2<f64>,
    pub radius: f64,
}

#[derive(Copy, Clone)]
pub struct Line {
    pub normal: Vector2<f64>,
    pub factor: f64,
}

pub enum ExtendedCircle {
    Circle(Circle),
    Line(Line),
}

impl ExtendedCircle {
    pub fn is_vertex_in(&self, vert: &Vector2<f64>) -> bool {
        match self {
            ExtendedCircle::Circle(circle) => circle.is_vertex_in(vert),
            ExtendedCircle::Line(line) => line.is_vertex_in(vert),
        }
    }
}

impl Circle {
    pub fn new(center: Vector2<f64>, radius: f64) -> Circle {
        Circle { center, radius }
    }
    pub fn is_vertex_in(&self, vert: &Vector2<f64>) -> bool {
        (self.center - vert).norm() - self.radius <= 0.
    }
}

impl Line {
    pub fn new(normal: Vector2<f64>, factor: f64) -> Line {
        Line { normal, factor }
    }
    pub fn is_vertex_in(&self, vert: &Vector2<f64>) -> bool {
        self.normal.dot(&vert) - self.factor <= 0.
    }
}

pub fn get_extended_circle(
    delaunay_struct_2d: &DelaunayStructure2D,
    ind: usize,
) -> Result<ExtendedCircle> {
    let ext_tri = delaunay_struct_2d.get_extended_triangle(ind)?;

    let res = match ext_tri {
        ExtendedTriangle::Triangle(tri) => {
            let pt1 = Vector2::new(tri[0][0], tri[0][1]);
            let pt2 = Vector2::new(tri[1][0], tri[1][1]);
            let pt3 = Vector2::new(tri[2][0], tri[2][1]);
            let (ctr, rad) = circle_center_and_radius(&pt1, &pt2, &pt3)
                .ok_or(anyhow::Error::msg("Could not compute circle"))?;

            ExtendedCircle::Circle(Circle::new(ctr, rad))
        }
        ExtendedTriangle::Segment(lin) => {
            let pt1 = Vector2::new(lin[0][0], lin[0][1]);
            let pt2 = Vector2::new(lin[1][0], lin[1][1]);

            let (nor, fac) = line_normal_and_factor(&pt1, &pt2);

            ExtendedCircle::Line(Line::new(nor, fac))
        }
    };

    Ok(res)
}

pub fn draw_triangle(
    document: Document,
    pts: &[[f64; 2]; 3],
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

pub fn draw_circle(document: Document, ctr: &Vector2<f64>, rad: f64) -> Document {
    let circle = element::Circle::new()
        .set("cx", ctr[0])
        .set("cy", ctr[1])
        .set("r", rad)
        .set("stroke", "green")
        .set("stroke-width", 0.1)
        .set("fill", "none");

    document.add(circle)
}

pub fn draw_svg(
    delaunay: &DelaunayStructure2D,
    name: String,
    highlight: Option<HashSet<usize>>,
    draw_circles: bool,
) -> Result<()> {
    let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));
    let highlight = highlight.unwrap_or(HashSet::new());

    let rect = element::Rectangle::new()
        .set("x", -50)
        .set("y", -50)
        .set("width", 1100)
        .set("height", 1100)
        .set("fill", "white");
    document = document.add(rect);

    for ind_triangle in 0..delaunay.get_simplicial().get_nb_triangles() {
        let tri = delaunay.get_simplicial().get_triangle(ind_triangle)?;

        let [h1, h2, h3] = tri.halfedges();
        let ind_pt1 = h1.first_node();
        let ind_pt2 = h2.first_node();
        let ind_pt3 = h3.first_node();

        if let (Node::Value(val1), Node::Value(val2), Node::Value(val3)) =
            (ind_pt1, ind_pt2, ind_pt3)
        {
            let pt1 = delaunay.get_vertices()[val1];
            let pt2 = delaunay.get_vertices()[val2];
            let pt3 = delaunay.get_vertices()[val3];
            if highlight.contains(&ind_triangle) {
                document = draw_triangle(
                    document,
                    &[
                        [pt1[0] * 1000., pt1[1] * 1000.],
                        [pt2[0] * 1000., pt2[1] * 1000.],
                        [pt3[0] * 1000., pt3[1] * 1000.],
                    ],
                    "red",
                    2.0,
                );
            } else {
                document = draw_triangle(
                    document,
                    &[
                        [pt1[0] * 1000., pt1[1] * 1000.],
                        [pt2[0] * 1000., pt2[1] * 1000.],
                        [pt3[0] * 1000., pt3[1] * 1000.],
                    ],
                    "black",
                    1.0,
                );
            }
        }
    }

    if draw_circles {
        for ind_triangle in 0..delaunay.get_simplicial().get_nb_triangles() {
            if let Ok(extended_cir) = get_extended_circle(delaunay, ind_triangle) {
                if let ExtendedCircle::Circle(circle) = extended_cir {
                    document =
                        draw_circle(document, &(circle.center * 1000.), circle.radius * 1000.);
                }
            }
        }
    }

    svg::save(name, &document).unwrap();
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let mut rng = rand::thread_rng();

    let mut vec_pts: Vec<[f64; 2]> = Vec::new();
    let mut vec_inds: Vec<usize> = Vec::new();
    for ind in 0..10000 {
        let (x, y): (f64, f64) = rng.gen();
        vec_pts.push([x, y]);
        vec_inds.push(ind);
    }
    // for ind in 0..10000 {
    //     let ind1 = ind % 100;
    //     let ind2 = ind / 100;

    //     let x = (ind1 as f64) / 100.;
    //     let y = (ind2 as f64) / 100.;
    //     vec_pts.push(Vector2::new(x, y));
    //     vec_inds.push(ind);
    // }

    let now = Instant::now();
    let mut del_struct = DelaunayStructure2D::new();
    del_struct.add_vertices_to_insert(&vec_pts);
    del_struct.update_delaunay()?;
    let duration = now.elapsed();
    let milli = duration.as_millis();
    log::info!("Delaunay computed in {}ms", milli);

    log::info!("Checking delaunay");
    if del_struct.is_valid()? {
        log::info!("Delaunay is valid");
    } else {
        log::info!("Non valid delaunay!");
    }

    if del_struct.get_vertices().len() <= 10000 {
        draw_svg(&del_struct, "delaunay.svg".to_string(), None, false)?;

        let vec_inds = build_hilbert_curve(&vec_pts, &vec_inds);

        let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));

        let rect = element::Rectangle::new()
            .set("x", -50)
            .set("y", -50)
            .set("width", 1100)
            .set("height", 1100)
            .set("fill", "white");
        document = document.add(rect);

        for ind in 0..(vec_pts.len() - 1) {
            let pt0_x = vec_pts[vec_inds[ind]][0] * 1000.;
            let pt0_y = vec_pts[vec_inds[ind]][1] * 1000.;
            let pt1_x = vec_pts[vec_inds[ind + 1]][0] * 1000.;
            let pt1_y = vec_pts[vec_inds[ind + 1]][1] * 1000.;
            let data = Data::new()
                .move_to((pt0_x, pt0_y))
                .line_by((pt1_x - pt0_x, pt1_y - pt0_y));

            let path = element::Path::new()
                .set("fill", "none")
                .set("stroke", "black")
                .set("stroke-width", 1.0)
                .set("d", data);

            document = document.add(path)
        }

        svg::save("hilbert_path.svg", &document).unwrap();
    }
    Ok(())
}
