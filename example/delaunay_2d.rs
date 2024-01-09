use anyhow::Result;
use env_logger;
use log;
use nalgebra::base::*;
use rand::Rng;
use std::time::Instant;

use svg::node::element;
use svg::node::element::path::Data;
use svg::Document;

use delaunay_lib::delaunay_2d::delaunay_struct_2d::{DelaunayStructure2D, ExtendedTriangle};
use delaunay_lib::delaunay_2d::simplicial_struct_2d::Node;

fn circle_center_and_radius(
    pt1: &Vector2<f64>,
    pt2: &Vector2<f64>,
    pt3: &Vector2<f64>,
) -> Option<(Vector2<f64>, f64)> {
    // I'm sure there is a justified equation behind those lines, but I don't remember it...

    let mat = Matrix3x2::new(
        pt2[0] - pt1[0],
        pt2[1] - pt1[1],
        pt3[0] - pt2[0],
        pt3[1] - pt2[1],
        pt1[0] - pt3[0],
        pt1[1] - pt3[1],
    );

    let b = Vector3::new(
        0.5 * (pt2 - pt1).norm_squared() + (pt2 - pt1).dot(pt1),
        0.5 * (pt3 - pt2).norm_squared() + (pt3 - pt2).dot(pt2),
        0.5 * (pt1 - pt3).norm_squared() + (pt1 - pt3).dot(pt3),
    );

    let mat_mod = mat.transpose() * mat;
    let b_mod = mat.transpose() * b;

    let opt_center = mat_mod.lu().solve(&b_mod);
    if let Some(center) = opt_center {
        let radius = (center - pt1).norm();
        Some((center, radius))
    } else {
        None
    }
}

fn get_if_circle(
    delaunay_struct_2d: &DelaunayStructure2D,
    ind: usize,
) -> Result<Option<(Vector2<f64>, f64)>> {
    let ext_tri = delaunay_struct_2d.get_extended_triangle(ind)?;

    if let ExtendedTriangle::Triangle(tri) = ext_tri {
        let pt1 = Vector2::new(tri[0][0], tri[0][1]);
        let pt2 = Vector2::new(tri[1][0], tri[1][1]);
        let pt3 = Vector2::new(tri[2][0], tri[2][1]);
        let (ctr, rad) = circle_center_and_radius(&pt1, &pt2, &pt3)
            .ok_or(anyhow::Error::msg("Could not compute circle"))?;

        Ok(Some((ctr, rad)))
    } else {
        Ok(None)
    }
}

fn draw_triangle(
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

fn draw_circle(document: Document, ctr: &Vector2<f64>, rad: f64) -> Document {
    let circle = element::Circle::new()
        .set("cx", ctr[0])
        .set("cy", ctr[1])
        .set("r", rad)
        .set("stroke", "green")
        .set("stroke-width", 0.1)
        .set("fill", "none");

    document.add(circle)
}

fn draw_svg(delaunay: &DelaunayStructure2D, name: String, draw_circles: bool) -> Result<()> {
    let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));

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

    if draw_circles {
        for ind_triangle in 0..delaunay.get_simplicial().get_nb_triangles() {
            if let Some((center, radius)) = get_if_circle(delaunay, ind_triangle)? {
                document = draw_circle(document, &(center * 1000.), radius * 1000.);
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
    for ind in 0..1000 {
        let (x, y): (f64, f64) = rng.gen();
        vec_pts.push([x, y]);
        vec_inds.push(ind);
    }

    let now = Instant::now();
    let mut del_struct = DelaunayStructure2D::new();
    del_struct.insert_vertices(&vec_pts, true)?;
    let duration = now.elapsed();
    let milli = duration.as_millis();
    log::info!("Delaunay computed in {}ms", milli);

    log::info!("Checking delaunay");
    if del_struct.is_valid()? {
        log::info!("Delaunay is valid");
    } else {
        log::info!("Non valid delaunay!");
    }

    draw_svg(&del_struct, "delaunay_2d.svg".to_string(), false)?;
    Ok(())
}
