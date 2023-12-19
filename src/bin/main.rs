use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::time::Instant;

use svg::node::element::path::Data;
use svg::node::element::Circle;
use svg::node::element::Path;
use svg::node::element::Rectangle;
use svg::Document;

use delaunay_lib::delaunay::delaunay_2d::delaunay_struct_2d;
use delaunay_lib::delaunay::delaunay_2d::geometry_operations_2d::build_hilbert_curve;
use delaunay_lib::delaunay::delaunay_2d::simplicial_struct_2d::Node;

pub fn draw_triangle(document: Document, pts: &[Vector2<f64>; 3]) -> Document {
    let [pt0, pt1, pt2] = pts;
    let data = Data::new()
        .move_to((pt0[0], pt0[1]))
        .line_by((pt1[0] - pt0[0], pt1[1] - pt0[1]))
        .line_by((pt2[0] - pt1[0], pt2[1] - pt1[1]))
        .close();

    let path = Path::new()
        .set("fill", "none")
        .set("stroke", "black")
        .set("stroke-width", 1.0)
        .set("d", data);

    document.add(path)
}

pub fn draw_circle(document: Document, ctr: &Vector2<f64>, rad: f64) -> Document {
    let circle = Circle::new()
        .set("cx", ctr[0])
        .set("cy", ctr[1])
        .set("r", rad)
        .set("stroke", "green")
        .set("stroke-width", 0.1)
        .set("fill", "none");

    document.add(circle)
}

fn main() -> Result<()> {
    let mut rng = rand::thread_rng();

    let mut vec_pts: Vec<Vector2<f64>> = Vec::new();
    let mut vec_inds: Vec<usize> = Vec::new();
    for ind in 0..10000 {
        let (x, y): (f64, f64) = rng.gen();
        vec_pts.push(Vector2::new(x, y));
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
    let mut del_struct = delaunay_struct_2d::DelaunayStructure2D::new();
    del_struct.add_vertices_to_insert(&vec_pts);
    del_struct.update_delaunay()?;
    let duration = now.elapsed();
    let milli = duration.as_millis();
    println!("Delaunay computed in {}ms", milli);

    println!("Checking delaunay");
    if del_struct.is_valid()? {
        println!("Delaunay is valid");
    } else {
        println!("Non valid delaunay!");
    }

    if del_struct.get_vertices().len() <= 10000 {
        let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));

        let rect = Rectangle::new()
            .set("x", -50)
            .set("y", -50)
            .set("width", 1100)
            .set("height", 1100)
            .set("fill", "white");
        document = document.add(rect);

        for ind_triangle in 0..del_struct.get_simplicial().get_nb_triangles() {
            let tri = del_struct.get_simplicial().get_triangle(ind_triangle)?;

            let [h1, h2, h3] = tri.halfedges();
            let ind_pt1 = h1.first_node();
            let ind_pt2 = h2.first_node();
            let ind_pt3 = h3.first_node();

            if let (Node::Value(val1), Node::Value(val2), Node::Value(val3)) =
                (ind_pt1, ind_pt2, ind_pt3)
            {
                let pt1 = del_struct.get_vertices()[val1];
                let pt2 = del_struct.get_vertices()[val2];
                let pt3 = del_struct.get_vertices()[val3];
                document = draw_triangle(document, &[pt1 * 1000., pt2 * 1000., pt3 * 1000.]);
            }
        }

        // for ind_triangle in del_struct.get_simplicial().get_triangle_indices() {
        //     let extended_cir = del_struct.get_extended_circle(ind_triangle)?;
        //     if let geometry_2d::ExtendedCircle::Circle(circle) = extended_cir {
        //         document = draw_circle(document, &(circle.center * 1000.), circle.radius * 1000.);
        //     }
        // }

        svg::save("delaunay.svg", &document).unwrap();

        let vec_inds = build_hilbert_curve(&vec_pts, &vec_inds);

        let mut document = Document::new().set("viewBox", (-50, -50, 1100, 1100));

        let rect = Rectangle::new()
            .set("x", -50)
            .set("y", -50)
            .set("width", 1100)
            .set("height", 1100)
            .set("fill", "white");
        document = document.add(rect);

        for ind in 0..(vec_pts.len() - 1) {
            let pt0 = vec_pts[vec_inds[ind]] * 1000.;
            let pt1 = vec_pts[vec_inds[ind + 1]] * 1000.;
            let data = Data::new()
                .move_to((pt0.x, pt0.y))
                .line_by((pt1.x - pt0.x, pt1.y - pt0.y));

            let path = Path::new()
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
