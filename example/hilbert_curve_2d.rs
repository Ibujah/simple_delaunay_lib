use anyhow::Result;
use env_logger;
use rand::Rng;

use svg::node::element;
use svg::node::element::path::Data;
use svg::Document;

use simple_delaunay_lib::delaunay_2d::geometry_operations_2d::build_hilbert_curve;

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
    Ok(())
}
