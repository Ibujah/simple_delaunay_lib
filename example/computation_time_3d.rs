use anyhow::Result;
use env_logger;
use rand::Rng;
use std::time::Instant;

use delaunay_lib::delaunay_3d::delaunay_struct_3d::DelaunayStructure3D;

fn generate_random_vertices(nb_vert: usize) -> Vec<[f64; 3]> {
    let mut rng = rand::thread_rng();

    let mut vec_pts: Vec<[f64; 3]> = Vec::new();
    for _ in 0..nb_vert {
        let (x, y, z): (f64, f64, f64) = rng.gen();
        vec_pts.push([x, y, z]);
    }
    vec_pts
}

fn main() -> Result<()> {
    env_logger::init();

    let nb_vert_array = [10, 100, 1000, 10000, 100000, 1000000];

    for nb_vert in nb_vert_array {
        let vec_pts = generate_random_vertices(nb_vert);

        let now = Instant::now();
        let mut del_struct = DelaunayStructure3D::new();
        del_struct.insert_vertices(&vec_pts, true)?;
        let duration = now.elapsed();
        let milli = duration.as_millis();

        println!("{} vertices: {}ms", nb_vert, milli);
    }

    Ok(())
}
