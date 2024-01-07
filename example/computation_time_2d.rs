use anyhow::Result;
use env_logger;
use rand::Rng;
use std::time::Instant;

use delaunay_lib::delaunay_2d::delaunay_struct_2d::DelaunayStructure2D;

fn generate_random_vertices(nb_vert: usize) -> Vec<[f64; 2]> {
    let mut rng = rand::thread_rng();

    let mut vec_pts: Vec<[f64; 2]> = Vec::new();
    for _ in 0..nb_vert {
        let (x, y): (f64, f64) = rng.gen();
        vec_pts.push([x, y]);
    }
    vec_pts
}

fn main() -> Result<()> {
    env_logger::init();

    let nb_vert_array = [10, 100, 1000, 10000, 100000, 1000000];

    for nb_vert in nb_vert_array {
        let vec_pts = generate_random_vertices(nb_vert);

        let now = Instant::now();
        let mut del_struct = DelaunayStructure2D::new();
        del_struct.add_vertices_to_insert(&vec_pts);
        del_struct.update_delaunay()?;
        let duration = now.elapsed();
        let milli = duration.as_millis();

        println!("{} vertices: {}ms", nb_vert, milli);
    }

    Ok(())
}
