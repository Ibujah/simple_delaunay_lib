#[cfg(test)]
mod delaunay_3d_test {
    use anyhow::Result;
    use env_logger;
    use rand::Rng;
    use std::time::Instant;

    use delaunay_lib::delaunay_3d::delaunay_struct_3d;

    #[ctor::ctor]
    fn init() {
        env_logger::init();
    }

    fn create_and_check_delaunay(vec_pts: &Vec<[f64; 3]>) -> Result<()> {
        let now = Instant::now();
        let mut del_struct = delaunay_struct_3d::DelaunayStructure3D::new();
        del_struct.add_vertices_to_insert(&vec_pts);
        del_struct.update_delaunay()?;
        let duration = now.elapsed();
        let milli = duration.as_millis();

        log::info!("Delaunay computed in {}ms", milli);

        log::info!("Checking delaunay");
        assert!(del_struct.is_valid()?);
        Ok(())
    }

    #[test]
    fn test_random() -> Result<()> {
        let mut rng = rand::thread_rng();

        let mut vec_pts: Vec<[f64; 3]> = Vec::new();
        for _ in 0..1000 {
            let (x, y, z): (f64, f64, f64) = rng.gen();
            vec_pts.push([x, y, z]);
        }
        create_and_check_delaunay(&vec_pts)?;
        Ok(())
    }

    #[test]
    fn test_regular() -> Result<()> {
        let mut vec_pts: Vec<[f64; 3]> = Vec::new();
        let mut vec_inds: Vec<usize> = Vec::new();
        let nb = 10;
        for ind in 0..(nb * nb * nb) {
            let ind1 = ind % nb;
            let ind2 = (ind / nb) % nb;
            let ind3 = ind / (nb * nb);

            let x = (ind1 as f64) / (nb as f64);
            let y = (ind2 as f64) / (nb as f64);
            let z = (ind3 as f64) / (nb as f64);
            vec_pts.push([x, y, z]);
            vec_inds.push(ind);
        }
        create_and_check_delaunay(&vec_pts)?;
        Ok(())
    }
}
