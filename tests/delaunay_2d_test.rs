#[cfg(test)]
mod delaunay_2d_test {
    use anyhow::Result;
    use env_logger;
    use rand::Rng;
    use simple_delaunay_lib::delaunay_2d::delaunay_struct_2d;
    use std::time::Instant;

    #[ctor::ctor]
    fn init() {
        env_logger::init();
    }

    fn create_and_check_delaunay(vec_pts: &Vec<[f64; 2]>) -> Result<()> {
        let now = Instant::now();
        let mut del_struct = delaunay_struct_2d::DelaunayStructure2D::new();
        del_struct.insert_vertices(&vec_pts, true)?;
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

        let mut vec_pts: Vec<[f64; 2]> = Vec::new();
        for _ in 0..1000 {
            let (x, y): (f64, f64) = rng.gen();
            vec_pts.push([x, y]);
        }
        create_and_check_delaunay(&vec_pts)?;
        Ok(())
    }

    #[test]
    fn test_regular() -> Result<()> {
        let mut vec_pts: Vec<[f64; 2]> = Vec::new();
        let mut vec_inds: Vec<usize> = Vec::new();
        for ind in 0..1000 {
            let ind1 = ind % 100;
            let ind2 = ind / 100;

            let x = (ind1 as f64) / 100.;
            let y = (ind2 as f64) / 100.;
            vec_pts.push([x, y]);
            vec_inds.push(ind);
        }
        create_and_check_delaunay(&vec_pts)?;
        Ok(())
    }

    #[test]
    fn test_update() -> Result<()> {
        let mut rng = rand::thread_rng();

        let mut vec_pts: Vec<[f64; 2]> = Vec::new();
        let mut vec_inds: Vec<usize> = Vec::new();
        for ind in 0..1000 {
            let (x, y): (f64, f64) = rng.gen();
            vec_pts.push([x, y]);
            vec_inds.push(ind);
        }
        let now = Instant::now();
        let mut del_struct = delaunay_struct_2d::DelaunayStructure2D::new();
        del_struct.insert_vertices(&vec_pts, true)?;
        let duration = now.elapsed();
        let milli = duration.as_millis();

        log::info!("Delaunay computed in {}ms", milli);

        log::info!("Checking delaunay");
        assert!(del_struct.is_valid()?);

        let mut vec_pts: Vec<[f64; 2]> = Vec::new();
        let mut vec_inds: Vec<usize> = Vec::new();
        for ind in 0..1000 {
            let (x, y): (f64, f64) = rng.gen();
            vec_pts.push([x, y]);
            vec_inds.push(ind);
        }
        let now = Instant::now();
        del_struct.insert_vertices(&vec_pts, true)?;
        let duration = now.elapsed();
        let milli = duration.as_millis();

        log::info!("Delaunay update computed in {}ms", milli);

        log::info!("Checking delaunay");
        assert!(del_struct.is_valid()?);
        Ok(())
    }
}
