#[cfg(test)]
mod essa_tests {
    use delaunay_lib::exact_computation::{float_ops, geometry_2d};
    use rand::Rng;

    #[test]
    fn test_circle() {
        let mut rng = rand::thread_rng();
        let (x1, y1): (f32, f32) = rng.gen();
        let (x2, y2): (f32, f32) = rng.gen();
        let (x3, y3): (f32, f32) = rng.gen();

        let sign1 = geometry_2d::incircle([[x1, y1], [x2, y2], [x3, y3]], [x1, y1]);
        let sign2 = geometry_2d::incircle([[x1, y1], [x2, y2], [x3, y3]], [x2, y2]);
        let sign3 = geometry_2d::incircle([[x1, y1], [x2, y2], [x3, y3]], [x3, y3]);

        assert_eq!(sign1, 0);
        assert_eq!(sign2, 0);
        assert_eq!(sign3, 0);
    }

    #[test]
    fn test_sum() {
        let vals = vec![1.0, 2.0, -1.0, -2.0];

        let sign = float_ops::sign_of_a_sum_f64(&vals);
        assert_eq!(sign, 0);
    }
}
