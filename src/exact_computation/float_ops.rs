use num_traits::Float;

pub fn split_f64(v: f64) -> (f32, f32) {
    let v1 = v as f32;
    let v2 = (v - (v1 as f64)) as f32;
    (v1, v2)
}

pub fn product_f32(a: f32, b: f32) -> f64 {
    let a = a as f64;
    let b = b as f64;
    a * b
}

pub fn product_f64(a: f64, b: f64) -> Vec<f64> {
    let (a1, a2) = split_f64(a);
    let (b1, b2) = split_f64(b);
    vec![
        product_f32(a1, b1),
        product_f32(a1, b2),
        product_f32(a2, b1),
        product_f32(a2, b2),
    ]
}

pub fn essa_f64(pts: &Vec<f64>) -> i32 {
    let (mut vec_pos, mut vec_neg) = pts.iter().fold(
        (Vec::new(), Vec::new()),
        |(mut vec_pos, mut vec_neg), &val| {
            if val > 0. {
                vec_pos.push(val)
            } else if val < 0. {
                vec_neg.push(-val)
            }
            (vec_pos, vec_neg)
        },
    );

    vec_pos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vec_neg.sort_by(|a, b| a.partial_cmp(b).unwrap());

    loop {
        let m = vec_pos.len();
        let n = vec_neg.len();

        if m == 0 && n == 0 {
            return 0;
        }
        if m > 0 && n == 0 {
            return 1;
        }
        if m == 0 && n > 0 {
            return -1;
        }

        let a = vec_pos.pop().unwrap();
        let b = vec_neg.pop().unwrap();

        let (_, e, _) = Float::integer_decode(a);
        let (_, f, _) = Float::integer_decode(b);

        let e = (e + 53) as i32;
        let f = (f + 53) as i32;

        if a >= (n as f64) * 2.0.powi(f) {
            return 1;
        }
        if b >= (m as f64) * 2.0.powi(e) {
            return -1;
        }

        if e == f {
            if a > b {
                let ap = a - b;
                vec_pos.push(ap);
                vec_pos.sort_by(|a, b| a.partial_cmp(b).unwrap());
            } else if a < b {
                let bp = b - a;
                vec_neg.push(bp);
                vec_neg.sort_by(|a, b| a.partial_cmp(b).unwrap());
            }
        } else if e > f {
            let val = 2.0.powi(f - 1);
            let u = if b == val { val } else { 2.0.powi(f) };
            let ap = a - u;
            let app = u - b;
            vec_pos.push(ap);
            vec_pos.push(app);
            vec_pos.sort_by(|a, b| a.partial_cmp(b).unwrap());
        } else {
            let val = 2.0.powi(e - 1);
            let v = if a == val { val } else { 2.0.powi(e) };
            let bp = b - v;
            let bpp = v - a;
            vec_neg.push(bp);
            vec_neg.push(bpp);
            vec_neg.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }
    }
}

pub fn interval_sum_f64(vals: &Vec<f64>) -> Option<i32> {
    fn to_interval(val: f64) -> (f64, f64) {
        let eps_min = 0.9999999999999999;
        let eps_max = 1.0000000000000001;
        if val > 0. {
            (val * eps_min, val * eps_max)
        } else {
            (val * eps_max, val * eps_min)
        }
    }

    let (s_inf, s_sup) = vals.iter().fold((0., 0.), |(s_min, s_max), &v| {
        let (v_min, v_max) = to_interval(v);
        (s_min + v_min, s_max + v_max)
    });

    if s_inf > 0. {
        Some(1)
    } else if s_sup < 0. {
        Some(-1)
    } else {
        None
    }

    // let s = vals.iter().fold(0., |s, &v| s + v);

    // if s > 1e-60 {
    //     Some(1)
    // } else if s < -1e-60 {
    //     Some(-1)
    // } else {
    //     None
    // }
}

pub fn sign_of_a_sum_f64(vals: &Vec<f64>) -> i32 {
    if let Some(sign) = interval_sum_f64(&vals) {
        sign
    } else {
        essa_f64(&vals)
    }
}
