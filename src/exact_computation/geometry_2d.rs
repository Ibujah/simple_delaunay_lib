use super::float_ops::{product_f32, product_f64, sign_of_a_sum_f64};

pub fn ccw(pts: [[f32; 2]; 3]) -> i32 {
    let [[x1, y1], [x2, y2], [x3, y3]] = pts;

    let x1y2 = product_f32(x1, y2);
    let x2y3 = product_f32(x2, y3);
    let x3y1 = product_f32(x3, y1);
    let x1y3 = -product_f32(x1, y3);
    let x2y1 = -product_f32(x2, y1);
    let x3y2 = -product_f32(x3, y2);

    let vals = vec![x1y2, x2y3, x3y1, x1y3, x2y1, x3y2];

    let sign = sign_of_a_sum_f64(&vals);

    sign
}

pub fn incircle(pts: [[f32; 2]; 3], pt: [f32; 2]) -> i32 {
    fn sub_det(
        sign: i8,
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        x3: f32,
        y3: f32,
        x4: f32,
        y4: f32,
    ) -> Vec<f64> {
        let x1x1 = if sign >= 0 {
            product_f32(x1, x1)
        } else {
            -product_f32(x1, x1)
        };
        let y1y1 = if sign >= 0 {
            product_f32(y1, y1)
        } else {
            -product_f32(y1, y1)
        };

        let x3y4 = product_f32(x3, y4);
        let x4y3 = -product_f32(x4, y3);

        let x4y2 = product_f32(x4, y2);
        let x2y4 = -product_f32(x2, y4);

        let x2y3 = product_f32(x2, y3);
        let x3y2 = -product_f32(x3, y2);

        let mut vals = Vec::new();
        vals.extend(product_f64(x1x1, x3y4));
        vals.extend(product_f64(x1x1, x4y3));
        vals.extend(product_f64(x1x1, x4y2));
        vals.extend(product_f64(x1x1, x2y4));
        vals.extend(product_f64(x1x1, x2y3));
        vals.extend(product_f64(x1x1, x3y2));
        vals.extend(product_f64(y1y1, x3y4));
        vals.extend(product_f64(y1y1, x4y3));
        vals.extend(product_f64(y1y1, x4y2));
        vals.extend(product_f64(y1y1, x2y4));
        vals.extend(product_f64(y1y1, x2y3));
        vals.extend(product_f64(y1y1, x3y2));
        vals
    }

    let [[x1, y1], [x2, y2], [x3, y3]] = pts;
    let [x4, y4] = pt;

    let mut vals = Vec::new();
    vals.extend(sub_det(1, x1, y1, x2, y2, x3, y3, x4, y4));
    vals.extend(sub_det(-1, x2, y2, x1, y1, x3, y3, x4, y4));
    vals.extend(sub_det(1, x3, y3, x1, y1, x2, y2, x4, y4));
    vals.extend(sub_det(-1, x4, y4, x1, y1, x2, y2, x3, y3));

    sign_of_a_sum_f64(&vals)
}
