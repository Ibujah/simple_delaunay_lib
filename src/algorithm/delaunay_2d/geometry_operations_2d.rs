use nalgebra::base::*;

pub fn circle_center_and_radius(
    pt1: &Vector2<f32>,
    pt2: &Vector2<f32>,
    pt3: &Vector2<f32>,
) -> Option<(Vector2<f32>, f32)> {
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

pub fn circle_center_and_radius_type2(
    pt1: &Vector2<f32>,
    pt2: &Vector2<f32>,
    pt3: &Vector2<f32>,
) -> Option<(Vector2<f32>, f32)> {
    // order points
    let vec12 = pt2 - pt1;
    let vec13 = pt3 - pt1;
    let (p1, p2, p3) = if vec12.dot(&vec13) > 0. {
        if vec12.norm() > vec13.norm() {
            (pt1, pt3, pt2)
        } else {
            (pt1, pt2, pt3)
        }
    } else {
        (pt2, pt1, pt3)
    };

    let mid12 = (p1 + p2) / 2.0;
    let mid13 = (p1 + p3) / 2.0;
    let nor12 = Vector2::new(p1.y - p2.y, p2.x - p1.x).normalize();
    let nor13 = Vector2::new(p1.y - p3.y, p3.x - p1.x).normalize();

    let mat = Matrix2::new(-nor12.x, nor13.x, -nor12.y, nor13.y);
    let b = mid12 - mid13;

    let opt_fac = mat.lu().solve(&b);

    if let Some(fac) = opt_fac {
        let ctr1 = mid12 + fac[0] * nor12;
        let ctr2 = mid13 + fac[1] * nor13;
        let center = (ctr1 + ctr2) / 2.0;
        let radius1 = (center - pt1).norm();
        let radius2 = (center - pt2).norm();
        let radius3 = (center - pt3).norm();
        let radius = if radius1 > radius2 {
            if radius3 > radius1 {
                radius3
            } else {
                radius1
            }
        } else {
            if radius3 > radius2 {
                radius3
            } else {
                radius2
            }
        };
        Some((center, radius))
    } else {
        None
    }
}

pub fn line_normal_and_factor(pt1: &Vector2<f32>, pt2: &Vector2<f32>) -> (Vector2<f32>, f32) {
    let vec = (pt1 - pt2).normalize();
    let normal = Vector2::new(-vec[1], vec[0]);
    let factor = normal.dot(pt1);
    (normal, factor)
}

pub fn line_normal_and_factor_excluding(
    pt1: &Vector2<f32>,
    pt2: &Vector2<f32>,
    pt_excl: &Vector2<f32>,
) -> (Vector2<f32>, f32) {
    let (normal, factor) = line_normal_and_factor(pt1, pt2);
    if normal.dot(pt_excl) + factor <= 0.0 {
        (-normal, -factor)
    } else {
        (normal, factor)
    }
}

pub fn build_hilbert_curve(
    vertices: &Vec<Vector2<f64>>,
    indices_to_add: &Vec<usize>,
) -> Vec<usize> {
    let mut curve_order = Vec::new();

    let mut pt_min = vertices[indices_to_add[0]];
    let mut pt_max = vertices[indices_to_add[0]];

    for &ind in indices_to_add.iter() {
        if pt_min.x > vertices[ind].x {
            pt_min.x = vertices[ind].x;
        }
        if pt_min.y > vertices[ind].y {
            pt_min.y = vertices[ind].y;
        }
        if pt_max.x < vertices[ind].x {
            pt_max.x = vertices[ind].x;
        }
        if pt_max.y < vertices[ind].y {
            pt_max.y = vertices[ind].y;
        }
    }

    let mut to_subdiv = Vec::new();
    let indices: Vec<usize> = indices_to_add.iter().map(|&x| x).collect();
    to_subdiv.push((0, pt_min, pt_max, indices));

    loop {
        if let Some((rot, pt_min, pt_max, indices_to_add)) = to_subdiv.pop() {
            if indices_to_add.len() > 1 {
                let sep_x = (pt_min.x + pt_max.x) / 2.0;
                let sep_y = (pt_min.y + pt_max.y) / 2.0;

                let mut ind_a = Vec::new();
                let mut ind_b = Vec::new();
                let mut ind_c = Vec::new();
                let mut ind_d = Vec::new();

                for &ind in indices_to_add.iter() {
                    let vert = vertices[ind];
                    if vert.x < sep_x {
                        if vert.y < sep_y {
                            ind_a.push(ind);
                        } else {
                            ind_b.push(ind);
                        }
                    } else {
                        if vert.y < sep_y {
                            ind_d.push(ind);
                        } else {
                            ind_c.push(ind);
                        }
                    }
                }

                let pt_a_min = pt_min;
                let pt_a_max = Vector2::new(sep_x, sep_y);

                let pt_b_min = Vector2::new(pt_min.x, sep_y);
                let pt_b_max = Vector2::new(sep_x, pt_max.y);

                let pt_c_min = Vector2::new(sep_x, sep_y);
                let pt_c_max = pt_max;

                let pt_d_min = Vector2::new(sep_x, pt_min.y);
                let pt_d_max = Vector2::new(pt_max.x, sep_y);

                if rot == 0 {
                    to_subdiv.push((3, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((0, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((0, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((7, pt_d_min, pt_d_max, ind_d));
                } else if rot == 1 {
                    to_subdiv.push((6, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((1, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((1, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((2, pt_a_min, pt_a_max, ind_a));
                } else if rot == 2 {
                    to_subdiv.push((5, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((2, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((2, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((1, pt_a_min, pt_a_max, ind_a));
                } else if rot == 3 {
                    to_subdiv.push((0, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((3, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((3, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((4, pt_b_min, pt_b_max, ind_b));
                } else if rot == 4 {
                    to_subdiv.push((7, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((4, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((4, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((3, pt_b_min, pt_b_max, ind_b));
                } else if rot == 5 {
                    to_subdiv.push((2, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((5, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((5, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((6, pt_c_min, pt_c_max, ind_c));
                } else if rot == 6 {
                    to_subdiv.push((1, pt_d_min, pt_d_max, ind_d));
                    to_subdiv.push((6, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((6, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((5, pt_c_min, pt_c_max, ind_c));
                } else if rot == 7 {
                    to_subdiv.push((4, pt_c_min, pt_c_max, ind_c));
                    to_subdiv.push((7, pt_b_min, pt_b_max, ind_b));
                    to_subdiv.push((7, pt_a_min, pt_a_max, ind_a));
                    to_subdiv.push((0, pt_d_min, pt_d_max, ind_d));
                }
            } else if indices_to_add.len() == 1 {
                curve_order.push(indices_to_add[0]);
            }
        } else {
            break;
        }
    }

    curve_order
}
