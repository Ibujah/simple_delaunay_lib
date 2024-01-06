/// Sorts vertices along 3D Hilbert curve
pub fn build_hilbert_curve_3d(vertices: &Vec<[f64; 3]>, indices_to_add: &Vec<usize>) -> Vec<usize> {
    let mut curve_order = Vec::new();

    let mut pt_min = vertices[indices_to_add[0]];
    let mut pt_max = vertices[indices_to_add[0]];

    for &ind in indices_to_add.iter() {
        if pt_min[0] > vertices[ind][0] {
            pt_min[0] = vertices[ind][0];
        }
        if pt_min[1] > vertices[ind][1] {
            pt_min[1] = vertices[ind][1];
        }
        if pt_min[2] > vertices[ind][2] {
            pt_min[2] = vertices[ind][2];
        }
        if pt_max[0] < vertices[ind][0] {
            pt_max[0] = vertices[ind][0];
        }
        if pt_max[1] < vertices[ind][1] {
            pt_max[1] = vertices[ind][1];
        }
        if pt_max[2] < vertices[ind][2] {
            pt_max[2] = vertices[ind][2];
        }
    }

    let mut to_subdiv = Vec::new();
    let indices: Vec<usize> = indices_to_add.iter().map(|&x| x).collect();
    to_subdiv.push(([0, 0, 0], 0, pt_min, pt_max, indices));

    loop {
        if let Some((start, dir, pt_min, pt_max, indices_to_add)) = to_subdiv.pop() {
            if indices_to_add.len() > 1 {
                let sep_x = (pt_min[0] + pt_max[0]) / 2.0;
                let sep_y = (pt_min[1] + pt_max[1]) / 2.0;
                let sep_z = (pt_min[2] + pt_max[2]) / 2.0;

                let mut sep_ind = [
                    [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
                    [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
                ];

                for &ind in indices_to_add.iter() {
                    let vert = vertices[ind];
                    let xind = if vert[0] < sep_x { 0 } else { 1 } as usize;
                    let yind = if vert[1] < sep_y { 0 } else { 1 } as usize;
                    let zind = if vert[2] < sep_z { 0 } else { 1 } as usize;
                    sep_ind[xind][yind][zind].push(ind);
                }

                let pt_x = [pt_min[0], sep_x, pt_max[0]];
                let pt_y = [pt_min[1], sep_y, pt_max[1]];
                let pt_z = [pt_min[2], sep_z, pt_max[2]];

                let (next_modif, dir) = match (dir, start[dir]) {
                    (0, 0) => Some(([1, 2, 1, 0, 1, 2, 1, 0], [1, 2, 2, 0, 0, 2, 2, 1])),
                    (0, 1) => Some(([2, 1, 2, 0, 2, 1, 2, 0], [2, 1, 1, 0, 0, 1, 1, 2])),
                    (1, 0) => Some(([2, 0, 2, 1, 2, 0, 2, 1], [2, 0, 0, 1, 1, 0, 0, 2])),
                    (1, 1) => Some(([0, 2, 0, 1, 0, 2, 0, 1], [0, 2, 2, 1, 1, 2, 2, 0])),
                    (2, 0) => Some(([0, 1, 0, 2, 0, 1, 0, 2], [0, 1, 1, 2, 2, 1, 1, 0])),
                    (2, 1) => Some(([1, 0, 1, 2, 1, 0, 1, 2], [1, 0, 0, 2, 2, 0, 0, 1])),
                    (_, _) => None,
                }
                .unwrap();

                let mut sep_subind = start;
                let mut start_ind = start;
                for i in 0..8 {
                    let mut vec_inds = Vec::new();
                    vec_inds.append(&mut sep_ind[sep_subind[0]][sep_subind[1]][sep_subind[2]]);
                    to_subdiv.push((
                        start_ind,
                        dir[i],
                        [
                            pt_x[sep_subind[0]],
                            pt_y[sep_subind[1]],
                            pt_z[sep_subind[2]],
                        ],
                        [
                            pt_x[sep_subind[0] + 1],
                            pt_y[sep_subind[1] + 1],
                            pt_z[sep_subind[2] + 1],
                        ],
                        vec_inds,
                    ));
                    sep_subind[next_modif[i]] = 1 - sep_subind[next_modif[i]];
                    start_ind[next_modif[i]] = 1 - start_ind[next_modif[i]];
                    start_ind[dir[i]] = 1 - start_ind[dir[i]];
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
