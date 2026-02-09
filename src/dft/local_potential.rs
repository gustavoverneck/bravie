use ndarray::Array3;
use nalgebra::Vector3;
use std::collections::HashMap;

use crate::core::structure::Structure;
use crate::core::fft::FftGrid;
use crate::io::upf::Pseudopotential;
use crate::utils::math::interpolate_linear; // Use a função refatorada

pub fn calculate_local_potential(
    structure: &Structure,
    fft_grid: &FftGrid,
    pseudos: &HashMap<usize, Pseudopotential>
) -> Array3<f64> {
    let (nx, ny, nz) = (fft_grid.size[0], fft_grid.size[1], fft_grid.size[2]);
    let mut v_local = Array3::<f64>::zeros((nx, ny, nz));

    // Pré-calcula inversa para PBC (Minimum Image Convention)
    let lattice_inv = structure.lattice.vectors.try_inverse().expect("Singular lattice");

    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let frac_pos = Vector3::new(
                    i as f64 / nx as f64,
                    j as f64 / ny as f64,
                    k as f64 / nz as f64,
                );
                
                let r_grid = structure.lattice.vectors * frac_pos;
                let mut v_val = 0.0;

                for atom in &structure.atoms {
                    let pseudo = pseudos.get(&atom.species_id)
                        .expect("Pseudopotencial não encontrado");

                    // Distância com PBC (Minimum Image Convention)
                    let diff = r_grid - atom.position;
                    let mut d_frac = lattice_inv * diff;
                    d_frac.x -= d_frac.x.round();
                    d_frac.y -= d_frac.y.round();
                    d_frac.z -= d_frac.z.round();
                    
                    let d_cart = structure.lattice.vectors * d_frac;
                    let dist = d_cart.norm();

                    // Interpola o Potencial Local
                    // Diferente da densidade, aqui não dividimos por 4*pi*r^2.
                    // O UPF já traz V_loc(r) pronto (em Ry).
                    v_val += interpolate_linear(dist, &pseudo.mesh, &pseudo.local);
                }

                v_local[[i, j, k]] = v_val;
            }
        }
    }

    v_local
}