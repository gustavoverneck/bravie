use std::f64::consts::PI;
use ndarray::{Array1, Array3};
use num_complex::{Complex, Complex64};
use crate::core::fft::FftGrid;
use crate::core::basis::PlaneWaveBasis;
use crate::core::structure::Structure;

/// Resolve a Equação de Poisson para obter o Potencial de Hartree.
/// Retorna (V_hartree_real, Energia_Hartree)
pub fn solve_hartree(
    rho_real: &Array3<f64>,
    fft_grid: &mut FftGrid,
    basis: &PlaneWaveBasis,
    structure: &Structure,
) -> (Array3<f64>, f64) {
    let n_g = basis.g_vectors.len();
    let volume = structure.lattice.volume(); // Certifique-se que 'volume' é pub em Structure ou Lattice

    // 1. Rho(r) -> Rho(G)
    // Precisamos colocar rho no buffer complexo do FFT Grid
    // Mapeamento seguro elemento a elemento
    let (nx, ny, nz) = (fft_grid.size[0], fft_grid.size[1], fft_grid.size[2]);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                fft_grid.buffer[[i, j, k]] = Complex::new(rho_real[[i, j, k]], 0.0);
            }
        }
    }

    // FFT Forward
    let mut rho_g = Array1::<Complex64>::zeros(n_g);
    fft_grid.to_recip_space(&mut rho_g);

    // 2. Aplica o Kernel de Poisson: 4*pi / G^2
    let mut v_hartree_g = Array1::<Complex64>::zeros(n_g);
    let prefactor = 4.0 * PI;

    for (idx, &g2) in basis.g_norm_sq.iter().enumerate() {
        if g2 < 1e-8 {
            // Singularidade G=0
            // O termo G=0 do potencial depende da convenção. 
            // Para sistemas neutros em 3D, V(G=0) = 0 é o padrão (fundo de carga positiva uniforme).
            v_hartree_g[idx] = Complex::new(0.0, 0.0);
        } else {
            v_hartree_g[idx] = rho_g[idx] * (prefactor / g2);
        }
    }

    // 3. V_H(G) -> V_H(r) (Inverse FFT)
    fft_grid.to_real_space(&v_hartree_g);

    // Extrai a parte real do potencial
    let mut v_hartree_real = Array3::<f64>::zeros((nx, ny, nz));
    let mut integral_eh = 0.0;
    let dvol = volume / (nx * ny * nz) as f64;

    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let v_val = fft_grid.buffer[[i, j, k]].re;
                let rho_val = rho_real[[i, j, k]];
                
                v_hartree_real[[i, j, k]] = v_val;
                
                // Energia de Hartree = 0.5 * integral( rho(r) * V_H(r) ) dr
                integral_eh += 0.5 * rho_val * v_val;
            }
        }
    }

    let energy_hartree = integral_eh * dvol;

    (v_hartree_real, energy_hartree)
}