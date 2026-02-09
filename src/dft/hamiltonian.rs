use ndarray::{Array1, Array3};
use num_complex::{Complex, Complex64};
use crate::core::fft::FftGrid;
use crate::core::basis::PlaneWaveBasis;
use crate::core::simulation::HamiltonianModel;

// Velocidade da luz em unidades atômicas (aproximada)
const SPEED_OF_LIGHT: f64 = 137.035999;

/// Aplica o Hamiltoniano Local (H_loc = T + V_eff) em uma função de onda.
///
/// # Argumentos
/// * `psi_g` - Coeficientes da função de onda no espaço recíproco.
/// * `v_eff` - Potencial efetivo total no espaço real.
/// * `fft_grid` - Grid FFT (usado como buffer temporário).
/// * `basis` - Base de ondas planas (contém |k+G|^2).
/// * `model` - Enum definindo se usamos Schrödinger ou Dirac.
///
/// Retorna: O vetor resultante da aplicação H * psi (no espaço recíproco).
pub fn apply_hamiltonian_local(
    psi_g: &Array1<Complex64>,
    v_eff: &Array3<f64>,
    fft_grid: &mut FftGrid,
    basis: &PlaneWaveBasis,
    model: HamiltonianModel
) -> Array1<Complex64> {
    let n_g = basis.g_vectors.len();
    let mut h_psi = Array1::<Complex64>::zeros(n_g);

    // 1. Aplica Energia Cinética (Operador T no Espaço Recíproco)
    // O operador T é diagonal em G: T|psi> = T(G) * psi(G)
    for i in 0..n_g {
        // g2 = |k + G|^2
        let g2 = basis.g_norm_sq[i];

        let kinetic_energy = match model {
            HamiltonianModel::Schrodinger => {
                // Modelo Clássico (Não-Relativístico)
                // T = p^2 (em unidades Rydberg: p^2 = k^2)
                g2
            },
            HamiltonianModel::DiracScalarRelativistic => {
                // Aproximação Escalar Relativística
                // T = (c^2 / 2) * (sqrt(1 + 4*k^2/c^2) - 1)
                let c = SPEED_OF_LIGHT;
                let c2 = c * c;
                
                // Fórmula para T relativístico em Ry
                0.5 * c2 * ((1.0 + 4.0 * g2 / c2).sqrt() - 1.0)
            }
        };

        h_psi[i] = Complex::new(kinetic_energy, 0.0) * psi_g[i];
    }

    // 2. Aplica Potencial Local (Operador V no Espaço Real)
    // V(r) é diagonal no espaço real. Usamos FFT para aplicar.
    // Operação: FFT_inv(psi) -> V_eff(r) * psi(r) -> FFT_fwd -> soma em h_psi
    
    // A. Transformada Inversa: G -> r
    // Leva a função de onda para o grid real
    fft_grid.to_real_space(psi_g);
    
    // B. Multiplicação Ponto-a-Ponto: psi(r) = psi(r) * V_eff(r)
    let (nx, ny, nz) = (fft_grid.size[0], fft_grid.size[1], fft_grid.size[2]);
    
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // psi(r) atual (no buffer do grid)
                let psi_r = fft_grid.buffer[[i, j, k]];
                
                // Potencial efetivo V(r)
                let v_r = v_eff[[i, j, k]];
                
                // Aplicação do operador V: psi_new = V * psi
                fft_grid.buffer[[i, j, k]] = psi_r * v_r;
            }
        }
    }
    
    // C. Transformada Direta: r -> G
    // Traz o resultado V*psi de volta para o espaço recíproco
    let mut v_psi_g = Array1::<Complex64>::zeros(n_g);
    fft_grid.to_recip_space(&mut v_psi_g);

    // 3. Soma Final: H|psi> = T|psi> + V|psi>
    // h_psi já contém a parte cinética, agora somamos a parte do potencial
    // Aplica o fator de escala apenas na parte do potencial que veio da FFT
    for i in 0..n_g {
        // v_psi_g já está na escala correta compatível com a energia cinética
        h_psi[i] = h_psi[i] + v_psi_g[i]; 
    }

    h_psi
}

/// Helper para atualizar o Potencial Efetivo Total
/// V_eff = V_local + V_Hartree + V_xc
pub fn update_v_eff(
    v_eff: &mut Array3<f64>,
    v_local: &Array3<f64>,
    v_hartree: &Array3<f64>,
    v_xc: &Array3<f64>
) {
    // Soma ponto a ponto usando iteradores para performance (evita bounds checking excessivo)
    v_eff.iter_mut()
        .zip(v_local.iter())
        .zip(v_hartree.iter())
        .zip(v_xc.iter())
        .for_each(|(((ve, vl), vh), vx)| {
            *ve = vl + vh + vx;
        });
}