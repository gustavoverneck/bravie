use ndarray::{Array1, Array2, Axis};
use num_complex::{Complex, Complex64};
use crate::core::fft::FftGrid;
use crate::core::basis::PlaneWaveBasis;
use crate::core::simulation::HamiltonianModel;
use crate::dft::hamiltonian::apply_hamiltonian_local;

/// Resultado da Diagonalização
pub struct BandSolverResult {
    pub eigenvalues: Vec<f64>,        // Energias (Ry)
    pub eigenvectors: Vec<Array1<Complex64>>, // Funções de onda (Coeffs G)
}

/// Resolve H * psi = E * psi para N bandas usando Steepest Descent pré-condicionado.
/// 
/// Algoritmo Simplificado (Band-by-Band):
/// 1. Gera chute inicial aleatório
/// 2. Ortogonaliza contra bandas inferiores já convergidas (Gram-Schmidt)
/// 3. Loop de Iteração:
///    a. Calcula H|psi>
///    b. Energia E = <psi|H|psi>
///    c. Resíduo R = (H - E)|psi>
///    d. Verifica convergência (|R| < tol)
///    e. Passo de descida: |psi_new> = |psi> - alpha * K * R
///       (Onde K é o precondicionador ~ 1/G^2)
///    f. Normaliza e Ortogonaliza
pub fn solve_bands(
    num_bands: usize,
    v_eff: &ndarray::Array3<f64>,
    fft_grid: &mut FftGrid,
    basis: &PlaneWaveBasis,
    model: HamiltonianModel
) -> BandSolverResult {
    let n_pw = basis.g_vectors.len();
    let max_iter = 200;
    let tol = 1e-10; // Tolerância de convergência (Ry)
    let alpha = 0.3; // Tamanho do passo de mistura

    let mut eigenvalues: Vec<f64> = Vec::new();
    let mut eigenvectors: Vec<ndarray::ArrayBase<ndarray::OwnedRepr<Complex<f64>>, ndarray::Dim<[usize; 1]>, Complex<f64>>> = Vec::new();

    println!("Iniciando solver iterativo para {} bandas...", num_bands);

    for b in 0..num_bands {
        // 1. Chute Inicial (Aleatório ou G=0 se b=0)
        let mut psi = Array1::<Complex64>::zeros(n_pw);
        if b == 0 {
            psi[0] = Complex::new(1.0, 0.0); // Estado Gamma
        } else {
            // Adiciona ruído nos primeiros componentes para quebrar simetria
             for i in 0..10.min(n_pw) {
                psi[i] = Complex::new(0.1 * (b as f64), 0.0);
             }
        }
        normalize(&mut psi);

        let mut energy = 0.0;
        
        // Loop SCF do Solver (Iterative Diagonalization)
        for iter in 0..max_iter {
            // A. Ortogonalização (Gram-Schmidt) contra bandas anteriores
            for prev_psi in &eigenvectors {
                let overlap = prev_psi.dot(&psi); // <prev|curr>
                // |curr> = |curr> - <prev|curr> * |prev>
                psi = psi - prev_psi.mapv(|x| x * overlap);
            }
            normalize(&mut psi);

            // B. Aplica Hamiltoniano
            let h_psi = apply_hamiltonian_local(&psi, v_eff, fft_grid, basis, model);

            // C. Calcula Energia (Valor Esperado)
            // E = <psi | H | psi>
            // Conjugado não é necessário se usarmos dot product hermitiano corretamente,
            // mas ndarray::dot não conjuga automaticamente o primeiro argumento.
            // Correção: sum( conj(psi_i) * h_psi_i )
            let e_val = psi.iter().zip(h_psi.iter())
                .map(|(p, hp)| p.conj() * hp)
                .sum::<Complex64>().re;
            
            energy = e_val;

            // D. Calcula Resíduo: R = H|psi> - E|psi>
            let residue = &h_psi - &psi.mapv(|x| x * Complex::new(e_val, 0.0));
            let error = residue.mapv(|x| x.norm_sqr()).sum().sqrt();

            if error < tol {
                if b == 0 || iter % 10 == 0 {
                    println!("  Banda {}: Conv. em {} iters. E = {:.6} Ry (Err: {:.1e})", b+1, iter, energy, error);
                }
                break;
            }

            // E. Passo de Descida com Precondicionamento (K)
            // K ~ 1 / (1 + |G|^2)  (Amortece altas frequências)
            // psi_new = psi - alpha * K * R
            for i in 0..n_pw {
                let g2 = basis.g_norm_sq[i];
                // Precondicionador diagonal simples Teter-Payne-Allan
                // Evita divisão por zero se g2=0. Adiciona shift na energia cinética.
                let preconditioner = 1.0 / (1.0 + g2); 
                
                psi[i] = psi[i] - residue[i] * Complex::new(alpha * preconditioner, 0.0);
            }
        }

        eigenvectors.push(psi);
        eigenvalues.push(energy);
    }

    BandSolverResult {
        eigenvalues,
        eigenvectors,
    }
}

fn normalize(psi: &mut Array1<Complex64>) {
    let norm = psi.mapv(|x| x.norm_sqr()).sum().sqrt();
    *psi /= Complex::new(norm, 0.0);
}