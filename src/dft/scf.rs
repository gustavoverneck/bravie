use ndarray::Array3;
use crate::Simulation;
use crate::dft::local_potential::calculate_local_potential;
use crate::dft::potentials::solve_hartree;
use crate::dft::xc::calculate_xc_lda;
use crate::dft::hamiltonian::update_v_eff;
use crate::dft::solver::solve_bands;
use crate::dft::density::compute_density_from_wavefunctions;
use crate::dft::mixing::AndersonMixer;

/// Parâmetros de controle do ciclo SCF
pub struct ScfParameters {
    pub max_iter: usize,
    pub tol_energy: f64, // Ry
    pub tol_rho: f64,    // Diferença na densidade
    pub mixing_beta: f64, // Fator de mistura (0.1 = 10% novo, 90% velho)
    pub mixing_history: usize,
}

impl Default for ScfParameters {
    fn default() -> Self {
        Self {
            max_iter: 200,
            tol_energy: 1e-6,
            tol_rho: 1e-5,
            mixing_beta: 0.1,
            mixing_history: 5,
        }
    }
}

pub fn run_scf_loop(sim: &mut Simulation, params: ScfParameters) {
    println!("\n=== Iniciando Ciclo Auto-Consistente (SCF) ===");
    println!("Max Iters: {}, Tol E: {:.1e} Ry, Beta: {:.2}", params.max_iter, params.tol_energy, params.mixing_beta);

    // 1. Inicializa Densidade (SAD) - Já deve ter sido feito fora ou garantimos aqui
    // sim.initialize_density(); // Assumindo que já foi chamado

    // 2. Pré-calcula V_loc (não muda durante o SCF)
    let v_loc = calculate_local_potential(&sim.structure, &sim.fft_grid, &sim.pseudos);
    
    let mut prev_energy = 0.0;
    
    // Número de elétrons total
    let mut n_electrons = 0.0;
    for atom in &sim.structure.atoms {
        if let Some(p) = sim.pseudos.get(&atom.species_id) {
            n_electrons += p.header.z_valence;
        }
    }
    // Ocupações: Assumindo isolante/semicondutor spin-degenerado (f=2.0)
    // Bandas ocupadas = N_el / 2
    let n_bands_occ = (n_electrons / 2.0).ceil() as usize;
    let n_bands_total = n_bands_occ + 4; // Calcula algumas vazias extra
    
    let mut occupations = vec![0.0; n_bands_total];
    for i in 0..n_bands_occ {
        occupations[i] = 2.0; // 2 elétrons por banda
    }

    // Inicializa o Mixer
    let mut mixer = AndersonMixer::new(params.mixing_beta, params.mixing_history);

    for iter in 1..=params.max_iter {
        // A. Calcula Potenciais Dependentes de Rho
        let (v_h, e_h) = solve_hartree(&sim.rho, &mut sim.fft_grid, &sim.bases[0], &sim.structure);
        let (v_xc, e_xc) = calculate_xc_lda(&sim.rho, sim.structure.lattice.volume()); // Volume

        // B. Atualiza V_eff
        update_v_eff(&mut sim.v_eff, &v_loc, &v_h, &v_xc);

        // C. Diagonaliza (Solver)
        // Passamos o V_eff atual e obtemos novos autovetores
        let bands = solve_bands(
            n_bands_total, 
            &sim.v_eff, 
            &mut sim.fft_grid, 
            &sim.bases[0],
            sim.hamiltonian_model
        );

        // Ordena bandas (IMPORTANTE!)
        let mut indices: Vec<usize> = (0..n_bands_total).collect();
        indices.sort_by(|&i, &j| bands.eigenvalues[i].partial_cmp(&bands.eigenvalues[j]).unwrap());
        
        // D. Calcula Energia Total (Harris-Foulkes ou Kohn-Sham direto)
        // E_total = Sum(epsilon_occ) - E_H - E_xc + E_H_rho + E_xc_rho + E_ewald...
        // Forma simplificada: Soma dos autovalores ocupados - dupla contagem
        // E_band = sum(occ * epsilon)
        let mut e_band = 0.0;
        for i in 0..n_bands_occ {
             // Usa índice ordenado
             let idx = indices[i];
             e_band += 2.0 * bands.eigenvalues[idx];
        }
        
        // Termos de correção de dupla contagem (Double Counting)
        // E_tot = E_band - integral(V_H * rho)/2 - integral(V_xc * rho) + E_xc(rho) + E_ewald
        // Por simplicidade agora, vamos monitorar E_band que deve diminuir e convergir.
        let current_energy = e_band; // Placeholder para métrica de convergência

        // E. Calcula Nova Densidade
        // Reordena vetores para passar para a função de densidade corretamente
        let sorted_eigenvectors: Vec<_> = indices.iter().map(|&i| bands.eigenvectors[i].clone()).collect();
        
        let mut rho_new = compute_density_from_wavefunctions(
            &sorted_eigenvectors, 
            &mut sim.fft_grid, 
            &occupations
        );

        // Renormaliza Rho_new para ter exatamente N_electrons
        let vol = sim.structure.lattice.volume();
        let n_grid = sim.rho.len() as f64;
        let dvol = vol / n_grid;
        let charge_new = rho_new.sum() * dvol;
        if charge_new > 1e-6 {
            rho_new.mapv_inplace(|x| x * (n_electrons / charge_new));
        }

        // F. Mistura (Linear Mixing)
        // rho_next = beta * rho_new + (1-beta) * rho_old
        let beta = params.mixing_beta;

        let rho_mixed = mixer.mix(&sim.rho, &rho_new);
        // Calcula erro da densidade (RMS)
        let rho_diff = &rho_new - &sim.rho;
        let rho_err = (rho_diff.mapv(|x| x*x).sum() * dvol).sqrt();

        // Atualiza sim.rho
        sim.rho = rho_mixed;

        // G. Verifica Convergência
        let e_diff = (current_energy - prev_energy).abs();
        
        println!("SCF {:2} | E_band: {:.6} Ry | dE: {:.1e} | dRho: {:.1e}", 
            iter, current_energy, e_diff, rho_err);

        if iter > 1 && e_diff < params.tol_energy && rho_err < params.tol_rho {
            println!("Convergência alcançada em {} iterações!", iter);
            break;
        }
        
        prev_energy = current_energy;
    }
}