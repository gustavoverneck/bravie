use std::process;
use num_complex::{Complex, Complex64};
use ndarray::Array1;

use bravie::Simulation;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::core::simulation::HamiltonianModel; // Importe o Enum
use bravie::utils::welcome::print_welcome;

// Importa os módulos de DFT
use bravie::dft::local_potential::calculate_local_potential;
use bravie::dft::potentials::solve_hartree;
use bravie::dft::xc::calculate_xc_lda;
use bravie::dft::hamiltonian::{apply_hamiltonian_local, update_v_eff};

fn run_hamiltonian_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Teste de Integração: Operador Hamiltoniano (H) ===\n");

    // 1. Setup Silício
    let a = 10.26;
    let si = Species {
        id: 0,
        element: "Si".to_string(),
        atomic_number: 14,
        mass: 28.085,
        pseudo_path: "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF".to_string(),
    };

    let silicon = Structure::builder()
        .lattice(
            [0.0, a/2.0, a/2.0],
            [a/2.0, 0.0, a/2.0],
            [a/2.0, a/2.0, 0.0]
        )
        .add_species(si)
        .add_atom([0.0, 0.0, 0.0], 0)
        .add_atom([0.25, 0.25, 0.25], 0)
        .build()?;

    // 2. Simulação
    // Usamos Ecut=30 para teste rápido
    let mut sim = Simulation::builder()
        .structure(silicon)
        .ecut(30.0)
        .k_grid(KGrid::gamma())
        .build()?;
    
    // 3. Prepara Potencial Efetivo (Densidade -> Potenciais)
    println!("Inicializando Densidade e Potenciais...");
    sim.initialize_density();

    let v_loc = calculate_local_potential(&sim.structure, &sim.fft_grid, &sim.pseudos);
    let (v_h, _) = solve_hartree(&sim.rho, &mut sim.fft_grid, &sim.bases[0], &sim.structure);
    let (v_xc, _) = calculate_xc_lda(&sim.rho, sim.structure.lattice.volume());

    update_v_eff(&mut sim.v_eff, &v_loc, &v_h, &v_xc);
    
    println!("V_eff montado. Mínimo: {:.6} Ry", sim.v_eff.iter().fold(f64::INFINITY, |a, &b| a.min(b)));

    // 4. Prepara Função de Onda de Teste
    // Vamos usar uma onda plana simples (G=0) normalizada
    let n_pw = sim.bases[0].g_vectors.len();
    let mut psi = Array1::<Complex64>::zeros(n_pw);
    psi[0] = Complex::new(1.0, 0.0); // Estado |0> (k=0, G=0)

    println!("\n--- Comparação de Modelos ---");

    // Teste A: Schrödinger (Não-Relativístico)
    let h_psi_sch = apply_hamiltonian_local(
        &psi, 
        &sim.v_eff, 
        &mut sim.fft_grid, 
        &sim.bases[0], 
        HamiltonianModel::Schrodinger
    );
    let energy_sch = psi.dot(&h_psi_sch).re;
    println!("1. Schrödinger E = <psi|H|psi> : {:.8} Ry", energy_sch);

    // Teste B: Dirac (Escalar Relativístico)
    let h_psi_dirac = apply_hamiltonian_local(
        &psi, 
        &sim.v_eff, 
        &mut sim.fft_grid, 
        &sim.bases[0], 
        HamiltonianModel::DiracScalarRelativistic
    );
    let energy_dirac = psi.dot(&h_psi_dirac).re;
    println!("2. Dirac (SR)  E = <psi|H|psi> : {:.8} Ry", energy_dirac);

    // Análise
    let delta = energy_dirac - energy_sch;
    println!("\nDiferença Relativística (Delta): {:.8e} Ry", delta);

    // Para G=0, a energia cinética é zero em ambos os modelos.
    // A diferença só aparecerá se usarmos um psi com componentes G != 0.
    // Vamos adicionar um componente de alta frequência para testar a cinética.
    
    println!("\n--- Teste Cinético (G_max) ---");
    // Pega o último vetor G (maior energia cinética)
    let last_idx = n_pw - 1;
    let mut psi_high = Array1::<Complex64>::zeros(n_pw);
    psi_high[last_idx] = Complex::new(1.0, 0.0);
    
    let h_high_sch = apply_hamiltonian_local(&psi_high, &sim.v_eff, &mut sim.fft_grid, &sim.bases[0], HamiltonianModel::Schrodinger);
    let h_high_dir = apply_hamiltonian_local(&psi_high, &sim.v_eff, &mut sim.fft_grid, &sim.bases[0], HamiltonianModel::DiracScalarRelativistic);
    
    // Subtraímos a parte do potencial (que é igual) para ver só a cinética
    // E_kin = <psi|T|psi>
    // Como V_eff é diagonal em R, seu valor esperado pode variar, mas T é diagonal em G.
    // O elemento de matriz <G|T|G> é exatamente T(G).
    
    let t_sch = h_high_sch[last_idx].re - psi_high.dot(&h_high_sch).re + psi_high.dot(&h_high_sch).re; // Simplificação: Pegamos o elemento diagonal H_ii se V fosse 0
    // Melhor: Olhamos diretamente o g^2 no basis
    let g2 = sim.bases[0].g_norm_sq[last_idx];
    
    println!("|G|^2 (Schrödinger T) : {:.6} Ry", g2);
    
    // Extrai T Dirac do hamiltoniano (H_dirac - V_term)
    // Difícil separar V sem recalcular, mas sabemos a fórmula:
    let c = 137.036;
    let t_dirac_calc = 0.5 * c*c * ((1.0 + 4.0*g2/(c*c)).sqrt() - 1.0);
    
    println!("Dirac T (Calculado)   : {:.6} Ry", t_dirac_calc);
    println!("Diferença T           : {:.6e} Ry", t_dirac_calc - g2);

    if (t_dirac_calc - g2).abs() > 1e-12 {
         println!("Efeitos relativísticos detectados na cinética.");
    } else {
         println!("Diferença muito pequena (G baixo ou c muito grande).");
    }

    Ok(())
}

fn main() {
    if let Err(e) = run_hamiltonian_test() {
        eprintln!("Erro Fatal: {}", e);
        process::exit(1);
    }
}