use std::process;
use bravie::Simulation;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::dft::local_potential::calculate_local_potential;
use bravie::dft::potentials::solve_hartree;
use bravie::dft::xc::calculate_xc_lda;
use bravie::dft::hamiltonian::update_v_eff;
use bravie::dft::solver::solve_bands; // Importe o novo módulo

fn run_solver_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Teste de Integração: Solver de Bandas ===\n");

    // 1. Setup Silício (Mesmo de antes)
    let a = 10.26;
    let si = Species { id: 0, element: "Si".to_string(), atomic_number: 14, mass: 28.085, pseudo_path: "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF".to_string() };
    let silicon = Structure::builder()
        .lattice([0.0, a/2.0, a/2.0], [a/2.0, 0.0, a/2.0], [a/2.0, a/2.0, 0.0])
        .add_species(si).add_atom([0.0, 0.0, 0.0], 0).add_atom([0.25, 0.25, 0.25], 0)
        .build()?;

    let mut sim = Simulation::builder()
        .structure(silicon).ecut(30.0).k_grid(KGrid::gamma()).build()?;
    
    // 2. Prepara Potencial Fixo (SAD)
    println!("Preparando Potencial (SAD)...");
    sim.initialize_density();
    let v_loc = calculate_local_potential(&sim.structure, &sim.fft_grid, &sim.pseudos);
    let (v_h, _) = solve_hartree(&sim.rho, &mut sim.fft_grid, &sim.bases[0], &sim.structure);
    let (v_xc, _) = calculate_xc_lda(&sim.rho, sim.structure.lattice.volume());
    update_v_eff(&mut sim.v_eff, &v_loc, &v_h, &v_xc);

    // 3. Diagonaliza!
    // Silício tem 8 elétrons de valência (2 átomos * 4 e-).
    // Bandas ocupadas = 8 / 2 (spin) = 4 bandas.
    // Vamos calcular 6 para ver algumas vazias.
    let num_bands = 6;
    let result = solve_bands(
        num_bands, 
        &sim.v_eff, 
        &mut sim.fft_grid, 
        &sim.bases[0],
        sim.hamiltonian_model
    );

    println!("\n=== Espectro de Energia (Ponto Gamma) ===");
    for (i, e) in result.eigenvalues.iter().enumerate() {
        let occ = if i < 4 { "Ocupada" } else { "Vazia  " };
        println!("Banda {}: {:.6} Ry  [{}]", i+1, e, occ);
    }
    
    // Cálculo do GAP
    let homo = result.eigenvalues[3];
    let lumo = result.eigenvalues[4];
    println!("\nGap Estimado (Gamma): {:.4} eV", (lumo - homo) * 13.605); // 1 Ry = 13.605 eV

    Ok(())
}

fn main() {
    if let Err(e) = run_solver_test() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}