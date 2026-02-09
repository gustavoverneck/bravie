use std::process;
use bravie::Simulation;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::dft::scf::{run_scf_loop, ScfParameters};

fn run_scf_full() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    
    // 1. Setup Silício
    let a = 10.26;
    let si = Species { id: 0, element: "Si".to_string(), atomic_number: 14, mass: 28.085, pseudo_path: "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF".to_string() };
    let silicon = Structure::builder()
        .lattice([0.0, a/2.0, a/2.0], [a/2.0, 0.0, a/2.0], [a/2.0, a/2.0, 0.0])
        .add_species(si).add_atom([0.0, 0.0, 0.0], 0).add_atom([0.25, 0.25, 0.25], 0)
        .build()?;

    // Ecut razoável para convergência
    let mut sim = Simulation::builder()
        .structure(silicon)
        .ecut(30.0)
        .k_grid(KGrid::gamma())
        .build()?;
    
    // Inicializa
    sim.initialize_density();

    // 2. Roda SCF
    let mut params = ScfParameters::default();

    params.mixing_beta = 0.3;
    params.mixing_history = 5;    
    params.tol_energy = 1e-6;
    params.max_iter = 100;

    run_scf_loop(&mut sim, params);

    println!("\nSimulação SCF Concluída.");
    Ok(())
}

fn main() {
    if let Err(e) = run_scf_full() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}