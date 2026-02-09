use std::process;
use bravie::Simulation;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::dft::scf::{run_scf_loop, ScfParameters};

fn run_ecut_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Validação Física: Convergência de Ecut ===\n");

    let a = 10.26; // Fixo no experimental
    let si_species = Species {
        id: 0, element: "Si".to_string(), atomic_number: 14, mass: 28.085,
        pseudo_path: "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF".to_string(),
    };

    let ecuts = vec![15.0, 20.0, 30.0, 40.0, 50.0];

    println!("Ecut (Ry) | Energia Total (Ry) | Diferença");
    println!("----------|--------------------|----------");

    let mut prev_energy = 0.0;

    for &ecut in &ecuts {
        let silicon = Structure::builder()
            .lattice([0.0, a/2.0, a/2.0], [a/2.0, 0.0, a/2.0], [a/2.0, a/2.0, 0.0])
            .add_species(si_species.clone())
            .add_atom([0.0, 0.0, 0.0], 0)
            .add_atom([0.25, 0.25, 0.25], 0)
            .build()?;

        let mut sim = Simulation::builder()
            .structure(silicon)
            .ecut(ecut)
            .k_grid(KGrid::gamma())
            .build()?;
        
        sim.initialize_density();

        let mut params = ScfParameters::default();
        params.tol_energy = 1e-5;
        
        let energy = run_scf_loop(&mut sim, params);
        
        let delta = if prev_energy == 0.0 { 0.0 } else { energy - prev_energy };
        println!("RESULTADO_ECUT: {:.1} , {:.8} , {:.1e}", ecut, energy, delta);
        
        prev_energy = energy;
    }

    Ok(())
}

fn main() {
    if let Err(e) = run_ecut_test() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}