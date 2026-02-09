use std::process;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::Simulation;
use bravie::dft::xc::calculate_xc_lda; // Importe o novo módulo

fn run_xc_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Teste de Integração: LDA Exchange-Correlation ===\n");

    // 1. Setup (Copia do anterior)
    let a = 10.26;
    let si = Species { id: 0, element: "Si".to_string(), atomic_number: 14, mass: 28.085, pseudo_path: "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF".to_string() };
    let silicon = Structure::builder()
        .lattice([0.0, a/2.0, a/2.0], [a/2.0, 0.0, a/2.0], [a/2.0, a/2.0, 0.0])
        .add_species(si).add_atom([0.0, 0.0, 0.0], 0).add_atom([0.25, 0.25, 0.25], 0)
        .build()?;

    let mut sim = Simulation::builder()
        .structure(silicon).ecut(30.0).k_grid(KGrid::gamma()).build()?;
    
    sim.initialize_density();

    // 2. Calcula XC
    println!("Calculando XC (LDA-PZ81)...");
    let (v_xc, e_xc) = calculate_xc_lda(&sim.rho, sim.structure.lattice.volume());

    // 3. Resultados
    println!("\n=== RESULTADOS XC ===");
    println!("Energia XC Total:   {:.6} Ry", e_xc);
    println!("Potencial Médio:    {:.6} Ry", v_xc.mean().unwrap());
    println!("Potencial Mínimo:   {:.6} Ry", v_xc.iter().fold(f64::INFINITY, |a, &b| a.min(b)));
    
    if e_xc > 0.0 {
        return Err("ERRO: Energia XC deve ser negativa (ligação).".into());
    }

    Ok(())
}

fn main() {
    if let Err(e) = run_xc_test() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}