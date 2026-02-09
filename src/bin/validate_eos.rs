use std::process;
use bravie::Simulation;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::dft::scf::{run_scf_loop, ScfParameters};

fn run_eos_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Validação Física: Equação de Estado (Energia vs Parâmetro de Rede) ===\n");

    // 1. Configuração Básica
    let si_pseudo = "pp/Si.pbe-n-kjpaw_psl.1.0.0.UPF";
    let si_species = Species {
        id: 0,
        element: "Si".to_string(),
        atomic_number: 14,
        mass: 28.085,
        pseudo_path: si_pseudo.to_string(),
    };

    // Parâmetros de teste
    let lattice_params = vec![
        9.8, 10.0, 10.1, 10.2, 10.26, 10.3, 10.4, 10.6, 10.8
    ]; // Bohr (10.26 é o experimental aproximado)
    
    println!("Param_Rede (Bohr) | Energia Total (Ry) | Convergência");
    println!("-------------------|--------------------|-------------");

    for &a in &lattice_params {
        // 2. Constrói estrutura com o parâmetro 'a' variável
        let silicon = Structure::builder()
            .lattice(
                [0.0, a/2.0, a/2.0],
                [a/2.0, 0.0, a/2.0],
                [a/2.0, a/2.0, 0.0]
            )
            .add_species(si_species.clone())
            .add_atom([0.0, 0.0, 0.0], 0)
            .add_atom([0.25, 0.25, 0.25], 0)
            .build()?;

        // 3. Simulação
        let mut sim = Simulation::builder()
            .structure(silicon)
            .ecut(30.0) // Ecut fixo
            .k_grid(KGrid::gamma())
            .build()?;
        
        sim.initialize_density();

        // 4. Roda SCF (Modo Silencioso/Rápido se possível)
        // Para não poluir o terminal, poderíamos silenciar o print do SCF,
        // mas vamos deixar para ver se algo explode.
        let mut params = ScfParameters::default();
        params.max_iter = 60;
        params.tol_energy = 1e-5;
        params.mixing_beta = 0.2; // Seguro
        params.mixing_history = 10;
        
        // Captura pânico se houver, para não parar o loop inteiro
        let energy = run_scf_loop(&mut sim, params);

        println!("RESULTADO_EOS: {:.4} , {:.8}", a, energy);
    }

    Ok(())
}

fn main() {
    if let Err(e) = run_eos_test() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}