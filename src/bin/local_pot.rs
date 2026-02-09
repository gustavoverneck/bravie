use std::process;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::Simulation;

// Importe o novo módulo
use bravie::dft::local_potential::calculate_local_potential;

fn run_local_pot_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Teste de Integração: Potencial Local (V_loc) ===\n");

    // 1. Configurar Estrutura (Silício FCC)
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

    // 2. Construir Simulação
    let ecut = 30.0;
    println!("Construindo simulação (Ecut={:.1} Ry)...", ecut);
    let mut sim = Simulation::builder()
        .structure(silicon)
        .ecut(ecut)
        .k_grid(KGrid::gamma())
        .build()?;

    // 3. Calcular V_loc
    println!("\nCalculando Potencial Local...");
    let v_local = calculate_local_potential(
        &sim.structure,
        &sim.fft_grid,
        &sim.pseudos
    );

    // 4. Análise Estatística
    let v_max = v_local.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let v_min = v_local.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let v_mean = v_local.mean().unwrap();

    println!("\n=== RESULTADOS V_LOC ===");
    println!("Potencial Mínimo (Perto do núcleo): {:.6} Ry", v_min);
    println!("Potencial Máximo (Interstício):     {:.6} Ry", v_max);
    println!("Potencial Médio:                    {:.6} Ry", v_mean);

    // 5. Validação
    if v_min > 0.0 {
        return Err("ERRO: V_loc deve ser atrativo (negativo) perto dos núcleos!".into());
    }
    
    // Check de ponto específico (Origem = Posição do Átomo)
    let v_origin = v_local[[0, 0, 0]];
    println!("Valor na Origem (0,0,0):            {:.6} Ry", v_origin);

    println!("\n✅ SUCESSO: Potencial Local calculado.");
    Ok(())
}

fn main() {
    if let Err(e) = run_local_pot_test() {
        eprintln!("Erro Fatal: {}", e);
        process::exit(1);
    }
}