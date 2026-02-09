use std::process;
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::Simulation;

// Importa o solver de Hartree. 
// Certifique-se de que "pub mod dft" está no lib.rs e "pub mod potentials" no dft/mod.rs
use bravie::dft::potentials::solve_hartree; 

fn run_potentials_test() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Teste de Integração: Potencial de Hartree ===\n");

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
    // Ecut maior melhora a precisão da densidade perto do núcleo
    let ecut = 30.0; 
    println!("Construindo simulação (Ecut={:.1} Ry)...", ecut);

    let mut sim = Simulation::builder()
        .structure(silicon)
        .ecut(ecut)
        .k_grid(KGrid::gamma())
        .build()?;

    // 3. Inicializar Densidade (SAD)
    // Isso popula sim.rho com a superposição atômica e corrige a carga
    sim.initialize_density();

    // 4. Calcular Potencial de Hartree
    println!("\nCalculando Potencial de Hartree...");
    
    // O solver precisa: 
    // - Rho (entrada)
    // - FFT Grid (buffer mutável para cálculo)
    // - Base (para obter |G|^2)
    // - Estrutura (para obter o Volume da célula)
    let (v_hartree, e_hartree) = solve_hartree(
        &sim.rho,
        &mut sim.fft_grid,
        &sim.bases[0],
        &sim.structure
    );

    // 5. Resultados e Análise
    println!("\n=== RESULTADOS HARTREE ===");
    println!("Energia de Hartree Total: {:.6} Ry", e_hartree);
    
    // Estatísticas do Potencial (V_H deve ser puramente repulsivo/positivo para elétrons?)
    // Nota: O sinal depende da convenção V_el-el. Geralmente é positivo (repulsão).
    let v_mean = v_hartree.mean().unwrap();
    
    // Encontrar min/max manualmente (ndarray iter)
    let v_max = v_hartree.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let v_min = v_hartree.iter().fold(f64::INFINITY, |a, &b| a.min(b));

    println!("Potencial Médio (V_H):    {:.6} Ry", v_mean);
    println!("Potencial Máximo:         {:.6} Ry", v_max);
    println!("Potencial Mínimo:         {:.6} Ry", v_min);

    // Validação Básica
    if e_hartree < 0.0 {
        return Err("ERRO: Energia de Hartree negativa! A auto-interação deve ser positiva.".into());
    }
    
    println!("\n✅ SUCESSO: Cálculo do potencial concluído sem erros.");

    Ok(())
}

fn main() {
    if let Err(e) = run_potentials_test() {
        eprintln!("Erro Fatal: {}", e);
        process::exit(1);
    }
}