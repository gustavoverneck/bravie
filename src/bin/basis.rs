use std::process;
use ndarray::Array1;
use num_complex::{Complex, Complex64};

// Imports do Bravie
use bravie::core::structure::{Structure, Species};
use bravie::core::kpoints::KGrid;
use bravie::utils::welcome::print_welcome;
use bravie::Simulation;

fn run_basis_demo() -> Result<(), Box<dyn std::error::Error>> {
    print_welcome();
    println!("=== Demonstração: Basis Set e FFT ===\n");

    // 1. Configurar Estrutura (Exemplo: Silício FCC)
    let a = 10.26; // Parâmetro de rede (Bohr)
    
    // Espécie dummy (o pseudopotencial precisa existir ou o mock deve aceitar)
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

    // 2. Construir a Simulação
    // Isso dispara internamente a criação das Bases e do Grid FFT
    let ecut = 20.0; // Ry
    println!("Construindo simulação com Ecut = {:.1} Ry...", ecut);

    let mut sim = Simulation::builder()
        .structure(silicon)
        .ecut(ecut)
        .k_grid(KGrid::gamma()) // Usando apenas ponto Gamma para o exemplo
        .build()?;

    // 3. Inspecionar os Motores Numéricos
    // Acessamos a primeira base (k=0) e o grid FFT
    let basis = &sim.bases[0]; 
    let fft = &mut sim.fft_grid;

    println!("\n[Análise da Base]");
    println!("  -> K-Point: {:?}", basis.k_point);
    println!("  -> Número de Ondas Planas (NPW): {}", basis.g_vectors.len());
    println!("  -> Cutoff Função de Onda: {:.1} Ry", basis.ecut);
    println!("  -> Cutoff Densidade (4x): {:.1} Ry", basis.ecut_rho);

    println!("\n[Análise do Grid FFT]");
    println!("  -> Dimensões: {} x {} x {}", fft.size[0], fft.size[1], fft.size[2]);
    println!("  -> Total de Pontos Reais: {}", fft.size[0] * fft.size[1] * fft.size[2]);

    // 4. Demonstração Prática: Ciclo FFT
    println!("\n[Teste Prático: Ida e Volta]");
    println!("Criando uma 'função de onda' sintética no espaço recíproco...");

    let n_pw = basis.g_vectors.len();
    let mut c_in = Array1::<Complex64>::zeros(n_pw);

    // Preenche com valores arbitrários (ex: 1/|G|)
    for i in 0..n_pw {
        c_in[i] = Complex::new(1.0 / (i as f64 + 1.0), 0.0);
    }

    // A. Recíproco -> Real (Transformada Inversa)
    // Isso popula o buffer interno 'sim.fft_grid.buffer'
    fft.to_real_space(&c_in);
    
    // Vamos espiar o valor em um ponto do espaço real
    let val_real = fft.buffer[[0, 0, 0]];
    println!("  -> Valor no ponto r=(0,0,0): {:.4e} + {:.4e}i", val_real.re, val_real.im);

    // B. Real -> Recíproco (Transformada Direta)
    let mut c_out = Array1::<Complex64>::zeros(n_pw);
    fft.to_recip_space(&mut c_out);

    // C. Comparar Entrada vs Saída
    let mut max_diff = 0.0;
    for i in 0..n_pw {
        let diff = (c_in[i] - c_out[i]).norm();
        if diff > max_diff {
            max_diff = diff;
        }
    }

    println!("  -> Erro máximo de reconstrução: {:.4e}", max_diff);

    if max_diff < 1e-9 {
        println!("\nSUCESSO: O motor numérico está consistente!");
    } else {
        println!("\nALERTA: Diferença numérica detectada.");
    }

    Ok(())
}

fn main() {
    if let Err(e) = run_basis_demo() {
        eprintln!("Erro: {}", e);
        process::exit(1);
    }
}