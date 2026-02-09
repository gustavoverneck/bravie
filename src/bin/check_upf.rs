use std::env;
use std::path::Path;
use bravie::io::upf::Pseudopotential;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Uso: cargo run --bin check_upf <caminho_para_arquivo.upf>");
        return Ok(());
    }

    let path = Path::new(&args[1]);
    let pseudo = Pseudopotential::from_file(path)?;

    println!("--- Análise do Pseudopotencial: {} ---", pseudo.header.element);
    println!("Carga de Valência (Z_valence): {:.4}", pseudo.header.z_valence);
    println!("Tamanho da Malha Radial: {}", pseudo.header.mesh_size);

    // Integração Numérica Radial
    // Integral = sum( rho_atom[i] * rab[i] )
    // Onde 'rab' é o peso de integração radial fornecido pelo UPF (dr/dx)
    let mut integral_charge = 0.0;
    
    // Verifica tamanhos para evitar pânico
    let n_points = pseudo.mesh.r.len().min(pseudo.rho_atom.len()).min(pseudo.mesh.rab.len());

    for i in 0..n_points {
        // Regra do Trapézio simples usando os pesos rab
        // UPF define: int f(r) dr = sum_i f(r_i) * rab(i)
        if i < n_points {
             integral_charge += pseudo.rho_atom[i].abs() * pseudo.mesh.rab[i];
        }
    }

    println!("Integral Radial Calculada (Soma direta): {:.4}", integral_charge);
    
    let diff = (integral_charge - pseudo.header.z_valence).abs();
    println!("Erro Absoluto: {:.4}", diff);
    
    if diff < 0.1 {
        println!("A integral radial do UPF está consistente.");
    } else {
        println!("AVISO: A integral radial difere significativamente de Z_valence.");
        println!("Isso sugere que 'rho_atom' no UPF pode ter uma definição diferente");
        println!("(ex: densidade de core incluída ou normalização diferente).");
    }

    Ok(())
}