use std::f64::consts::PI;
use ndarray::Array3;
use nalgebra::Vector3;
use crate::core::structure::Structure;
use crate::core::fft::FftGrid;
use crate::io::upf::Pseudopotential;
use std::collections::HashMap;

/// Calcula a densidade inicial (SAD) e aplica renormalização de carga.
pub fn calculate_initial_density(
    structure: &Structure,
    fft_grid: &FftGrid,
    pseudos: &HashMap<usize, Pseudopotential>
) -> Array3<f64> {
    let (nx, ny, nz) = (fft_grid.size[0], fft_grid.size[1], fft_grid.size[2]);
    let mut rho = Array3::<f64>::zeros((nx, ny, nz));

    // Matriz inversa para condições de contorno periódicas
    // Usa .vectors conforme sua estrutura atual
    let lattice_inv = structure.lattice.vectors.try_inverse().expect("Lattice matrix singular");

    // 1. Superposição das Densidades Atômicas
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // Posição fracionária [0, 1]
                let frac_pos = Vector3::new(
                    i as f64 / nx as f64,
                    j as f64 / ny as f64,
                    k as f64 / nz as f64,
                );
                
                // Posição cartesiana real
                let r_grid = structure.lattice.vectors * frac_pos;
                let mut rho_val = 0.0;

                for atom in &structure.atoms {
                    let pseudo = pseudos.get(&atom.species_id)
                        .expect("Pseudopotencial não encontrado");

                    // Vetor diferença
                    let diff = r_grid - atom.position;

                    // Minimum Image Convention (MIC)
                    let mut d_frac = lattice_inv * diff;
                    d_frac.x -= d_frac.x.round();
                    d_frac.y -= d_frac.y.round();
                    d_frac.z -= d_frac.z.round();

                    let d_cart = structure.lattice.vectors * d_frac;
                    let dist = d_cart.norm();

                    // Interpola valor do UPF
                    rho_val += interpolate_rho_atom(dist, pseudo);
                }

                rho[[i, j, k]] = rho_val;
            }
        }
    }

    // 2. Renormalização de Carga (CRÍTICO)
    // Corrige erros de amostragem do grid garantindo neutralidade exata.
    
    let volume = structure.lattice.volume();
    let n_points = (nx * ny * nz) as f64;
    let dvol = volume / n_points; // Volume por voxel
    
    // Carga integrada numericamente (geralmente menor que o esperado)
    let current_charge: f64 = rho.sum() * dvol;
    
    // Carga total esperada (soma dos Z_valence dos átomos)
    let mut target_charge = 0.0;
    for atom in &structure.atoms {
        if let Some(p) = pseudos.get(&atom.species_id) {
            target_charge += p.header.z_valence;
        }
    }

    println!("   > Carga SAD Calculada: {:.6} e", current_charge);
    println!("   > Carga Alvo (Z_val) : {:.6} e", target_charge);

    if current_charge.abs() > 1e-9 {
        let scale = target_charge / current_charge;
        println!("   > Aplicando Fator de Renormalização: {:.6}", scale);
        
        // Multiplica todo o grid pelo fator de correção
        // rho *= scale (ndarray suporta ops escalares)
        rho.mapv_inplace(|v| v * scale);
    } else {
        println!("   > AVISO: Carga zero detectada, pulando renormalização.");
    }

    rho
}

fn interpolate_rho_atom(r: f64, pseudo: &Pseudopotential) -> f64 {
    let mesh = &pseudo.mesh;
    let rho_data = &pseudo.rho_atom; // Lembre-se: UPF armazena 4*pi*r^2 * rho

    if r > mesh.r.last().copied().unwrap_or(0.0) {
        return 0.0;
    }

    // Tratamento de singularidade r -> 0
    // Densidade deve ser finita, mas UPF guarda r^2 * rho.
    // Usamos o segundo ponto para evitar divisão por zero
    if r < 1e-6 {
        if mesh.r.len() > 1 {
             let r_safe = mesh.r[1];
             let val = rho_data[1];
             return val / (4.0 * PI * r_safe * r_safe);
        }
        return 0.0;
    }

    // Busca Binária
    let idx = match mesh.r.binary_search_by(|val| val.partial_cmp(&r).unwrap()) {
        Ok(i) => i,
        Err(i) => if i > 0 { i - 1 } else { 0 },
    };

    if idx >= mesh.r.len() - 1 {
        return 0.0;
    }

    // Interpolação Linear
    let r1 = mesh.r[idx];
    let r2 = mesh.r[idx+1];
    let y1 = rho_data[idx];
    let y2 = rho_data[idx+1];

    let t = (r - r1) / (r2 - r1);
    let y_interp = y1 + t * (y2 - y1);

    // Converte de Radial Charge (UPF) para Volumetric Charge
    y_interp / (4.0 * PI * r * r)
}