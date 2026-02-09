use nalgebra::Vector3;

use crate::core::structure::Structure;

#[derive(Debug, Clone)]
pub struct KPoint {
    pub coord: [f64; 3], // Coordenadas fracionárias (em relação à recíproca)
    pub weight: f64,     // Peso para integração na Zona de Brillouin
}

#[derive(Debug, Clone)]
pub struct KGrid {
    pub k_points: Vec<KPoint>,
    // Futuro: Simetrias para reduzir o número de pontos
}

impl KGrid {
    /// Cria um grid contendo apenas o ponto Gamma [0,0,0]
    pub fn gamma() -> Self {
        Self {
            k_points: vec![KPoint {
                coord: [0.0, 0.0, 0.0],
                weight: 1.0,
            }],
        }
    }

    /// Gera uma malha Monkhorst-Pack (uniforme)
    /// `grid`: Número de pontos em cada direção [nkx, nky, nkz]
    /// `shift`: Deslocamento [sx, sy, sz] (geralmente 0.0 ou 0.5)
    pub fn monkhorst_pack(grid: [usize; 3], shift: [f64; 3]) -> Self {
        let mut k_points = Vec::new();
        let nk = grid[0] * grid[1] * grid[2];
        let weight = 1.0 / nk as f64;

        for i in 0..grid[0] {
            for j in 0..grid[1] {
                for k in 0..grid[2] {
                    // Fórmula MP: u_r = (2r - n + 1) / 2n  (centrado)
                    // Ou simples: (i + shift) / n
                    let kx = (i as f64 + shift[0]) / grid[0] as f64;
                    let ky = (j as f64 + shift[1]) / grid[1] as f64;
                    let kz = (k as f64 + shift[2]) / grid[2] as f64;

                    // Normaliza para o intervalo [-0.5, 0.5) (Primeira BZ)
                    let wrap = |x: f64| if x >= 0.5 { x - 1.0 } else { x };

                    k_points.push(KPoint {
                        coord: [wrap(kx), wrap(ky), wrap(kz)],
                        weight,
                    });
                }
            }
        }

        Self { k_points }
    }

    pub fn band_path(points: Vec<[f64; 3]>, points_per_segment: usize) -> Self {
        let mut k_points = Vec::new();
        let weight = 0.0; // Bandas não têm peso no cálculo de densidade (só geometria)

        // Itera sobre pares de pontos: (P0 -> P1), (P1 -> P2), etc.
        for i in 0..(points.len() - 1) {
            let start = Vector3::from(points[i]);
            let end = Vector3::from(points[i+1]);
            
            let vector = end - start; // Vetor direção
            
            // Adiciona os pontos intermediários
            for step in 0..points_per_segment {
                let t = step as f64 / points_per_segment as f64; // Fração de 0.0 a 1.0
                let k_vec = start + vector * t;
                
                k_points.push(KPoint {
                    coord: [k_vec.x, k_vec.y, k_vec.z],
                    weight, 
                });
            }
        }
        
        // Adiciona o último ponto final para fechar o caminho
        k_points.push(KPoint {
            coord: points.last().unwrap().clone(),
            weight,
        });

        Self { k_points }
    }
}