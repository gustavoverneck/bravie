use nalgebra::Vector3;
use crate::core::structure::Structure;

/// Representa a base de ondas planas para um ponto K específico.
/// Responsável por determinar a geometria do grid e listar os vetores G ativos.
pub struct PlaneWaveBasis {
    /// Energia de corte para as funções de onda (Ry)
    pub ecut: f64,
    
    /// Energia de corte para a densidade de carga (geralmente 4 * ecut)
    pub ecut_rho: f64,
    
    /// Dimensões do grid FFT (nx, ny, nz)
    pub fft_grid: [usize; 3],
    
    /// Lista de índices (i, j, k) dos vetores G onde |k + G|^2 <= Ecut
    pub g_vectors: Vec<(i32, i32, i32)>,
    
    /// Ponto K associado a esta base (coordenadas fracionárias)
    pub k_point: Vector3<f64>,
}

impl PlaneWaveBasis {
    /// Cria uma nova base para um dado Structure e Ecut.
    /// Se k_point for None, assume Gamma (0, 0, 0).
    pub fn new(structure: &Structure, ecut: f64, k_point: Option<[f64; 3]>) -> Self {
        // 1. Definição do Dual Grid (Densidade requer 4x a energia)
        // Isso garante que a convolução |psi|^2 seja exata no grid.
        let ecut_rho = 4.0 * ecut;
        
        let k_vec = if let Some(k) = k_point {
            Vector3::from(k)
        } else {
            Vector3::zeros()
        };

        // 2. Calcula as dimensões ótimas do grid FFT (baseado em Ecut_rho)
        let fft_grid = Self::calculate_optimal_fft_grid(&structure.lattice.reciprocal(), ecut_rho);

        // 3. Gera os vetores G ativos para este k-point (baseado em Ecut)
        let g_vectors = Self::generate_g_vectors(structure, fft_grid, ecut, k_vec);

        println!(
            "    Basis Init: Ecut={:.1} Ry | Grid=[{}, {}, {}] | NG={} (k={:?})",
            ecut, fft_grid[0], fft_grid[1], fft_grid[2], g_vectors.len(), k_vec.as_slice()
        );

        Self {
            ecut,
            ecut_rho,
            fft_grid,
            g_vectors,
            k_point: k_vec,
        }
    }

    /// Calcula tamanho do grid para evitar aliasing (Shannon-Nyquist).
    /// Grid deve cobrir 2 * G_max_rho.
    fn calculate_optimal_fft_grid(recip_lattice: &nalgebra::Matrix3<f64>, ecut_rho: f64) -> [usize; 3] {
        // G_max é o raio da esfera de densidade no espaço recíproco
        let g_max = ecut_rho.sqrt();

        // Normas dos vetores da rede recíproca (colunas da matriz transposta ou linhas da direta?)
        // Lattice::reciprocal() retorna B onde B^T * A = 2pi * I. 
        // Assumindo que recip_lattice colunas são b1, b2, b3.
        let b1 = recip_lattice.column(0).norm();
        let b2 = recip_lattice.column(1).norm();
        let b3 = recip_lattice.column(2).norm();

        // Número de pontos necessários: 2 * g_max / |bi|
        // O fator 2 vem do teorema de Nyquist (amostragem)
        let nx = (2.0 * g_max / b1).ceil() as usize;
        let ny = (2.0 * g_max / b2).ceil() as usize;
        let nz = (2.0 * g_max / b3).ceil() as usize;

        [
            Self::next_fft_size(nx),
            Self::next_fft_size(ny),
            Self::next_fft_size(nz),
        ]
    }

    /// Gera a lista de vetores G inteiros (i, j, k) dentro da esfera de energia cinética.
    fn generate_g_vectors(
        structure: &Structure, 
        grid_dim: [usize; 3], 
        ecut: f64, 
        k_point: Vector3<f64>
    ) -> Vec<(i32, i32, i32)> {
        let mut g_vecs = Vec::new();
        let recip = structure.lattice.reciprocal(); // Matriz onde colunas são b1, b2, b3

        // Define a caixa de busca baseada no grid FFT.
        // O grid FFT é calculado para 4*Ecut, então cobrir metade dele (frequência de Nyquist)
        // é mais que suficiente para encontrar todos os vetores de Ecut.
        // Usamos índices com sinal (i32).
        let nx_search = (grid_dim[0] / 2) as i32;
        let ny_search = (grid_dim[1] / 2) as i32;
        let nz_search = (grid_dim[2] / 2) as i32;

        for i in -nx_search..=nx_search {
            for j in -ny_search..=ny_search {
                for k in -nz_search..=nz_search {
                    let g_int = Vector3::new(i as f64, j as f64, k as f64);
                    
                    // Vetor K + G em coordenadas fracionárias
                    let kg_frac = k_point + g_int;

                    // Vetor Cartesiano: q = B * (k+G)
                    let q_cart = recip * kg_frac;

                    // Energia cinética (unidades atômicas/Ry) = |q|^2
                    if q_cart.norm_squared() <= ecut {
                        g_vecs.push((i, j, k));
                    }
                }
            }
        }
        
        // Otimização opcional: Ordenar vetores por energia (ajuda no pré-condicionamento depois)
        // g_vecs.sort_by(|a, b| ...);
        
        g_vecs
    }

    /// Encontra o próximo tamanho de grid que é produto de primos pequenos (2, 3, 5, 7).
    /// Isso é crítico para performance O(N log N) da FFT.
    fn next_fft_size(n: usize) -> usize {
        let mut size = n;
        while !Self::is_smooth_number(size) {
            size += 1;
        }
        size
    }

    fn is_smooth_number(mut n: usize) -> bool {
        if n == 0 { return false; }
        for p in &[2, 3, 5, 7] {
            while n % p == 0 {
                n /= p;
            }
        }
        n == 1
    }
}