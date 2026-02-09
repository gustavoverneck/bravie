// src/core/simulation.rs

use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;
use ndarray::Array3;

// Imports dos seus módulos
use crate::core::kpoints::KGrid;
use crate::core::structure::Structure;
use crate::io::upf::{Pseudopotential, UpfError};
use crate::utils::welcome::print_welcome;
use crate::core::basis::PlaneWaveBasis;
use crate::core::fft::FftGrid;         
use crate::dft::density::calculate_initial_density;

#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("A estrutura cristalina não foi definida.")]
    MissingStructure,

    #[error("A energia de corte (Ecut) não foi definida.")]
    MissingEcut,
    
    #[error("K-Grid vazio ou inválido.")]
    InvalidKGrid,

    #[error("Arquivo de pseudopotencial não encontrado para espécie '{0}': {1}")]
    PseudoFileNotFound(String, String),

    #[error("Erro ao carregar pseudopotencial: {0}")]
    UpfLoadError(#[from] UpfError),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HamiltonianModel {
    Schrodinger,            // Clássico (p^2 / 2m)
    DiracScalarRelativistic // Relativístico (sqrt(p^2 c^2 + m^2 c^4))
}

pub struct Simulation {
    // Inputs Físicos
    pub structure: Structure,
    pub ecut: f64,
    pub k_grid: KGrid,
    pub pseudos: HashMap<usize, Pseudopotential>,

    // Motores de Cálculo (Adicionados)
    pub bases: Vec<PlaneWaveBasis>, // Bases de ondas planas (uma por k-point)
    pub fft_grid: FftGrid,          // Gerenciador da FFT e memória
    pub rho: Array3<f64>,           // Densidade de carga no espaço real
    pub v_eff: Array3<f64>,         // Potencial Efetivo Total no espaço real
    pub hamiltonian_model: HamiltonianModel,
}

impl Simulation {
    pub fn builder() -> SimulationBuilder {
        SimulationBuilder::new()
    }
    
    pub fn run(&mut self) {
        print_welcome();
        
        let natoms = self.structure.atoms.len();
        let nk = self.k_grid.k_points.len();
        let (nx, ny, nz) = (self.fft_grid.size[0], self.fft_grid.size[1], self.fft_grid.size[2]);

        println!("--- Inicialização Completa ---");
        println!("Sistema: {} átomos, {} espécies", natoms, self.pseudos.len());
        println!("K-Points: {} pontos na ZB", nk);
        println!("Grid FFT: {} x {} x {} (Total: {})", nx, ny, nz, nx*ny*nz);
        println!("Cutoffs: WFC={:.1} Ry, Rho={:.1} Ry", self.ecut, self.bases[0].ecut_rho);
        print!( "{}", self.structure.clone());
        
        // Aqui começaria o loop SCF:
        // 1. Inicializar densidade aleatória ou superposição atômica
        // 2. Loop { V_eff -> Diagonalização -> Rho_new -> Mix -> Check Convergência }
    }

    /// Preenche o grid rho com a superposição das densidades atômicas
    pub fn initialize_density(&mut self) {
        println!("Calculando densidade inicial (SAD)...");
        
        let rho_sad = calculate_initial_density(
            &self.structure, 
            &self.fft_grid, 
            &self.pseudos
        );
        
        // Atualiza o estado da simulação
        self.rho = rho_sad;
        
        // Check de Carga Total (Integral)
        // Carga = sum(rho) * volume_voxel
        // volume_voxel = det(Lattice) / N_grid
        let volume = self.structure.lattice.volume();
        let n_grid = self.fft_grid.size[0] * self.fft_grid.size[1] * self.fft_grid.size[2];
        let dvol = volume / n_grid as f64;
        
        let total_charge: f64 = self.rho.sum() * dvol;
        
        println!("Densidade inicial calculada.");
        println!("  - Carga Total Integrada: {:.4} e", total_charge);
        
        // Verifica neutralidade (soma dos eletrons de valencia)
        let mut expected_charge = 0.0;
        for atom in &self.structure.atoms {
            if let Some(p) = self.pseudos.get(&atom.species_id) {
                expected_charge += p.header.z_valence;
            }
        }
        println!("  - Carga Esperada (Zval): {:.4} e", expected_charge);
    }
}

pub struct SimulationBuilder {
    structure: Option<Structure>,
    ecut: Option<f64>,
    k_grid: Option<KGrid>,
    pub hamiltonian_model: HamiltonianModel,
}

impl SimulationBuilder {
    pub fn new() -> Self {
        Self {
            structure: None,
            ecut: None,
            k_grid: None,
            hamiltonian_model: HamiltonianModel::Schrodinger,
        }
    }

    pub fn structure(mut self, structure: Structure) -> Self {
        self.structure = Some(structure);
        self
    }

    pub fn ecut(mut self, ecut: f64) -> Self {
        self.ecut = Some(ecut);
        self
    }

    // Corrigido typo: k_drid -> k_grid
    pub fn k_grid(mut self, k_grid: KGrid) -> Self {
        self.k_grid = Some(k_grid);
        self
    }

    pub fn relativistic(mut self, enabled: bool) -> Self {
        self.hamiltonian_model = if enabled {
            HamiltonianModel::DiracScalarRelativistic
        } else {
            HamiltonianModel::Schrodinger
        };
        self
    }

    pub fn build(self) -> Result<Simulation, SimulationError> {
        // 1. Validações Básicas
        let structure = self.structure.ok_or(SimulationError::MissingStructure)?;
        let ecut = self.ecut.ok_or(SimulationError::MissingEcut)?;
        
        // Se K-Grid não for definido, assume Gamma Point
        let k_grid = self.k_grid.unwrap_or_else(|| KGrid::gamma());
        
        if k_grid.k_points.is_empty() {
            return Err(SimulationError::InvalidKGrid);
        }

        // 2. Carregamento de Pseudopotenciais
        let mut pseudos = HashMap::new();
        println!("Carregando pseudopotenciais...");
        
        for species in &structure.species {
            let path_str = &species.pseudo_path;
            let path = Path::new(path_str);

            if !path.exists() {
                return Err(SimulationError::PseudoFileNotFound(
                    species.element.clone(),
                    path_str.clone()
                ));
            }

            let upf = Pseudopotential::from_file(path)?;
            pseudos.insert(species.id, upf);
            println!("  [OK] {} -> {}", species.element, path_str);
        }

        // 3. Inicialização dos Motores Numéricos (Basis e FFT)
        println!("Inicializando grids e bases...");
        
        // Gera uma base de ondas planas para CADA ponto K
        // Precisamos acessar .coord do KPoint
        let bases: Vec<PlaneWaveBasis> = k_grid.k_points.iter()
            .map(|kp| {
                PlaneWaveBasis::new(&structure, ecut, Some(kp.coord))
            })
            .collect();

        // O Grid FFT é geométrico, independe do k-point (exceto para algoritmos avançados).
        // Usamos a primeira base para definir as dimensões (nx, ny, nz).
        let fft_grid = FftGrid::new(&bases[0]);

        // 4. Alocação da Densidade (Rho)
        let (nx, ny, nz) = (fft_grid.size[0], fft_grid.size[1], fft_grid.size[2]);
        let rho = Array3::<f64>::zeros((nx, ny, nz));
        let v_eff = Array3::<f64>::zeros((nx, ny, nz));

        Ok(Simulation {
            structure,
            ecut,
            k_grid,
            pseudos,
            bases,
            fft_grid,
            rho,
            v_eff,
            hamiltonian_model: self.hamiltonian_model,
        })
    }

}