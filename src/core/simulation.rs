use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;
use ndarray::Array3;

// Imports dos seus módulos
use crate::core::kpoints::KGrid;
use crate::core::structure::Structure;
use crate::io::upf::{Pseudopotential, UpfError};
use crate::utils::welcome::print_welcome;
use crate::core::basis::PlaneWaveBasis; // Novo
use crate::core::fft::FftGrid;         // Novo

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
}

pub struct SimulationBuilder {
    structure: Option<Structure>,
    ecut: Option<f64>,
    k_grid: Option<KGrid>,
}

impl SimulationBuilder {
    pub fn new() -> Self {
        Self {
            structure: None,
            ecut: None,
            k_grid: None,
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

        Ok(Simulation {
            structure,
            ecut,
            k_grid,
            pseudos,
            bases,
            fft_grid,
            rho,
        })
    }
}