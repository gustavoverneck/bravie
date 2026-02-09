use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;
use crate::core::structure::Structure;
use crate::io::upf::{Pseudopotential, UpfError}; //
use crate::utils::welcome::print_welcome;


pub struct Simulation {
    pub structure: Structure,                 // A rede cristalina e átomos
    pub ecut: f64,                            // Energia de corte (Kinetic Energy Cutoff) em Ry ou Ha
    pub pseudos: HashMap<usize, Pseudopotential>, // Mapa: ID da Espécie -> Pseudopotencial carregado
}

#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("A estrutura cristalina não foi definida.")]
    MissingStructure,

    #[error("A energia de corte (Ecut) não foi definida.")]
    MissingEcut,

    #[error("Arquivo de pseudopotencial não encontrado para espécie '{0}': {1}")]
    PseudoFileNotFound(String, String),

    #[error("Erro ao carregar pseudopotencial: {0}")]
    UpfLoadError(#[from] UpfError), // Encapsula automaticamente erros do UPF
}

impl Simulation {
    pub fn builder() -> SimulationBuilder {
        SimulationBuilder::new()
    }
    
    /// Inicia o ciclo de cálculo
    pub fn run(&self) {
        print_welcome();
        println!("Iniciando simulação para {} átomos.", self.structure.atoms.len());
        println!("Energia de corte: {:.2}", self.ecut);
        // Lógica do SCF viria aqui...
    }
}


pub struct SimulationBuilder {
    structure: Option<Structure>,
    ecut: Option<f64>,
}

impl SimulationBuilder {
    /// Cria um novo builder vazio
    pub fn new() -> Self {
        Self {
            structure: None,
            ecut: None, // Adicionar valor padrão?
        }
    }

    /// Define a estrutura cristalina
    pub fn structure(mut self, structure: Structure) -> Self {
        self.structure = Some(structure);
        self
    }

    /// Define a energia de corte (ecut)
    pub fn ecut(mut self, ecut: f64) -> Self {
        self.ecut = Some(ecut);
        self
    }

    pub fn build(self) -> Result<Simulation, SimulationError> {
        // Validações básicas
        let structure = self.structure.ok_or(SimulationError::MissingStructure)?;
        let ecut = self.ecut.ok_or(SimulationError::MissingEcut)?;

        // Carregamento automático dos Pseudos
        let mut pseudos = HashMap::new();
        
        println!("Verificando pseudopotenciais...");
        
        for species in &structure.species {
            let path_str = &species.pseudo_path;
            let path = Path::new(path_str);

            // Verifica existência do arquivo antes de tentar abrir
            if !path.exists() {
                return Err(SimulationError::PseudoFileNotFound(
                    species.element.clone(),
                    path_str.clone()
                ));
            }

            // Tenta carregar e converte UpfError -> SimulationError automaticamente graças ao #[from]
            let upf = Pseudopotential::from_file(path)?;
            
            pseudos.insert(species.id, upf);
            println!("   OK {} -> {}", species.element, path_str);
        }

        Ok(Simulation {
            structure,
            ecut,
            pseudos,
        })
    }

}