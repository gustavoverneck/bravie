use nalgebra::{Matrix3, Vector3};
use std::fmt;

#[derive(Debug, Clone)]
pub struct Species {
    pub id: usize,          // Índice único
    pub element: String,    // Símbolo "H", "C", "Si"
    pub atomic_number: u8,  // Z
    pub mass: f64,          // Massa atômica (opcional por enquanto)
    pub pseudo_path: String,// Caminho para arquivo .upf
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub species_id: usize,
    pub position: Vector3<f64>,
}

#[derive(Debug, Clone)]
pub struct Lattice {
    pub vectors: Matrix3<f64>,
}

impl Lattice {
    pub fn new(a1: Vector3<f64>,a2: Vector3<f64>, a3: Vector3<f64>) -> Self {
        Self {
            vectors: Matrix3::from_columns(&[a1, a2, a3])
        }
    }

    pub fn volume(&self) -> f64 {
        self.vectors.determinant().abs()
    }

    pub fn reciprocal(&self) -> Matrix3<f64> {
        let vol = self.volume();
        let a1 = self.vectors.column(0);
        let a2 = self.vectors.column(0);
        let a3 = self.vectors.column(0);
        let factor = 2.0 * std::f64::consts::PI / vol;

        let b1 = a2.cross(&a3) * factor;
        let b2 = a3.cross(&a1) * factor;
        let b3 = a1.cross(&a2) * factor;

        Matrix3::from_columns(&[b1, b2, b3])
    }

}
#[derive(Debug, Clone)]
pub struct Structure {
    pub lattice: Lattice,
    pub species: Vec<Species>,
    pub atoms: Vec<Atom>,
}

pub struct StructureBuilder {
    pub lattice: Option<Lattice>,
    pub species: Vec<Species>,
    pub atoms: Vec<Atom>,   
}

impl StructureBuilder {
    pub fn new() -> Self{
        Self { 
            lattice: None, 
            species: Vec::new(),
            atoms: Vec::new() 
        }
    }

    pub fn fcc(mut self, a: f64) -> Self {
        let half = a / 2.0;
        self.lattice = Some(Lattice::new(
            Vector3::new(half, 0.0, 0.0),
            Vector3::new(0.0, half, 0.0),
            Vector3::new(0.0, 0.0, half),
        ));
        self
    }

    pub fn cubic(mut self, a: f64) -> Self {
        self.lattice = Some(Lattice::new(
            Vector3::new(a, 0.0, 0.0),
            Vector3::new(0.0, a, 0.0),
            Vector3::new(0.0, 0.0, a),
        ));
        self
    }

    pub fn lattice(mut self, a1: [f64; 3], a2: [f64; 3], a3: [f64; 3]) -> Self {
        let vectors = Matrix3::from_columns(&[
            Vector3::new(a1[0], a1[1], a1[2]),
            Vector3::new(a2[0], a2[1], a2[2]),
            Vector3::new(a3[0], a3[1], a3[2]),
        ]);

        self.lattice = Some(Lattice { vectors });
        self
    }

    pub fn add_atom(mut self, pos: [f64; 3], species_id: usize) -> Self {
        self.atoms.push(Atom {
        position: Vector3::from(pos),
        species_id,
        });
        self
    }

    pub fn add_species(mut self, species: Species) -> Self {
        self.species.push(species);
        self
    }


    pub fn build(self) -> Result<Structure, String> {
        let lattice = self.lattice.ok_or("Erro: Rede não definida!")?;
        if self.atoms.is_empty() {
            return Err("Erro: Estrutura vazia.".to_string());
        }
        Ok(Structure { 
            lattice, 
            species: self.species,
            atoms: self.atoms 
        })
    }
}

impl Structure {
    pub fn builder() -> StructureBuilder {
        StructureBuilder::new()
    }
}

impl fmt::Display for Structure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Structure (Volume: {:.4} a.u.^3)", self.lattice.volume())?;
        writeln!(f, "Lattice Vectors:")?;
        
        // Acessa as colunas e formata manualmente
        let a1 = self.lattice.vectors.column(0);
        let a2 = self.lattice.vectors.column(1);
        let a3 = self.lattice.vectors.column(2);
        
        writeln!(f, "  a1: [{:.4}, {:.4}, {:.4}]", a1[0], a1[1], a1[2])?;
        writeln!(f, "  a2: [{:.4}, {:.4}, {:.4}]", a2[0], a2[1], a2[2])?;
        writeln!(f, "  a3: [{:.4}, {:.4}, {:.4}]", a3[0], a3[1], a3[2])?;
        
        writeln!(f, "Atoms ({}):", self.atoms.len())?;
        for (i, atom) in self.atoms.iter().enumerate() {
            let element = self.species.iter()
                .find(|s| s.id == atom.species_id)
                .map(|s| s.element.as_str())
                .unwrap_or("Unknown");
            writeln!(f, "  {}: {} at [{:.4}, {:.4}, {:.4}]", 
                i+1, element, atom.position[0], atom.position[1], atom.position[2])?;
        }
        Ok(())
    }
}