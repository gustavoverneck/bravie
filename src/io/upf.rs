use roxmltree::Document;
use std::fs;
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum UpfError {
    #[error("Erro de Leitura de Arquivo: {0}")]
    Io(#[from] std::io::Error),
    #[error("Erro de Parsing XML: {0}")]
    Xml(#[from] roxmltree::Error),
    #[error("Campo obrigatório faltando no UPF: {0}")]
    MissingField(String),
    #[error("Erro ao converter string para número")]
    ParseNumber,
}

#[derive(Debug, Clone)]
pub struct Pseudopotential {
    pub header: Header,
    pub mesh: RadialMesh,
    pub local: Vec<f64>,        // Potencial Local V_loc(r)
    pub nonlocal: Vec<BetaFunction>, // Projetores Não-Locais Beta(r)
    pub rho_atom: Vec<f64>,     // Densidade Atômica (para chute inicial)
    pub dij: Vec<f64>,          // Matriz de coeficientes D_ij (Opcional)
}

#[derive(Debug, Clone)]
pub struct Header {
    pub element: String,
    pub z_valence: f64,
    pub mesh_size: usize,
    pub functional: String,
    pub number_of_proj: usize,
}

#[derive(Debug, Clone)]
pub struct RadialMesh {
    pub r: Vec<f64>,   // Grid Radial
    pub rab: Vec<f64>, // Derivada dr/dx (para integração)
}

#[derive(Debug, Clone)]
pub struct BetaFunction {
    pub index: usize,
    pub angular_momentum: i32,
    pub cutoff_radius_index: usize, // Índice no mesh onde beta vai a zero
    pub data: Vec<f64>,             // O projetor em si
}

impl Pseudopotential {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, UpfError> {
        let content = fs::read_to_string(path)?;
        Self::from_str(&content)
    }

    pub fn from_str(xml_content: &str) -> Result<Self, UpfError> {
        // Parsear o XML
        let doc = Document::parse(xml_content)?;
        let root = doc.root_element();

        // 1. HEADER
        let header_node = root.children()
            .find(|n| n.has_tag_name("PP_HEADER"))
            .ok_or(UpfError::MissingField("PP_HEADER".into()))?;

        let header = Header {
            element: header_node.attribute("element").unwrap_or("X").to_string(),
            z_valence: header_node.attribute("z_valence").unwrap_or("0.0").parse().unwrap_or(0.0),
            mesh_size: header_node.attribute("mesh_size").unwrap_or("0").parse().unwrap_or(0),
            functional: header_node.attribute("functional").unwrap_or("unknown").to_string(),
            number_of_proj: header_node.attribute("number_of_proj").unwrap_or("0").parse().unwrap_or(0),
        };

        // 2. MESH (Grid Radial)
        let mesh_node = root.children()
            .find(|n| n.has_tag_name("PP_MESH"))
            .ok_or(UpfError::MissingField("PP_MESH".into()))?;

        let r_node = mesh_node.children().find(|n| n.has_tag_name("PP_R"))
            .ok_or(UpfError::MissingField("PP_R".into()))?;
        let rab_node = mesh_node.children().find(|n| n.has_tag_name("PP_RAB"))
            .ok_or(UpfError::MissingField("PP_RAB".into()))?;

        let mesh = RadialMesh {
            r: parse_numbers(r_node.text().unwrap_or(""))?,
            rab: parse_numbers(rab_node.text().unwrap_or(""))?,
        };

        // 3. POTENCIAL LOCAL
        let local_node = root.children()
            .find(|n| n.has_tag_name("PP_LOCAL"))
            .ok_or(UpfError::MissingField("PP_LOCAL".into()))?;
        let local = parse_numbers(local_node.text().unwrap_or(""))?;

        // 4. PROJETORES NÃO-LOCAIS (Betas)
        let mut nonlocal = Vec::new();
        if let Some(nl_node) = root.children().find(|n| n.has_tag_name("PP_NONLOCAL")) {
            for child in nl_node.children() {
                if child.tag_name().name().starts_with("PP_BETA") {
                    let l = child.attribute("angular_momentum")
                        .unwrap_or("0").parse().unwrap_or(0);
                    
                    let cutoff = child.attribute("cutoff_radius_index")
                        .unwrap_or("0").parse().unwrap_or(0); // Útil para otimização

                    let data = parse_numbers(child.text().unwrap_or(""))?;

                    nonlocal.push(BetaFunction {
                        index: nonlocal.len(),
                        angular_momentum: l,
                        cutoff_radius_index: cutoff,
                        data,
                    });
                }
            }
        }

        // 5. RHO ATOM (Densidade de Carga Inicial)
        // Essencial para o primeiro passo do SCF
        let rho_atom = if let Some(rho_node) = root.children().find(|n| n.has_tag_name("PP_RHOATOM")) {
            parse_numbers(rho_node.text().unwrap_or(""))?
        } else {
            // Se não tiver, retorna zeros (arriscado, mas evita crash)
            vec![0.0; header.mesh_size]
        };

        // 6. DIJ (Coeficientes de Energia Não-Local)
        // Alguns arquivos colocam isso explicito. 
        // Se não existir, assumimos vazio (trataremos como identidade ou zeros depois)
        let dij = if let Some(dij_node) = root.children().find(|n| n.has_tag_name("PP_DIJ")) {
             parse_numbers(dij_node.text().unwrap_or(""))?
        } else {
             Vec::new()
        };

        Ok(Pseudopotential {
            header,
            mesh,
            local,
            nonlocal,
            rho_atom,
            dij,
        })
    }
}

/// Helper: Converte string gigante de números separada por espaços/novas linhas em Vec<f64>
fn parse_numbers(text: &str) -> Result<Vec<f64>, UpfError> {
    text.split_whitespace()
        .map(|s| s.parse::<f64>().map_err(|_| UpfError::ParseNumber))
        .collect()
}