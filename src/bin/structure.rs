use bravie::core::structure::{Structure, Species};
use bravie::utils::welcome::{print_welcome};

fn main() {
    print_welcome();
    let a = 4.65; // Parâmetro de rede do grafeno em Bohr (aprox. 2.46 Angstrom)
    let c = 15.0; // Vácuo na direção Z

    let carbon = Species {
        id: 1,
        element: "C".to_string(),
        atomic_number: 6,
        mass: 12.011,
        pseudo_path: "C.upf".to_string(),
    };

    // Vetores da rede hexagonal
    let a1 = [a, 0.0, 0.0];
    let a2 = [-a / 2.0, a * (3.0f64).sqrt() / 2.0, 0.0];
    let a3 = [0.0, 0.0, c];

    let structure = Structure::builder()
        .lattice(a1, a2, a3)
        .add_species(carbon)
        .add_atom([0.0, 0.0, 0.0], 1) // Átomo de Carbono 1
        .add_atom([a / 2.0, a * (3.0f64).sqrt() / 6.0, 0.0], 1) // Átomo de Carbono 2 (1/3, 1/3 em coordenadas internas)
        .build()
        .expect("Falha ao construir a estrutura de grafeno");

    println!("Estrutura de Grafeno criada com sucesso:");
    println!("{}", structure);
}
