// src/bin/help.rs
use bravie::utils::welcome::print_welcome;

// FUTURE USAGE

// fn print_usage() {
//     println!("USO:");
//     println!("    bravie [COMANDO] [OPÇÕES] <ARQUIVO_INPUT>\n");
// }

// fn print_commands() {
//     println!("COMANDOS DISPONÍVEIS:");
//     println!("    run       Executa uma simulação completa (leitura -> scf -> output)");
//     println!("    scf       Executa apenas o ciclo de Autoconsistência (Self-Consistent Field)");
//     println!("    bands     Calcula a estrutura de bandas (requer densidade convergida)");
//     println!("    relax     Otimiza a geometria do sistema (minimiza forças)");
//     println!("    check     Verifica se o arquivo de input e pseudopotenciais são válidos");
//     println!("    help      Mostra esta mensagem de ajuda\n");
// }

// fn print_options() {
//     println!("OPÇÕES GERAIS:");
//     println!("    -i, --input <FILE>    Define o arquivo de entrada (.toml, .json)");
//     println!("    -o, --output <FILE>   Define o arquivo de saída (padrão: stdout)");
//     println!("    -v, --verbose         Ativa logs detalhados (debug mode)");
//     println!("    -p, --parallel <N>    Define número de threads (OpenMP/Rayon)");
//     println!("    --version             Mostra a versão atual do Bravie\n");
// }

// fn print_examples() {
//     println!("EXEMPLOS:");
//     println!("    # Rodar uma simulação simples de grafeno");
//     println!("    bravie run -i grafeno.toml -o grafeno.out");
//     println!("");
//     println!("    # Verificar se os pseudopotenciais existem");
//     println!("    bravie check -i silicio.toml");
//     println!("");
//     println!("    # Rodar relaxação com log detalhado");
//     println!("    bravie relax --input agua.toml --verbose\n");
// }

fn print_footer() {
    println!("------------------------------------------------------------");
    println!("Desenvolvido em Rust por Gustavo Verneck.");
    println!("Repositório: https://github.com/gustavoverneck/bravie");
    println!("------------------------------------------------------------");
}

fn main() {
    print_welcome();
    // print_usage();
    // print_commands();
    // print_options();
    // print_examples();
    // print_footer();
}