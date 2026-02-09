# Bravie
Um solver de DFT em Rust, em desenvolvimento.

## Roadmap & Progresso

### Fase 1: Fundação e Estrutura
- [x] **Sistema de Build:** Configuração do Cargo e estrutura de módulos (`core`, `io`, `utils`).
- [x] **Estrutura Cristalina:** Definição de `Atom`, `Species`, `Lattice` e `Structure`.
- [x] **Leitura de Pseudopotenciais:** Parser para formato UPF (Unified Pseudopotential Format) v2.
- [x] **Configuração da Simulação:** Builder Pattern para configurar o sistema (`SimulationBuilder`).
- [x] **Constantes Físicas:** Módulo de conversão de unidades (Bohr, Hartree, Ry, Angstrom).
- [ ] **Tratamento de Erros:** Implementação robusta com `thiserror` (Em andamento).

### Fase 2: Motor Físico (Basis Set)
- [ ] **Geração de K-Points:**
    - [ ] Malha Monkhorst-Pack (SCF).
    - [ ] Caminhos de Banda (Band Paths) para visualização.
- [ ] **Base de Ondas Planas:**
    - [ ] Geração de vetores G dentro da esfera de corte ($E_{cut}$).
    - [ ] Mapeamento G-Space <-> FFT Grid.
- [ ] **FFT 3D:** Integração com biblioteca FFT (`ndrustfft`).

### Fase 3: Hamiltoniano e Potenciais
- [ ] **Energia Cinética:** Cálculo do termo $\frac{1}{2}|k+G|^2$.
- [ ] **Potencial Local ($V_{loc}$):** Interpolação e Transformada de Fourier do potencial radial.
- [ ] **Potencial Não-Local ($V_{nl}$):** Projetores beta e coeficientes $D_{ij}$.
- [ ] **Potencial de Hartree ($V_{H}$):** Solução da equação de Poisson no espaço recíproco.
- [ ] **Exchange-Correlation ($V_{xc}$):** Implementação básica (LDA - Perdew-Zunger ou PZ81).

### Fase 4: O Solver (SCF Loop)
- [ ] **Diagonalização Iterativa:** Algoritmo de Davidson ou Conjugate Gradient (CG).
- [ ] **Cálculo de Densidade:** $\rho(r) = \sum |\psi(r)|^2$.
- [ ] **Mixing de Densidade:** Algoritmo de Broyden ou Anderson para acelerar convergência.
- [ ] **Ciclo Auto-Consistente:** Loop `Potencial -> Diagonalização -> Densidade -> Mixing`.
- [ ] **Cálculo de Energia Total:** Soma dos componentes (Cinética, Hartree, XC, Ewald, etc.).

### Fase 5: Funcionalidades Avançadas
- [ ] **Estrutura de Bandas:** Cálculo não-autoconsistente (NSCF) em caminho de alta simetria.
- [ ] **Exportação de Dados:** Arquivos `.cube` para visualizar densidade e orbitais (ex: VESTA).
- [ ] **Cálculo de Forças:** Teorema de Hellmann-Feynman.
- [ ] **Relaxação Geométrica:** Otimização da posição dos átomos (BFGS).