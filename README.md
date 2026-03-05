# Bravie
Um solver de Teoria do Funcional da Densidade (DFT) escrito em Rust, focado em alto desempenho e na clareza de implementação.

**Referência Geral para a Arquitetura do Código:**
* Martin, R. M. (2004). *Electronic Structure: Basic Theory and Practical Methods*. Cambridge University Press.
* Payne, M. C., et al. (1992). *Iterative minimization techniques for ab initio total-energy calculations: molecular dynamics and conjugate gradients*. Reviews of Modern Physics, 64(4), 1045. 
* Ashcroft, Mermin. (2022). *Física do Estado Sólido*. Cengage Learning.
* Oliveira, de Jesus. (2017). *Introdução à Física do Estado Sólido*. LF Editorial.
* Giustino (2014). Materials Modelling using Density Functional Theory. Oxford.
---

## Roadmap & Progresso

### Fase 1: Fundação e Estrutura
Estabelecimento das estruturas de dados base e leitura de inputs físicos.
- [x] **Sistema de Build e Módulos:** Estrutura do projeto orientada a alta performance em Rust (`core`, `io`, `utils`).
- [x] **Estrutura Cristalina:** Definição vetorial rigorosa da rede recíproca, `Atom`, `Species`, `Lattice` e `Structure`.
- [x] **Leitura de Pseudopotenciais:** Parser para formato UPF (Unified Pseudopotential Format) v2 extraindo malhas radiais, projetores não-locais e funções de onda atômicas.
- [x] **Configuração da Simulação:** Builder Pattern (`SimulationBuilder`) para gerenciar as dependências do estado do sistema.
- [x] **Constantes Físicas:** Sistema interno em Unidades Atômicas (Rydberg/Bohr) com conversões robustas.

### Fase 2: Motor Físico (Basis Set e Espaço Recíproco)
Construção da malha de integração e do motor de transformadas de Fourier, que são o coração de um código plane-wave.
- [x] **Geração de K-Points na Zona de Brillouin:**
    - Malha uniforme para SCF (Monkhorst-Pack) e gerador de caminhos de alta simetria para bandas.
    - *Ref: Monkhorst, H. J., & Pack, J. D. (1976). Special points for Brillouin-zone integrations. Physical Review B, 13(12), 5188.*
- [x] **Base de Ondas Planas (Plane Waves):**
    - Geração de vetores $\mathbf{G}$ limitados pela energia de corte ($E_{cut}$) e mapeamento direto entre os índices espaciais (G-Space) e a malha real.
- [x] **FFT 3D Otimizada:**
    - Uso eficiente de `ndrustfft` com mapeamento linear e gather/scatter paralelo via `rayon` para transição ultrarrápida entre $\psi(\mathbf{r})$ e $\psi(\mathbf{G})$.

### Fase 3: Construção do Hamiltoniano de Kohn-Sham
Implementação dos operadores que atuam sobre as funções de onda.
- [ ] **Energia Cinética e Equação de Poisson (Hartree):** - Atuação diagonal no espaço recíproco: $\hat{T} = \frac{1}{2}|\mathbf{k}+\mathbf{G}|^2$ e $V_H(\mathbf{G}) = 4\pi\rho(\mathbf{G})/|\mathbf{G}|^2$.
- [ ] **Potencial Local ($V_{loc}$):** - Transformada de Fourier esférica e interpolação spline do potencial radial de valência fornecido pelo UPF.
- [ ] **Potencial Não-Local ($V_{nl}$):** - Implementação da forma separável de Kleinman-Bylander, calculando os fatores de estrutura e o produto interno na base de ondas planas.
    - *Ref: Kleinman, L., & Bylander, D. M. (1982). Efficacious Form for Model Pseudopotentials. Physical Review Letters, 48(20), 1425.*
- [ ] **Troca e Correlação (Exchange-Correlation - $V_{xc}$):**
    - Implementação do funcional LDA paramétrico para começar a fechar o ciclo SCF.
    - *Ref: Perdew, J. P., & Zunger, A. (1981). Self-interaction correction to density-functional approximations for many-electron systems. Physical Review B, 23(10), 5048.*

### Fase 4: O Solver (Ciclo Auto-Consistente - SCF)
Resolução iterativa do problema de minimização de autovalores.
- [x] **Densidade Inicial (SAD):** - Chute inicial robusto baseado na Superposição de Densidades Atômicas reais (SAD), garantindo neutralidade e acelerando convergência.
- [ ] **Diagonalização Iterativa (Eigensolver):** - Implementação do método LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient) ou Davidson com precondicionamento de Payne/Teter focado na energia cinética.
    - *Ref: Knyazev, A. V. (2001). Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method. SIAM Journal on Scientific Computing, 23(2), 517-541.*
- [ ] **Mixing de Densidade:** - Implementação do esquema de Broyden (ou Pulay DIIS) para atualizar iterativamente $\rho(\mathbf{r})$ prevenindo divergências de sloshing de carga (comum em metais).
    - *Ref: Johnson, D. D. (1988). Modified Broyden’s method for calculating charge densities. Physical Review B, 38(18), 12807.*
- [ ] **Cálculo da Energia Total:**
    - Soma da energia de banda (autovalores), menos o duplo cômputo de Hartree e XC, somada à energia de Ewald (interação íon-íon de longo alcance).
    - *Ref: Ewald, P. P. (1921). Die Berechnung optischer und elektrostatischer Gitterpotentiale. Annalen der Physik, 369(3), 253-287.*

### Fase 5: Propriedades Físicas e Relaxação Estrutural
Ferramentas de pós-processamento e otimização geométrica.
- [ ] **Cálculos NSCF e Estrutura de Bandas:** - Congelamento da densidade (NSCF loop) para resolver autovalores em caminhos de alta simetria. Introdução de Smearing (Fermi-Dirac/Methfessel-Paxton) para sistemas metálicos.
    - *Ref: Methfessel, M., & Paxton, A. T. (1989). High-precision sampling for Brillouin-zone integration in metals. Physical Review B, 40(6), 3616.*
- [ ] **Teorema de Hellmann-Feynman (Cálculo de Forças):** - Derivação analítica das forças Ewald, Locais e Não-Locais agindo sobre os íons com base na densidade de estado fundamental.
    - *Ref: Feynman, R. P. (1939). Forces in Molecules. Physical Review, 56(4), 340.*
- [ ] **Otimização Geométrica (Relax):** - Implementação do algoritmo BFGS para atualizar as coordenadas atômicas visando o mínimo global de energia do sistema.
    - *Ref: Pfrommer, B. G., et al. (1997). Relaxation of crystals with the quasi-Newton method. Journal of Computational Physics, 131(1), 233-240.*
- [ ] **Exportação de Densidades e Orbitais:**
    - Suporte nativo para exportar grids no formato `.cube` compatível com softwares como VESTA e XCrySDen.
