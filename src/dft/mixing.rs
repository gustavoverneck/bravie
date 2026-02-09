use ndarray::Array3;
use nalgebra::{DMatrix, DVector};

/// Gerencia o histórico para Anderson Mixing (Pulay).
pub struct AndersonMixer {
    beta: f64,              // Fator de mistura (0.0 a 1.0)
    history_size: usize,    // Tamanho do histórico (M)
    
    // Histórico das densidades de ENTRADA (rho_in)
    rho_in_history: Vec<Array3<f64>>,
    
    // Histórico dos resíduos (R = rho_out - rho_in)
    residual_history: Vec<Array3<f64>>,
}

impl AndersonMixer {
    pub fn new(beta: f64, history_size: usize) -> Self {
        Self {
            beta,
            history_size,
            rho_in_history: Vec::new(),
            residual_history: Vec::new(),
        }
    }

    /// Calcula a próxima densidade de entrada (rho_next) baseada na saída atual (rho_out).
    pub fn mix(&mut self, rho_in: &Array3<f64>, rho_out: &Array3<f64>) -> Array3<f64> {
        let residual = rho_out - rho_in;

        // Gerencia histórico
        if self.rho_in_history.len() >= self.history_size {
            self.rho_in_history.remove(0);
            self.residual_history.remove(0);
        }

        self.rho_in_history.push(rho_in.clone());
        self.residual_history.push(residual.clone());

        let m = self.rho_in_history.len();

        if m <= 1 {
            return rho_in + &(residual.mapv(|x| x * self.beta));
        }

        // Monta Matriz A e Vetor B
        let mut a_mat = DMatrix::<f64>::zeros(m - 1, m - 1);
        let mut b_vec = DVector::<f64>::zeros(m - 1);
        let r_k = &self.residual_history[m - 1];

        for i in 0..(m - 1) {
            let dr_i = &self.residual_history[i] - r_k;
            b_vec[i] = -(&dr_i * r_k).sum();
            
            for j in i..(m - 1) {
                let dr_j = &self.residual_history[j] - r_k;
                let val = (&dr_i * &dr_j).sum();
                a_mat[(i, j)] = val;
                a_mat[(j, i)] = val;
            }
        }

        // --- CORREÇÃO CRÍTICA: REGULARIZAÇÃO ---
        // Adiciona um valor pequeno (1e-8) na diagonal para evitar singularidade
        // quando o sistema está convergindo (resíduos paralelos).
        for i in 0..(m - 1) {
            a_mat[(i, i)] += 1e-8; 
        }

        // Tenta inverter
        let alpha = match a_mat.try_inverse() {
            Some(inv) => inv * b_vec,
            None => {
                // Se falhar mesmo com regularização, limpa histórico e faz passo linear
                self.rho_in_history.clear();
                self.residual_history.clear();
                return rho_in + &(residual.mapv(|x| x * self.beta));
            }
        };

        // Constrói Rho Otimizado
        let mut rho_opt = self.rho_in_history[m - 1].clone();
        let mut res_opt = self.residual_history[m - 1].clone();

        for i in 0..(m - 1) {
            let coeff = alpha[i];
            // Proteção contra coeficientes explosivos
            if coeff.abs() > 10.0 {
                // Se o Anderson sugerir um passo gigante, aborte e use linear
                return rho_in + &(residual.mapv(|x| x * self.beta));
            }
            
            let d_rho = &self.rho_in_history[i] - &self.rho_in_history[m - 1];
            let d_res = &self.residual_history[i] - &self.residual_history[m - 1];
            rho_opt = rho_opt + d_rho.mapv(|x| x * coeff);
            res_opt = res_opt + d_res.mapv(|x| x * coeff);
        }

        rho_opt + res_opt.mapv(|x| x * self.beta)
    }
}