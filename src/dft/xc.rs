use std::f64::consts::PI;
use ndarray::Array3;

/// Calcula o potencial e a energia de Troca-Correlação (XC) usando LDA (Perdew-Zunger 81).
/// Retorna uma tupla: (V_xc_potential, E_xc_total_energy)
pub fn calculate_xc_lda(rho: &Array3<f64>, volume: f64) -> (Array3<f64>, f64) {
    let (nx, ny, nz) = rho.dim();
    let mut v_xc = Array3::<f64>::zeros((nx, ny, nz));
    let mut total_energy = 0.0;
    
    // Volume por voxel para integração
    let dvol = volume / (nx * ny * nz) as f64;

    for ((i, j, k), &n) in rho.indexed_iter() {
        if n < 1e-12 {
            // Densidade nula -> Potencial nulo
            v_xc[[i, j, k]] = 0.0;
            continue;
        }

        // 1. Raio de Wigner-Seitz (rs)
        // n = 3 / (4 * pi * rs^3)  =>  rs = (3 / (4 * pi * n))^(1/3)
        let rs = (3.0 / (4.0 * PI * n)).powf(1.0 / 3.0);

        // 2. Exchange (Slater)
        let (ex, vx) = lda_exchange(n);

        // 3. Correlation (Perdew-Zunger 81)
        let (ec, vc) = lda_correlation_pz81(rs);

        // Soma total
        let e_total_local = ex + ec; // Energia por partícula
        let v_total_local = vx + vc; // Potencial (derivada funcional)

        v_xc[[i, j, k]] = v_total_local;
        
        // Integral E_xc = sum( n(r) * eps_xc(r) ) * dVol
        total_energy += n * e_total_local;
    }

    (v_xc, total_energy * dvol)
}

/// LDA Exchange (Slater Dirac)
/// Retorna (epsilon_x, v_x) em Ry
fn lda_exchange(n: f64) -> (f64, f64) {
    // E_x = -3/4 * (3/pi)^(1/3) * n^(1/3)
    // V_x = - (3*n/pi)^(1/3)
    // Em Rydberg: multiplicar fator se necessário? 
    // As fórmulas padrão geralmente são em Hartree.
    // Hartree -> Ry: Multiplicar energia por 2.
    
    // Formula em Hartree (a.u.):
    // Ex = -3/4 * (3/pi)^(1/3) * n^(1/3)
    // Vx = 4/3 * Ex = -(3/pi)^(1/3) * n^(1/3)
    
    // Vamos calcular em Ry direto (fator 2 na energia):
    let third = 1.0 / 3.0;
    let factor = (3.0 / PI).powf(third);
    
    let v_x_hartree = -factor * n.powf(third);
    let e_x_hartree = 0.75 * v_x_hartree; // = 3/4 * Vx

    (2.0 * e_x_hartree, 2.0 * v_x_hartree)
}

/// LDA Correlation (Perdew-Zunger 1981)
/// Retorna (epsilon_c, v_c) em Ry
fn lda_correlation_pz81(rs: f64) -> (f64, f64) {
    // Constantes do paper PZ81 (para Paramagnetico/Spin-unpolarized)
    let a = 0.0311;
    let b = -0.048;
    let c = 0.0020;
    let d = -0.0116;
    
    let gamma = -0.1423;
    let beta1 = 1.0529;
    let beta2 = 0.3334;

    let e_c_hartree: f64;
    let v_c_hartree: f64;

    if rs < 1.0 {
        // Regime de alta densidade (ln)
        let ln_rs = rs.ln();
        e_c_hartree = a * ln_rs + b + c * rs * ln_rs + d * rs;
        
        // Derivada d(n*eps)/dn = eps - (rs/3) * d(eps)/d(rs)
        let de_drs = a / rs + c * ln_rs + c + d;
        v_c_hartree = e_c_hartree - (rs / 3.0) * de_drs;
    } else {
        // Regime de baixa densidade (raiz)
        let sqrt_rs = rs.sqrt();
        let denom = 1.0 + beta1 * sqrt_rs + beta2 * rs;
        e_c_hartree = gamma / denom;
        
        let de_drs = -gamma * (0.5 * beta1 / sqrt_rs + beta2) / (denom * denom);
        v_c_hartree = e_c_hartree - (rs / 3.0) * de_drs;
    }

    (2.0 * e_c_hartree, 2.0 * v_c_hartree)
}