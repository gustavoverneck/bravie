use crate::io::upf::RadialMesh;

/// Interpola linearmente uma função radial f(r) definida em um mesh.
pub fn interpolate_linear(r: f64, mesh: &RadialMesh, data: &[f64]) -> f64 {
    // Se r estiver fora do alcance, retorna 0.0 (assumindo decaimento)
    // Para potenciais de longo alcance (Coulomb), isso assume V_loc -> 0 no cutoff.
    if r > mesh.r.last().copied().unwrap_or(0.0) {
        return 0.0;
    }

    // Tratamento para r -> 0
    if r < 1e-6 {
        if !data.is_empty() { return data[0]; }
        return 0.0;
    }

    // Busca Binária
    let idx = match mesh.r.binary_search_by(|val| val.partial_cmp(&r).unwrap()) {
        Ok(i) => i,
        Err(i) => if i > 0 { i - 1 } else { 0 },
    };

    if idx >= mesh.r.len() - 1 {
        return data.last().copied().unwrap_or(0.0);
    }

    // Interpolação Linear
    let r1 = mesh.r[idx];
    let r2 = mesh.r[idx+1];
    let y1 = data[idx];
    let y2 = data[idx+1];

    let t = (r - r1) / (r2 - r1);
    y1 + t * (y2 - y1)
}