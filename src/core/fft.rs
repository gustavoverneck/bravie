use ndarray::{Array1, Array3};
use ndrustfft::{FftHandler, ndfft_par, ndifft_par};
use num_complex::Complex64;
use rayon::prelude::*; // Importante para o gather paralelo
use crate::core::basis::PlaneWaveBasis;

pub struct FftGrid {
    pub size: [usize; 3],
    
    // Buffers persistentes
    pub buffer: Array3<Complex64>, 
    scratch: Array3<Complex64>, 

    handler_x: FftHandler<f64>,
    handler_y: FftHandler<f64>,
    handler_z: FftHandler<f64>,

    // Mapeamento Linear (Flat Index)
    // Em vez de (u, v, w), guardamos o índice direto na memória linear do buffer.
    // Isso evita calcular (u * ny * nz + v * nz + w) milhões de vezes.
    map_g_to_flat_index: Vec<usize>,
}

impl FftGrid {
    pub fn new(basis: &PlaneWaveBasis) -> Self {
        let (nx, ny, nz) = (basis.fft_grid[0], basis.fft_grid[1], basis.fft_grid[2]);
                
        println!("    FFT Grid init: {}x{}x{}", nx, ny, nz);

        let buffer = Array3::zeros((nx, ny, nz));
        let scratch = Array3::zeros((nx, ny, nz));

        let handler_x = FftHandler::new(nx);
        let handler_y = FftHandler::new(ny);
        let handler_z = FftHandler::new(nz);

        // Pré-cálculo dos strides para indexação linear
        // O layout padrão do ndarray (C-order) é: idx = x*stride_x + y*stride_y + z*stride_z
        // Para Array3::zeros, strides são (ny*nz, nz, 1)
        let stride_x = ny * nz;
        let stride_y = nz;
        let stride_z = 1;

        let inx = nx as i32;
        let iny = ny as i32;
        let inz = nz as i32;

        let mut map_g_to_flat_index = Vec::with_capacity(basis.g_vectors.len());

        for &(ig, jg, kg) in &basis.g_vectors {
            // Wrap around (Periodic Boundary Conditions)
            let u = ((ig % inx) + inx) as usize % nx;
            let v = ((jg % iny) + iny) as usize % ny;
            let w = ((kg % inz) + inz) as usize % nz;
            
            // Cálculo do índice linear (flat) uma única vez na vida
            let flat_idx = u * stride_x + v * stride_y + w * stride_z;
            map_g_to_flat_index.push(flat_idx);
        }

        Self {
            size: [nx, ny, nz],
            buffer,
            scratch,
            handler_x, handler_y, handler_z,
            map_g_to_flat_index,
        }
    }

    /// IFFT: Coeficientes -> Grid -> FFT Inversa -> Buffer Real
    pub fn to_real_space(&mut self, coeffs_recip: &Array1<Complex64>) {
        // Passo 1: Limpar buffer
        self.buffer.fill(Complex64::new(0.0, 0.0));
        
        // CORREÇÃO AQUI:
        // Convertemos ambos para "slices" brutos do Rust (&[T]).
        // Slices têm o método 'get_unchecked' e são mais leves que o ArrayView do ndarray.
        let raw_buffer = self.buffer.as_slice_mut().expect("Buffer deve ser contíguo na memória");
        let raw_coeffs = coeffs_recip.as_slice().expect("Coeffs deve ser contíguo");

        let n_coeffs = coeffs_recip.len();
        
        // Passo 2: Scatter (Loop Unsafe Otimizado)
        for (g_idx, &flat_pos) in self.map_g_to_flat_index.iter().enumerate() {
            if g_idx < n_coeffs {
                unsafe {
                    // Agora estamos chamando get_unchecked em primitivos slices do Rust
                    *raw_buffer.get_unchecked_mut(flat_pos) = *raw_coeffs.get_unchecked(g_idx);
                }
            }
        }
        
        // Passo 3: FFT 3D (Ping-Pong buffers)
        ndifft_par(&self.buffer, &mut self.scratch, &self.handler_x, 0);
        ndifft_par(&self.scratch, &mut self.buffer, &self.handler_y, 1);
        ndifft_par(&self.buffer, &mut self.scratch, &self.handler_z, 2);
        
        // Resultado em scratch -> buffer
        self.buffer.assign(&self.scratch);
    }

    /// FFT: Grid Real -> FFT Forward -> Extrair Coeficientes
    pub fn to_recip_space(&mut self, coeffs_out: &mut Array1<Complex64>) {
        // Passo 1: FFT 3D
        ndfft_par(&self.buffer, &mut self.scratch, &self.handler_x, 0);
        ndfft_par(&self.scratch, &mut self.buffer, &self.handler_y, 1);
        ndfft_par(&self.buffer, &mut self.scratch, &self.handler_z, 2);

        // Agora o resultado está em 'scratch'.
        let raw_scratch = self.scratch.as_slice().expect("Buffer deve ser contíguo");

        // OTIMIZAÇÃO 3: Gather Paralelo
        // Diferente da escrita, a leitura pode ser feita em paralelo trivialmente!
        // Usamos Rayon para preencher 'coeffs_out' em paralelo.
        
        // Se coeffs_out e map tiverem o mesmo tamanho (deveriam):
        coeffs_out.as_slice_mut().expect("Coeffs contíguo")
            .par_iter_mut()
            .zip(&self.map_g_to_flat_index) // Zipa com o índice de onde ler
            .for_each(|(out_val, &flat_idx)| {
                // Leitura unsafe também é válida e rápida, mas aqui o ganho maior é o paralelismo
                unsafe {
                    *out_val = *raw_scratch.get_unchecked(flat_idx);
                }
            });
    }
}