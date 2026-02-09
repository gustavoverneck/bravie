
pub const PI: f64 = std::f64::consts::PI;

// --- Distância ---
// Unidade interna do Bravie: BOHR (Raio de Bohr)
// 1 Bohr approx 0.529 Å
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;
pub const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_TO_ANGSTROM;

// SI (Metros) e CGS (Centímetros)
pub const BOHR_TO_METER: f64 = BOHR_TO_ANGSTROM * 1.0e-10;
pub const METER_TO_BOHR: f64 = 1.0 / BOHR_TO_METER;

pub const BOHR_TO_CM: f64 = BOHR_TO_METER * 100.0;
pub const CM_TO_BOHR: f64 = 1.0 / BOHR_TO_CM;

pub const BOHR_TO_NM: f64 = BOHR_TO_ANGSTROM * 0.1;
pub const NM_TO_BOHR: f64 = 1.0 / BOHR_TO_NM;

// --- Energia ---
// Unidade interna do Bravie: RYDBERG (Ry)
// 1 Hartree = 27.211386245988 eV
// 1 Hartree = 2.0 Rydberg
pub const HA_TO_EV: f64 = 27.211386245988;
pub const EV_TO_HA: f64 = 1.0 / HA_TO_EV;
pub const HA_TO_RY: f64 = 2.0;
pub const RY_TO_HA: f64 = 0.5;

// SI (Joules) e CGS (Ergs)
pub const HA_TO_JOULE: f64 = 4.3597447222071e-18;
pub const JOULE_TO_HA: f64 = 1.0 / HA_TO_JOULE;

pub const HA_TO_ERG: f64 = HA_TO_JOULE * 1.0e7;
pub const ERG_TO_HA: f64 = 1.0 / HA_TO_ERG;

// Espectroscopia (cm^-1 e THz)
pub const HA_TO_WAVENUMBER: f64 = 219474.63136320; // cm^-1
pub const HA_TO_THZ: f64 = 6579.683920502; // Terahertz

// MASSA (Base: Massa do Elétron / me)
pub const ELECTRON_MASS_SI: f64 = 9.10938356e-31; // kg
pub const PROTON_MASS_AU: f64 = 1836.15267343;    // me
pub const AMU_TO_AU: f64 = 1822.888486209;        // Unidade de Massa Atômica -> me

// Conversão direta de kg para unidades atômicas de massa (me)
pub const KG_TO_AU_MASS: f64 = 1.0 / ELECTRON_MASS_SI;
pub const AU_MASS_TO_KG: f64 = ELECTRON_MASS_SI;

// CARGA ELÉTRICA (Base: Carga Elementar / e)
pub const ELEMENTARY_CHARGE_SI: f64 = 1.602176634e-19; // Coulombs

// Conversão de Coulomb para unidade atômica de carga (e)
pub const COULOMB_TO_AU_CHARGE: f64 = 1.0 / ELEMENTARY_CHARGE_SI;
pub const AU_CHARGE_TO_COULOMB: f64 = ELEMENTARY_CHARGE_SI;

// CGS (Statcoulomb / esu)
// 1 e approx 4.803e-10 esu
pub const AU_CHARGE_TO_ESU: f64 = 4.803204712570263e-10;

// PRESSÃO (Base: Hartree / Bohr^3)
// 1 Atomic Unit of Pressure approx 29.4 Mbar
pub const AU_PRESSURE_TO_PASCAL: f64 = HA_TO_JOULE / (BOHR_TO_METER * BOHR_TO_METER * BOHR_TO_METER);
pub const PASCAL_TO_AU_PRESSURE: f64 = 1.0 / AU_PRESSURE_TO_PASCAL;

pub const AU_PRESSURE_TO_GPA: f64 = AU_PRESSURE_TO_PASCAL * 1.0e-9;
pub const GPA_TO_AU_PRESSURE: f64 = 1.0 / AU_PRESSURE_TO_GPA;

pub const AU_PRESSURE_TO_BAR: f64 = AU_PRESSURE_TO_PASCAL * 1.0e-5;

// Constante de estrutura fina
pub const FINE_STRUCTURE_CONST: f64 = 7.2973525693e-3;