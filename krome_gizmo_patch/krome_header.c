#include "../GIZMO_config.h"
const int krome_idx_E = 0; // E
const int krome_idx_Hk = 1; // H-
const int krome_idx_QE = 2; // QE
const int krome_idx_QHk = 3; // QH-
const int krome_idx_H = 4; // H
const int krome_idx_HE = 5; // HE
const int krome_idx_H2 = 6; // H2
const int krome_idx_QH = 7; // QH
const int krome_idx_QH2 = 8; // QH2
const int krome_idx_Hj = 9; // H+
const int krome_idx_HEj = 10; // HE+
const int krome_idx_H2j = 11; // H2+
const int krome_idx_QHj = 12; // QH+
const int krome_idx_QH2j = 13; // QH2+
const int krome_idx_QH3j = 14; // QH3+
const int krome_idx_HEjj = 15; // HE++
const int krome_idx_CR = 16; // CR
const int krome_idx_g = 17; // g
const int krome_idx_qg = 18; // qg
const int krome_idx_Tgas = 19; // Tgas
const int krome_idx_dummy = 20; // dummy
const char* krome_names[] = {
  "E",
  "H-",
  "QE",
  "QH-",
  "H",
  "HE",
  "H2",
  "QH",
  "QH2",
  "H+",
  "HE+",
  "H2+",
  "QH+",
  "QH2+",
  "QH3+",
  "HE++",
  "CR",
  "g",
  "qg",
  "Tgas",
  "dummy"
};

const int krome_idx_cool_h2 = 0;
const int krome_idx_cool_h2gp = 1;
const int krome_idx_cool_atomic = 2;
const int krome_idx_cool_cen = 2;
const int krome_idx_cool_hd = 3;
const int krome_idx_cool_z = 4;
const int krome_idx_cool_metal = 4;
const int krome_idx_cool_dh = 5;
const int krome_idx_cool_enthalpic = 5;
const int krome_idx_cool_dust = 6;
const int krome_idx_cool_compton = 7;
const int krome_idx_cool_cie = 8;
const int krome_idx_cool_continuum = 9;
const int krome_idx_cool_cont = 9;
const int krome_idx_cool_exp = 10;
const int krome_idx_cool_expansion = 10;
const int krome_idx_cool_ff = 11;
const int krome_idx_cool_bss = 11;
const int krome_idx_cool_custom = 12;
const int krome_idx_cool_co = 13;
const int krome_idx_cool_zcie = 14;
const int krome_idx_cool_zcienouv = 15;
const int krome_idx_cool_zextend = 16;
const int krome_idx_cool_gh = 17;
const int krome_idx_cool_oh = 18;
const int krome_idx_cool_h2o = 19;
const int krome_idx_cool_hcn = 20;
const int krome_idx_cool_darkmol = 21;
const int krome_idx_cool_darkatom = 22;
const int krome_idx_cool_adarkatom = 23;
const int krome_ncools = 24;

const int krome_idx_heat_chem = 0;
const int krome_idx_heat_compress = 1;
const int krome_idx_heat_compr = 1;
const int krome_idx_heat_photo = 2;
const int krome_idx_heat_dh = 3;
const int krome_idx_heat_enthalpic = 3;
const int krome_idx_heat_photoav = 4;
const int krome_idx_heat_av = 4;
const int krome_idx_heat_cr = 5;
const int krome_idx_heat_dust = 6;
const int krome_idx_heat_xray = 7;
const int krome_idx_heat_visc = 8;
const int krome_idx_heat_viscous = 8;
const int krome_idx_heat_custom = 9;
const int krome_idx_heat_zcie = 10;
const int krome_idx_heat_adarkatom = 11;
const int krome_idx_heat_darkmol = 12;
const int krome_nheats = 13;

const int krome_nrea=45;
const int krome_nmols=16;
const int krome_nspec=21;
const int krome_natoms=5;
const int krome_ndust=0;
const int krome_ndustTypes=0;
const int krome_nPhotoBins=0;
const int krome_nPhotoRates=0;

const double krome_boltzmann_eV = 8.617332478e-5; //eV / K
const double krome_boltzmann_J = 1.380648e-23; //J / K
const double krome_boltzmann_erg = 1.380648e-16; //erg / K
const double krome_iboltzmann_eV = 1e0/8.617332478e-5; //K / eV
const double krome_iboltzmann_erg = 1e0/1.380648e-16; //K / erg
const double krome_planck_eV = 4.135667516e-15; //eV s
const double krome_planck_J = 6.62606957e-34; //J s
const double krome_planck_erg = 6.62606957e-27; //erg s
const double krome_iplanck_eV = 1e0/4.135667516e-15; //1 / eV / s
const double krome_iplanck_J = 1e0/6.62606957e-34; //1 / J / s
const double krome_iplanck_erg = 1e0/6.62606957e-27; //1 / erg / s
const double krome_gravity = 6.674e-8; //cm3 / g / s2
const double krome_e_mass = 9.10938188e-28; //g
const double krome_p_mass = 1.67262158e-24; //g
const double krome_n_mass = 1.674920e-24; //g
const double krome_ip_mass = 1e0/1.67262158e-24; //1/g
const double krome_clight = 2.99792458e10; //cm/s
const double krome_pi = 3.14159265359e0; //#
const double krome_eV_to_erg = 1.60217646e-12; //eV -> erg
const double krome_ry_to_eV = 13.60569e0; //rydberg -> eV
const double krome_ry_to_erg = 2.179872e-11; //rydberg -> erg
const double krome_seconds_per_year = 365e0*24e0*3600e0; //yr -> s
const double krome_km_to_cm = 1e5; //km -> cm
const double krome_cm_to_Mpc = 1.e0/3.08e24; //cm -> Mpc
const double krome_kvgas_erg = 8.e0*1.380648e-16/3.14159265359e0/1.67262158e-24; //
const double krome_pre_kvgas_sqrt = 1.87504433599e-08; //
const double krome_pre_planck = 1.4744993357e-47; //erg/cm2*s3
const double krome_exp_planck = 6.62606957e-27 / 1.380648e-16; //s*K
const double krome_stefboltz_erg = 5.670373e-5; //erg/s/cm2/K4
const double krome_N_avogadro = 6.0221e23; //#
const double krome_Rgas_J = 8.3144621e0; //J/K/mol
const double krome_Rgas_kJ = 8.3144621e-3; //kJ/K/mol
const double krome_hubble = 0.704e0; //dimensionless
const double krome_Omega0 = 1.0e0; //dimensionless
const double krome_Omegab = 0.0456e0; //dimensionless
const double krome_Hubble0 = 1.e2*0.704e0*1e5*1.e0/3.08e24; //1/s
const double krome_keV_to_g = 1.782661907e-30; //keV/c^2 -> g
const double krome_g_to_keV = 1.e0/1.782661907e-30; //g -> keV/c^2
const double krome_fine_structure_constant = 7.2973525664e-3; //#
#ifdef KROME_RE
const double krome_qe_mass = KROME_RE * krome_e_mass; // !g
#else
const double krome_qe_mass = 1.000000e+01*9.10938188e-28 ; //9.10938188d-28 !g
#endif
#ifdef KROME_RP
const double krome_qp_mass = KROME_RP * krome_p_mass; // !g
const double krome_iqp_mass = 1.e0/(krome_qp_mass); //1/g
#else
const double krome_qp_mass = 1.000000e+02*1.67262158e-24 ; //1.67262158d-24 !g
const double krome_iqp_mass = 1.e0/(1.000000e+02*1.67262158e-24 ); //1/g
#endif
#ifdef KROME_RN
const double krome_qn_mass = KROME_RN * krome_n_mass; //  !g
#else
const double krome_qn_mass = 1.674920e-24; //g
#endif
#ifdef KROME_RA
const double krome_Dalpha = KROME_RA * krome_fine_structure_constant; 
#else
const double krome_Dalpha = 1.000000e+00*7.2973525664e-3 ; //7.2973525664d-3 !dark fine structure constant
#endif
#ifdef XI
const double krome_xi = XI; 
#else
const double krome_xi = 1.e-2; //Dark CMB temperature to Normal CMB temperature ratio
#endif
