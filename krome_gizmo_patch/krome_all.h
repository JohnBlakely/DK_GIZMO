/* c header for krome.f90/the krome_main module.
  variables for which krome will return an updated value need to be passed
  by reference using an argument pointer, e.g., "*x".
  arrays as input should also be passed by reference/argument pointer.
  passing 2-d arrays as arguments use an array of pointers to an array (i.e., "**x").
  to return an array from a function, it must return an argument pointer,
  e.g., "extern double *functionname()".
*/


extern void krome(double* x, double rhogas, double *tgas, double *dt);
extern void krome_equilibrium(double* x, double rhogas, double tgas, int* verbosity);
extern void krome_init();
extern void krome_get_coe(double* x, double *tgas, double* krome_get_coe_var);
extern void krome_get_coet(double *tgas, double* krome_get_coet_var);

/* c header for krome_user.f90/the krome_user module.
  variables for which krome will return an updated value need to be passed
  by reference using an argument pointer, e.g., "*x".
  arrays as input should also be passed by reference/argument pointer.
  passing 2-d arrays as arguments use an array of pointers to an array (i.e., "**x").
  to return an array from a function, it must return an argument pointer,
  e.g., "extern double *functionname()".
*/

extern const int krome_idx_E; // E
extern const int krome_idx_Hk; // H-
extern const int krome_idx_QE; // QE
extern const int krome_idx_QHk; // QH-
extern const int krome_idx_H; // H
extern const int krome_idx_HE; // HE
extern const int krome_idx_H2; // H2
extern const int krome_idx_QH; // QH
extern const int krome_idx_QH2; // QH2
extern const int krome_idx_Hj; // H+
extern const int krome_idx_HEj; // HE+
extern const int krome_idx_H2j; // H2+
extern const int krome_idx_QHj; // QH+
extern const int krome_idx_QH2j; // QH2+
extern const int krome_idx_QH3j; // QH3+
extern const int krome_idx_HEjj; // HE++
extern const int krome_idx_CR; // CR
extern const int krome_idx_g; // g
extern const int krome_idx_qg; // qg
extern const int krome_idx_Tgas; // Tgas
extern const int krome_idx_dummy; // dummy
extern const char* krome_names[];

extern const int krome_idx_cool_h2;
extern const int krome_idx_cool_h2gp;
extern const int krome_idx_cool_atomic;
extern const int krome_idx_cool_cen;
extern const int krome_idx_cool_hd;
extern const int krome_idx_cool_z;
extern const int krome_idx_cool_metal;
extern const int krome_idx_cool_dh;
extern const int krome_idx_cool_enthalpic;
extern const int krome_idx_cool_dust;
extern const int krome_idx_cool_compton;
extern const int krome_idx_cool_cie;
extern const int krome_idx_cool_continuum;
extern const int krome_idx_cool_cont;
extern const int krome_idx_cool_exp;
extern const int krome_idx_cool_expansion;
extern const int krome_idx_cool_ff;
extern const int krome_idx_cool_bss;
extern const int krome_idx_cool_custom;
extern const int krome_idx_cool_co;
extern const int krome_idx_cool_zcie;
extern const int krome_idx_cool_zcienouv;
extern const int krome_idx_cool_zextend;
extern const int krome_idx_cool_gh;
extern const int krome_idx_cool_oh;
extern const int krome_idx_cool_h2o;
extern const int krome_idx_cool_hcn;
extern const int krome_idx_cool_darkmol;
extern const int krome_idx_cool_darkatom;
extern const int krome_idx_cool_adarkatom;
extern const int krome_ncools;

extern const int krome_idx_heat_chem;
extern const int krome_idx_heat_compress;
extern const int krome_idx_heat_compr;
extern const int krome_idx_heat_photo;
extern const int krome_idx_heat_dh;
extern const int krome_idx_heat_enthalpic;
extern const int krome_idx_heat_photoav;
extern const int krome_idx_heat_av;
extern const int krome_idx_heat_cr;
extern const int krome_idx_heat_dust;
extern const int krome_idx_heat_xray;
extern const int krome_idx_heat_visc;
extern const int krome_idx_heat_viscous;
extern const int krome_idx_heat_custom;
extern const int krome_idx_heat_zcie;
extern const int krome_idx_heat_adarkatom;
extern const int krome_idx_heat_darkmol;
extern const int krome_nheats;

extern const int krome_nrea;
extern const int krome_nmols;
extern const int krome_nspec;
extern const int krome_natoms;
extern const int krome_ndust;
extern const int krome_ndustTypes;
extern const int krome_nPhotoBins;
extern const int krome_nPhotoRates;

extern const double krome_boltzmann_eV; //eV / K
extern const double krome_boltzmann_J; //J / K
extern const double krome_boltzmann_erg; //erg / K
extern const double krome_iboltzmann_eV; //K / eV
extern const double krome_iboltzmann_erg; //K / erg
extern const double krome_planck_eV; //eV s
extern const double krome_planck_J; //J s
extern const double krome_planck_erg; //erg s
extern const double krome_iplanck_eV; //1 / eV / s
extern const double krome_iplanck_J; //1 / J / s
extern const double krome_iplanck_erg; //1 / erg / s
extern const double krome_gravity; //cm3 / g / s2
extern const double krome_e_mass; //g
extern const double krome_p_mass; //g
extern const double krome_n_mass; //g
extern const double krome_ip_mass; //1/g
extern const double krome_clight; //cm/s
extern const double krome_pi; //#
extern const double krome_eV_to_erg; //eV -> erg
extern const double krome_ry_to_eV; //rydberg -> eV
extern const double krome_ry_to_erg; //rydberg -> erg
extern const double krome_seconds_per_year; //yr -> s
extern const double krome_km_to_cm; //km -> cm
extern const double krome_cm_to_Mpc; //cm -> Mpc
extern const double krome_kvgas_erg; //
extern const double krome_pre_kvgas_sqrt; //
extern const double krome_pre_planck; //erg/cm2*s3
extern const double krome_exp_planck; //s*K
extern const double krome_stefboltz_erg; //erg/s/cm2/K4
extern const double krome_N_avogadro; //#
extern const double krome_Rgas_J; //J/K/mol
extern const double krome_Rgas_kJ; //kJ/K/mol
extern const double krome_hubble; //dimensionless
extern const double krome_Omega0; //dimensionless
extern const double krome_Omegab; //dimensionless
extern const double krome_Hubble0; //1/s
extern const double krome_keV_to_g; //keV/c^2 -> g
extern const double krome_g_to_keV; //g -> keV/c^2
extern const double krome_fine_structure_constant; //#
extern const double krome_qe_mass; //9.10938188d-28 !g
extern const double krome_qp_mass; //1.67262158d-24 !g
extern const double krome_qn_mass; //g
extern const double krome_Dalpha; //7.2973525664d-3 !dark fine structure constant
extern const double krome_xi; //Dark CMB temperature to Normal CMB temperature ratio
extern const double krome_iqp_mass; //1/g


extern void krome_set_tcmb(double arg);
extern double krome_get_tcmb();
extern void krome_set_zredshift(double arg);
extern double krome_get_zredshift();
extern void krome_set_orthopararatio(double arg);
extern double krome_get_orthopararatio();
extern void krome_set_metallicity(double arg);
extern double krome_get_metallicity();
extern void krome_set_tfloor(double arg);
extern double krome_get_tfloor();


extern double krome_get_table_tdust(double *x,double *tgas);
extern double krome_num2col(double num, double *x, double tgas);
extern void krome_print_phys_variables();
extern void krome_set_mpi_rank(int rank);
extern void krome_store(double *x, double tgas, double dt);
extern void krome_restore(double *x, double *tgas, double *dt);
extern void krome_thermo_on();
extern void krome_thermo_off();
extern void krome_get_coef(double tgas, double* x, double* krome_get_coef_var);
extern double krome_get_mu_x(double *xin);
extern double krome_get_gamma_x(double *xin, double intgas);
extern void krome_consistent_x(double *x);
extern void krome_n2x(double *n, double rhogas, double* krome_n2x_var);
extern void krome_x2n(double *x, double rhogas, double* krome_x2n_var);
extern void krome_thermo(double *x, double *tgas, double dt);
extern double krome_get_heating(double *x, double intgas);
extern void krome_get_heating_array(double *x, double intgas, double* krome_get_heating_array_var);
extern double krome_get_cooling(double *x, double intgas);
extern void krome_get_cooling_array(double *x, double intgas, double* krome_get_cooling_array_var);
extern void krome_plot_cooling(double *n);
/* here is another example of a fortran subroutine which accepts optional arguments. */
extern void krome_dump_cooling(double *n, double tgas, int nfile_in);
extern void krome_conservelin_x(double *x, double *ref);
extern void krome_conservelingetref_x(double *x, double* krome_conservelingetref_var);
extern void krome_conserve(double *x, double *xi, double* krome_conserve_var);
extern double krome_get_gamma(double *x, double tgas);
extern void krome_get_zatoms(int* krome_get_zatoms_var);
extern double krome_get_mu(double *x);
// extern char **krome_get_rnames();
extern void krome_get_mass(double* krome_get_mass_var);
extern void krome_get_imass(double* krome_get_imass_var);
extern double krome_get_hnuclei(double *x);
extern void krome_get_charges(double* krome_get_charges_var);
// extern char **krome_get_names();
extern int krome_get_index(char *name);
extern double krome_get_rho(double *n);
extern void krome_scale_z(double *n, double z);
extern void krome_set_z(double xarg);
extern void krome_set_dust_to_gas(double xarg);
extern void krome_set_clump(double xarg);
extern double krome_get_electrons(double *x);
extern void krome_print_best_flux(double *xin, double tgas, int nbest);
extern void krome_print_best_flux_frac(double *xin, double tgas, double frac);
extern void krome_print_best_flux_spec(double *xin, double tgas, int nbest, int idx_find);
extern void krome_get_flux(double *n, double tgas, double* krome_get_flux_var);
extern void krome_explore_flux(double *x, double tgas, int ifile, double xvar);
extern void krome_get_qeff(double* krome_get_qeff_var);
extern void krome_dump_flux(double *n, double tgas, int nfile);
extern void krome_dump_rates(double intmin, double intmax, int imax, int funit);
extern void krome_get_info(double *x, double tgas);

