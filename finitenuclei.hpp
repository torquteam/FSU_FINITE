#ifndef finitenuclei_hpp
#define finitenuclei_hpp

#include <stdio.h>
#include <string>
using namespace std;
extern "C" {
double proton_scalardens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless);
double neutron_scalardens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless);
double neutron_vectordens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless);
double proton_vectordens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless);
double neutron_tensordens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless);
double proton_tensordens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless);
double dGdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double Fn_r_unitless, double Gn_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp);
double dFdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double Fn_r_unitless, double Gn_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp);
double dAdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp);
double dBdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp);
void rk4_n(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double FG_unitless[2], double** meson_fields_unitless, int nrows, char direction, double** &Fn_unitless, double** &Gn_unitless, double fw, double fp, int ncols_meson);
void rk4_p(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double AB_unitless[2], double** meson_fields_unitless, int nrows, char direction, double** &Ap_unitless, double** &Bp_unitless, double fw, double fp, int ncols_meson);
double neutronfield_Solve(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_mev, double** &Fn_unitless, double** &Gn_unitless, bool stitch, double r_init_fm, double r_final_fm, double fw, double fp, int ncols_meson);
double protonfield_Solve(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_mev, double** &Ap_unitless, double** &Bp_unitless, bool stitch, double r_init_fm, double r_final_fm, double fw, double fp, int ncols_meson);
void normalize(double** &AF_unitless, double** &BG_unitless, int nrows);
double energy_bisect_n(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Fn_unitless, double** &Gn_unitless, double fw, double fp, int ncols_meson);
double energy_bisect_p(double alpha,double** meson_fields_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Ap_unitless, double** &Bp_unitless, double fw, double fp, int ncols_meson);
void energy_spectrum_neutron(double** meson_fields_unitless, int nrows_meson, int A, double fw, double fp, int ncols_meson);
void energy_spectrum_proton(double** meson_fields_unitless,int nrows_meson, int A, double fw, double fp, int ncols_meson);
int shell_fill(int A, int Z, int en_col, int j_col);
void get_densities(double** meson_fields_unitless, int A, string energy_spectrum_neutron, string energy_spectrum_proton, double** &densities_svnp_unitless, int nrows_meson, int ncols_dens, double fw, double fp, int ncols_meson);
double greens_meson(double r_unitless, double rp_unitless, double meson_mass_mev);
double greens_coulomb(double r_unitless, double rp_unitless);
void get_nonlinear_meson_fields(double** &meson_fields_unitless, int npoints_meson, int A, double** densities, int nrows_dens, int ncols_dens, int sdens_n_col, int vdens_n_col, int sdens_p_col, int vdens_p_col, double gs2, double gw2, double gp2, double gd2, double kappa, double lambda, double zeta, double xi, double lambda_v, double lambda_s, double fw, double fp, double mSigma_mev, double mOmega_mev, double mRho_mev, double mDelta_mev, int gridsize_meson, int meson_iterations, int ncols_meson);
double get_BA(double** meson_field_unitless, double** densities_unitless, string n_energies, string p_energies, int npoints_meson, int npoints_densities, int ncols_density, int A, double kappa, double lambda, double zeta, double xi, double lambda_v, double lambda_s, double fw, double fp);
void get_nucleon_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_N_P_fm2[2]);
void subtract_restmass(string proton_spectrum, string neutron_spectrum);
int hartree_method(double fin_couplings[16], int A, int Z, int iterations, int gridsize, int meson_iterations, double Observables[7], double convergence_help, bool print_densities, bool print_meson_fields);
void EM_formfactors(double GE_pn[2], double GM_pn[2], double q);
void WEAK_formfactors(double WGE_pn[2], double WGM_pn[2], double q, double qwn);
void vector_formfactor(double** densities_svtnp_unitless, double q_unitless, int nrows, double vFF_pn[2]);
void tensor_formfactor(double** densities_svtnp_unitless, double q_unitless, int nrows, double tFF_pn[2]);
double charge_formfactor(double q_unitless, double** densities_svtnp_unitless, int nrows, int Z);
double weak_formfactor(double q_unitless, double** densities_svtnp_unitless, int nrows, int A, int Z, double qwn);
double get_WEAK_CHARGE_densities_v2(double** &densities_svtnp_unitless, int Z, int A, int ncols_dens, int nrows_dens, double qwn);
void get_weak_charge_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_C_W_fm2[2], double qwn);
}
#endif