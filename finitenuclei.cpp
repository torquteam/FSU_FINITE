#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "NumMethods.hpp"

using namespace std;

// Physical constants
const double pi = 4.0*atan(1.0);
const double mN_mev = 939.0; //939.56542052;
const double mP_mev = 939.0; //938.27208816;
const double mNuc_mev = (mN_mev+mP_mev)/2.0;

// conversions
const double e0_e2permevfm = 0.055263; // e^2/(MeV fm) 
const double hbar_mevfm = 197.32698; // MeV fm
const double r0_fm = 1.25; //fm (arbitrary)
const double enscale_mev = pow(hbar_mevfm,2)/(2.0*mNuc_mev*pow(r0_fm,2)); // arbitrary
const double fm_to_inversemev = 1.0/197.32698;
const double e0_unitless = e0_e2permevfm*r0_fm*enscale_mev;
const double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;

const double mNuc_unitless = mNuc_mev/enscale_mev;

// Generic?
//const double qwp = 0.0713;  // weak vector-charge of the proton
//const double qwn = -0.9821; // weak vector-charge of the neutron

// 208Pb
const double qwp = 0.0713;  // weak vector-charge of the proton
const double qwn_Pb = -0.9821; // weak vector-charge of the neutron
const double qwn_Ca = -0.9795; // weak vector-charge of the neutron

const double q_transfer_Ca = 0.8733;
const double q_transfer_Pb = 0.3977;

data2 dm2;
nummeth nm;


// proton scalar density given proton wave functions (unitless)
double proton_scalardens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless) {
    double sdens_p = pfrac*(2.0*jp+1.0)/(4.0*pi*pow(r_unitless,2))*(pow(Ap_r_unitless,2.0) - pow(Bp_r_unitless,2.0));
    return sdens_p;
}

// neutron scalar density given neutron wave functions (unitless)
double neutron_scalardens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless) {
    double sdens_n = nfrac*(2.0*jn+1.0)/(4.0*pi*pow(r_unitless,2.0))*(pow(Fn_r_unitless,2.0) - pow(Gn_r_unitless,2.0));
    return sdens_n;
}

// neutron vector density given neutron wave functions
double neutron_vectordens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless) {
    double vdens_n = nfrac*(2.0*jn+1.0)/(4.0*pi*pow(r_unitless,2.0))*(pow(Fn_r_unitless,2.0) + pow(Gn_r_unitless,2.0));
    return vdens_n;
}

// proton vector density given proton wave functions
double proton_vectordens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless) {
    double vdens_p = pfrac*(2.0*jp+1.0)/(4.0*pi*pow(r_unitless,2))*(pow(Ap_r_unitless,2.0) + pow(Bp_r_unitless,2.0));
    return vdens_p;
}

// neutron tensor density given neutron wave functions
double neutron_tensordens(double jn, double nfrac, double Fn_r_unitless, double Gn_r_unitless, double r_unitless) {
    double tdens_n = 2.0*nfrac*(2.0*jn+1.0)/(4.0*pi*pow(r_unitless,2.0))*(Fn_r_unitless*Gn_r_unitless);
    return tdens_n;
}

// proton tensor density given proton wave functions
double proton_tensordens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless) {
    double tdens_p = 2.0*pfrac*(2.0*jp+1.0)/(4.0*pi*pow(r_unitless,2))*(Ap_r_unitless*Bp_r_unitless);
    return tdens_p;
}

// obtain a grid of the initial meson field (r,....) (ouput is unitless)
void init_meson_r(int npoints, double** &array, double r_init_fm, double r_final_fm) {
    double r_init_unitless = r_init_fm/r0_fm;
    double r_final_unitless = r_final_fm/r0_fm;

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;

    for (int i=0; i<npoints; ++i) {
        array[i][0] = r_unitless;
        r_unitless = r_unitless+stp;
    }
}

// obtain a grid of the initial meson fields (ouput is unitless)
void init_meson(int npoints, double V_mev, double R_fm, double a_fm, double** &array, double r_init_fm, double r_final_fm, double coupling, int col) {
    double a_unitless = a_fm/r0_fm;
    double R_unitless = R_fm/r0_fm;
    double r_init_unitless = r_init_fm/r0_fm;
    double r_final_unitless = r_final_fm/r0_fm;
    double V_unitless = V_mev/enscale_mev;

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;

    V_unitless = V_unitless*sqrt(coupling);

    for (int i=0; i<npoints; ++i) {
        array[i][col] = V_unitless/(1.0+exp((r_unitless-R_unitless)/a_unitless));
        r_unitless = r_unitless+stp;
    }
}

// obtain a grid of the derivatives of the meson fields (2/r + d/dr)
void divergence_der(int npoints, double** &array, int ref_col, int col) {

    double dydx1;

    // get the derivative
    array[0][col] = (array[2][ref_col] - array[1][ref_col])/(array[2][0] - array[1][0]);
    array[1][col] = (array[2][ref_col] - array[1][ref_col])/(array[2][0] - array[1][0]);
    for (int i=2; i<(npoints-2); ++i) {
        dydx1 = (array[i+1][ref_col] - array[i-1][ref_col])/(array[i+1][0] - array[i-1][0]);
        array[i][col] = dydx1;
    }
    array[npoints-2][col] = (array[npoints-1][ref_col] - array[npoints-3][ref_col])/(array[npoints-1][0] - array[npoints-3][0]);
    array[npoints-1][col] = (array[npoints-1][ref_col] - array[npoints-2][ref_col])/(array[npoints-1][0] - array[npoints-2][0]);
    
    array[0][col] = array[0][col] + (2.0/array[1][0])*array[1][ref_col];
    for (int i=1; i<npoints; ++i) {
        array[i][col] = array[i][col] + (2.0/array[i][0])*array[i][ref_col];
    }
    
}

// obtain a grid of the derivatives of the meson fields
void meson_der(int npoints, double** &array, int ref_col, int col) {

    double dydx1, dydx2;

    // get the derivative
    array[0][col] = (array[1][ref_col] - array[0][ref_col])/(array[1][0] - array[0][0]);
    array[1][col] = (array[2][ref_col] - array[0][ref_col])/(array[2][0] - array[0][0]);
    for (int i=2; i<(npoints-2); ++i) {
        dydx1 = (array[i+1][ref_col] - array[i-1][ref_col])/(array[i+1][0] - array[i-1][0]);
        dydx2 = (array[i+2][ref_col] - array[i-2][ref_col])/(array[i+2][0] - array[i-2][0]);
        array[i][col] = 0.5*(dydx1 + dydx2);
    }
    array[npoints-2][col] = (array[npoints-1][ref_col] - array[npoints-3][ref_col])/(array[npoints-1][0] - array[npoints-3][0]);
    array[npoints-1][col] = (array[npoints-1][ref_col] - array[npoints-2][ref_col])/(array[npoints-1][0] - array[npoints-2][0]);

}


// obtain a grid of the initial coulomb field (r,V(r)) (ouput is unitless)
void init_coulomb(int npoints, double R_fm, double** &array, double r_init_fm, double r_final_fm, int Z, int col) {
    double R_unitless = R_fm/r0_fm;   // convert unitless
    double r_init_unitless = r_init_fm/r0_fm;     // convert unitless
    double r_final_unitless = r_final_fm/r0_fm;   // convert unitless

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;
    int i=0;

    // fill grid for coulomb field
    while (r_unitless<R_unitless) {
        array[i][col] = Z/(4.0*pi*e0_unitless)/(2.0*pow(R_unitless,3))*(3.0*pow(R_unitless,2) - pow(r_unitless,2));
        i = i+1;
        r_unitless = r_unitless+stp;
    }
    while (i<npoints) {
        array[i][col] = Z/(4.0*pi*e0_unitless)*1.0/(r_unitless);
        i = i+1;
        r_unitless = r_unitless+stp;
    }
}
// input is unitless (r/r0, en/Econv)
double dGdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double Fn_r_unitless, double Gn_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp) {
    double res;
    double S_r = gs_sigma_r_unitless - 0.5*gd_delta_r_unitless;
    double V_r = gw_omega_r_unitless - 0.5*gp_rho_r_unitless;
    double T_r = fw/(2.0*mNuc_unitless)*dgw_omega_dr_r_unitless - fp/(4.0*mNuc_unitless)*dgp_rho_dr_r_unitless;
    
    res = -conv_r0_en*(En_unitless - mNuc_unitless + S_r - V_r)*Fn_r_unitless - alpha/r_unitless*Gn_r_unitless + T_r*Gn_r_unitless;
    return res;
}

// input is unitless (r/r0, en/Econv)
double dFdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double Fn_r_unitless, double Gn_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp) {
    double res;
    double S_r = gs_sigma_r_unitless - 0.5*gd_delta_r_unitless;
    double V_r = gw_omega_r_unitless - 0.5*gp_rho_r_unitless;
    double T_r = fw/(2.0*mNuc_unitless)*dgw_omega_dr_r_unitless - fp/(4.0*mNuc_unitless)*dgp_rho_dr_r_unitless;
    
    res = conv_r0_en*(En_unitless + mNuc_unitless - S_r - V_r)*Gn_r_unitless + alpha/r_unitless*Fn_r_unitless - T_r*Fn_r_unitless;
    return res;
}

double dBdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp) {
    double res;
    double S_r = gs_sigma_r_unitless + 0.5*gd_delta_r_unitless;
    double V_r = gw_omega_r_unitless + 0.5*gp_rho_r_unitless + e_coulomb_r_unitless;
    double T_r = fw/(2.0*mNuc_unitless)*dgw_omega_dr_r_unitless + fp/(4.0*mNuc_unitless)*dgp_rho_dr_r_unitless;
    
    res = -conv_r0_en*(Ep_unitless - mNuc_unitless + S_r - V_r)*Ap_r_unitless - alpha/r_unitless*Bp_r_unitless + T_r*Bp_r_unitless;
    return res;
}

// input is unitless (r/r0, en/Econv)
double dAdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double gd_delta_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless, double dgw_omega_dr_r_unitless, double dgp_rho_dr_r_unitless, double fw, double fp) {
    double res;
    double S_r = gs_sigma_r_unitless + 0.5*gd_delta_r_unitless;
    double V_r = gw_omega_r_unitless + 0.5*gp_rho_r_unitless + e_coulomb_r_unitless;
    double T_r = fw/(2.0*mNuc_unitless)*dgw_omega_dr_r_unitless + fp/(4.0*mNuc_unitless)*dgp_rho_dr_r_unitless;

    res = conv_r0_en*(Ep_unitless + mNuc_unitless - S_r - V_r)*Bp_r_unitless + alpha/r_unitless*Ap_r_unitless - T_r*Ap_r_unitless;
    return res;
}

void rk4_n(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double FG_unitless[2], double** meson_fields_unitless, int nrows, char direction, double** &Fn_unitless, double** &Gn_unitless, double fw, double fp, int ncols_meson) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    double gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless;
    double h = (meson_fields_unitless[nrows-1][0]-meson_fields_unitless[0][0])/(nrows-1);

    // need initial conditions for Fn_r, Gn_r
    double Fn_r_unitless, Gn_r_unitless;
    double r_unitless = r_init_unitless;
    double r_mid;
    double* meson_mid_array;
    
    double mstar = mNuc_unitless - meson_fields_unitless[0][1] + 0.5*meson_fields_unitless[0][4];
    double estar = en_unitless - meson_fields_unitless[0][2] + 0.5*meson_fields_unitless[0][3];
    if (direction == 'r') {
        if (alpha>0) {
            Fn_r_unitless = 1e-8;
            Gn_r_unitless = Fn_r_unitless*r_unitless*(mstar-estar)/(2.0*alpha+1.0);
        } else {
            Gn_r_unitless = 1e-8;
            Fn_r_unitless = Gn_r_unitless*r_unitless*(mstar+estar)/(1.0-2.0*alpha);
        }
        Fn_unitless[0][0] = r_unitless; Fn_unitless[0][1] = Fn_r_unitless;
        Gn_unitless[0][0] = r_unitless; Gn_unitless[0][1] = Gn_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            r_mid = r_unitless+h/2.0;
            dm2.interpolate_multi(nrows,ncols_meson,meson_fields_unitless,r_mid,0,meson_mid_array,true);
            gs_sigma_r_mid_unitless = meson_mid_array[0];
            gw_omega_r_mid_unitless = meson_mid_array[1];
            gp_rho_r_mid_unitless = meson_mid_array[2];
            gd_delta_r_mid_unitless = meson_mid_array[3];
            dgw_omega_dr_mid_unitless = meson_mid_array[5];
            dgp_rho_dr_mid_unitless = meson_mid_array[6];
            k1 = h*dFdr(r_unitless, alpha, en_unitless, meson_fields_unitless[i-1][1], meson_fields_unitless[i-1][2], meson_fields_unitless[i-1][3], meson_fields_unitless[i-1][4], Fn_r_unitless, Gn_r_unitless, meson_fields_unitless[i-1][6], meson_fields_unitless[i-1][7], fw,fp);
            l1 = h*dGdr(r_unitless, alpha, en_unitless, meson_fields_unitless[i-1][1], meson_fields_unitless[i-1][2], meson_fields_unitless[i-1][3], meson_fields_unitless[i-1][4], Fn_r_unitless, Gn_r_unitless, meson_fields_unitless[i-1][6], meson_fields_unitless[i-1][7], fw,fp);
            k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            k4 = h*dFdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[i][1], meson_fields_unitless[i][2], meson_fields_unitless[i][3], meson_fields_unitless[i][4], Fn_r_unitless+k3, Gn_r_unitless+l3, meson_fields_unitless[i][6], meson_fields_unitless[i][7], fw,fp);
            l4 = h*dGdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[i][1], meson_fields_unitless[i][2], meson_fields_unitless[i][3], meson_fields_unitless[i][4], Fn_r_unitless+k3, Gn_r_unitless+l3, meson_fields_unitless[i][6], meson_fields_unitless[i][7], fw,fp);
            Fn_r_unitless = Fn_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Gn_r_unitless = Gn_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = meson_fields_unitless[i][0];
            Fn_unitless[i][0] = r_unitless; Fn_unitless[i][1] = Fn_r_unitless;
            Gn_unitless[i][0] = r_unitless; Gn_unitless[i][1] = Gn_r_unitless;
            delete meson_mid_array;
        }
    } else {
        mstar = mNuc_unitless - meson_fields_unitless[nrows-1][1] + 0.5*meson_fields_unitless[nrows-1][4];
        estar = en_unitless - meson_fields_unitless[nrows-1][2] + 0.5*meson_fields_unitless[nrows-1][3];
        double w = sqrt(mstar*mstar-estar*estar);
        Fn_r_unitless = 1e-30;
        Gn_r_unitless = -Fn_r_unitless*(conv_r0_en*w + alpha/r_unitless)/(estar+mstar);
        h = -h;
        Fn_unitless[0][0] = r_unitless; Fn_unitless[0][1] = Fn_r_unitless;
        Gn_unitless[0][0] = r_unitless; Gn_unitless[0][1] = Gn_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            r_mid = r_unitless+h/2.0;
            dm2.interpolate_multi(nrows,ncols_meson,meson_fields_unitless,r_mid,0,meson_mid_array,true);
            gs_sigma_r_mid_unitless = meson_mid_array[0];
            gw_omega_r_mid_unitless = meson_mid_array[1];
            gp_rho_r_mid_unitless = meson_mid_array[2];
            gd_delta_r_mid_unitless = meson_mid_array[3];
            dgw_omega_dr_mid_unitless = meson_mid_array[5];
            dgp_rho_dr_mid_unitless = meson_mid_array[6];
            k1 = h*dFdr(r_unitless, alpha, en_unitless, meson_fields_unitless[nrows-i][1], meson_fields_unitless[nrows-i][2], meson_fields_unitless[nrows-i][3], meson_fields_unitless[nrows-i][4], Fn_r_unitless, Gn_r_unitless, meson_fields_unitless[nrows-i][6], meson_fields_unitless[nrows-i][7], fw,fp);
            l1 = h*dGdr(r_unitless, alpha, en_unitless, meson_fields_unitless[nrows-i][1], meson_fields_unitless[nrows-i][2], meson_fields_unitless[nrows-i][3], meson_fields_unitless[nrows-i][4], Fn_r_unitless, Gn_r_unitless, meson_fields_unitless[nrows-i][6], meson_fields_unitless[nrows-i][7], fw,fp);
            k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless, fw,fp);
            k4 = h*dFdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[nrows-i-1][1], meson_fields_unitless[nrows-i-1][2], meson_fields_unitless[nrows-i-1][3], meson_fields_unitless[nrows-i-1][4], Fn_r_unitless+k3, Gn_r_unitless+l3, meson_fields_unitless[nrows-i-1][6], meson_fields_unitless[nrows-i-1][7], fw,fp);
            l4 = h*dGdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[nrows-i-1][1], meson_fields_unitless[nrows-i-1][2], meson_fields_unitless[nrows-i-1][3], meson_fields_unitless[nrows-i-1][4], Fn_r_unitless+k3, Gn_r_unitless+l3, meson_fields_unitless[nrows-i-1][6], meson_fields_unitless[nrows-i-1][7], fw,fp);
            Fn_r_unitless = Fn_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Gn_r_unitless = Gn_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = meson_fields_unitless[nrows-1-i][0];
            Fn_unitless[i][0] = r_unitless; Fn_unitless[i][1] = Fn_r_unitless;
            Gn_unitless[i][0] = r_unitless; Gn_unitless[i][1] = Gn_r_unitless;
            delete meson_mid_array;
        }
    }

    FG_unitless[0] = Fn_r_unitless;
    FG_unitless[1] = Gn_r_unitless;
}

void rk4_p(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double AB_unitless[2], double** meson_fields_unitless, int nrows, char direction, double** &Ap_unitless, double** &Bp_unitless, double fw, double fp, int ncols_meson) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    double gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless;
    double h = (meson_fields_unitless[nrows-1][0]-meson_fields_unitless[0][0])/(nrows-1);

    // need initial conditions for Ap_r, Bp_r
    double Ap_r_unitless, Bp_r_unitless;
    double r_unitless = r_init_unitless;
    double r_mid;
    double* meson_mid_array;
    
    double mstar = mNuc_unitless - meson_fields_unitless[0][1] - 0.5*meson_fields_unitless[0][4];
    double estar = en_unitless - meson_fields_unitless[0][2] - 0.5*meson_fields_unitless[0][3];
    if (direction == 'r') {
        if (alpha>0) {
            Ap_r_unitless = 1e-8;
            Bp_r_unitless = Ap_r_unitless*r_unitless*(mstar-estar)/(2.0*alpha+1.0);
        } else {
            Bp_r_unitless = 1e-8;
            Ap_r_unitless = Bp_r_unitless*r_unitless*(mstar+estar)/(1.0-2.0*alpha);
        }
        Ap_unitless[0][0] = r_unitless; Ap_unitless[0][1] = Ap_r_unitless;
        Bp_unitless[0][0] = r_unitless; Bp_unitless[0][1] = Bp_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            r_mid = r_unitless+h/2.0;
            dm2.interpolate_multi(nrows,ncols_meson,meson_fields_unitless,r_mid,0,meson_mid_array,true);
            gs_sigma_r_mid_unitless = meson_mid_array[0];
            gw_omega_r_mid_unitless = meson_mid_array[1];
            gp_rho_r_mid_unitless = meson_mid_array[2];
            gd_delta_r_mid_unitless = meson_mid_array[3];
            e_coulomb_r_mid_unitless = meson_mid_array[4];
            dgw_omega_dr_mid_unitless = meson_mid_array[5];
            dgp_rho_dr_mid_unitless = meson_mid_array[6];
            k1 = h*dAdr(r_unitless, alpha, en_unitless, meson_fields_unitless[i-1][1], meson_fields_unitless[i-1][2], meson_fields_unitless[i-1][3], meson_fields_unitless[i-1][4], meson_fields_unitless[i-1][5], Ap_r_unitless, Bp_r_unitless, meson_fields_unitless[i-1][6], meson_fields_unitless[i-1][7], fw,fp);
            l1 = h*dBdr(r_unitless, alpha, en_unitless, meson_fields_unitless[i-1][1], meson_fields_unitless[i-1][2], meson_fields_unitless[i-1][3], meson_fields_unitless[i-1][4], meson_fields_unitless[i-1][5], Ap_r_unitless, Bp_r_unitless, meson_fields_unitless[i-1][6], meson_fields_unitless[i-1][7], fw,fp);
            k2 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            l2 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            k3 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            l3 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            k4 = h*dAdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[i][1], meson_fields_unitless[i][2], meson_fields_unitless[i][3], meson_fields_unitless[i][4], meson_fields_unitless[i][5], Ap_r_unitless+k3, Bp_r_unitless+l3, meson_fields_unitless[i][6], meson_fields_unitless[i][7],fw,fp);
            l4 = h*dBdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[i][1], meson_fields_unitless[i][2], meson_fields_unitless[i][3], meson_fields_unitless[i][4], meson_fields_unitless[i][5], Ap_r_unitless+k3, Bp_r_unitless+l3, meson_fields_unitless[i][6], meson_fields_unitless[i][7],fw,fp);
            Ap_r_unitless = Ap_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Bp_r_unitless = Bp_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = meson_fields_unitless[i][0];
            Ap_unitless[i][0] = r_unitless; Ap_unitless[i][1] = Ap_r_unitless;
            Bp_unitless[i][0] = r_unitless; Bp_unitless[i][1] = Bp_r_unitless;
            delete meson_mid_array;
        }
    } else {
        mstar = mNuc_unitless - meson_fields_unitless[nrows-1][1] - 0.5*meson_fields_unitless[nrows-1][4];
        estar = en_unitless - meson_fields_unitless[nrows-1][2] - 0.5*meson_fields_unitless[nrows-1][3];
        double w = sqrt(mstar*mstar-estar*estar);
        Ap_r_unitless = 1e-30;
        Bp_r_unitless = -Ap_r_unitless*(conv_r0_en*w + alpha/r_unitless)/(estar+mstar);
        h = -h;
        Ap_unitless[0][0] = r_unitless; Ap_unitless[0][1] = Ap_r_unitless;
        Bp_unitless[0][0] = r_unitless; Bp_unitless[0][1] = Bp_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            r_mid = r_unitless+h/2.0;
            dm2.interpolate_multi(nrows,ncols_meson,meson_fields_unitless,r_mid,0,meson_mid_array,true);
            gs_sigma_r_mid_unitless = meson_mid_array[0];
            gw_omega_r_mid_unitless = meson_mid_array[1];
            gp_rho_r_mid_unitless = meson_mid_array[2];
            gd_delta_r_mid_unitless = meson_mid_array[3];
            e_coulomb_r_mid_unitless = meson_mid_array[4];
            dgw_omega_dr_mid_unitless = meson_mid_array[5];
            dgp_rho_dr_mid_unitless = meson_mid_array[6];
            k1 = h*dAdr(r_unitless, alpha, en_unitless, meson_fields_unitless[nrows-i][1], meson_fields_unitless[nrows-i][2], meson_fields_unitless[nrows-i][3], meson_fields_unitless[nrows-i][4], meson_fields_unitless[nrows-i][5], Ap_r_unitless, Bp_r_unitless, meson_fields_unitless[nrows-i][6], meson_fields_unitless[nrows-i][7], fw,fp);
            l1 = h*dBdr(r_unitless, alpha, en_unitless, meson_fields_unitless[nrows-i][1], meson_fields_unitless[nrows-i][2], meson_fields_unitless[nrows-i][3], meson_fields_unitless[nrows-i][4], meson_fields_unitless[nrows-i][5], Ap_r_unitless, Bp_r_unitless, meson_fields_unitless[nrows-i][6], meson_fields_unitless[nrows-i][7], fw,fp);
            k2 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            l2 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            k3 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            l3 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_mid_unitless, gw_omega_r_mid_unitless, gp_rho_r_mid_unitless, gd_delta_r_mid_unitless, e_coulomb_r_mid_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0, dgw_omega_dr_mid_unitless, dgp_rho_dr_mid_unitless,fw,fp);
            k4 = h*dAdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[nrows-i-1][1], meson_fields_unitless[nrows-i-1][2], meson_fields_unitless[nrows-i-1][3], meson_fields_unitless[nrows-i-1][4], meson_fields_unitless[nrows-i-1][5], Ap_r_unitless+k3, Bp_r_unitless+l3, meson_fields_unitless[nrows-i-1][6], meson_fields_unitless[nrows-i-1][7],fw,fp);
            l4 = h*dBdr(r_unitless+h, alpha, en_unitless, meson_fields_unitless[nrows-i-1][1], meson_fields_unitless[nrows-i-1][2], meson_fields_unitless[nrows-i-1][3], meson_fields_unitless[nrows-i-1][4], meson_fields_unitless[nrows-i-1][5], Ap_r_unitless+k3, Bp_r_unitless+l3, meson_fields_unitless[nrows-i-1][6], meson_fields_unitless[nrows-i-1][7],fw,fp);
            Ap_r_unitless = Ap_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Bp_r_unitless = Bp_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = meson_fields_unitless[nrows-1-i][0];
            Ap_unitless[i][0] = r_unitless; Ap_unitless[i][1] = Ap_r_unitless;
            Bp_unitless[i][0] = r_unitless; Bp_unitless[i][1] = Bp_r_unitless;
            delete meson_mid_array;
        }
    }

    AB_unitless[0] = Ap_r_unitless;
    AB_unitless[1] = Bp_r_unitless;
}

// Compute wave functions for given energy and alpha (returns the determinant: -G_ext_unitless*F_int_unitless + F_ext_unitless*G_int_unitless) Fn and Gn must have dimensions (npoints_rk4*2 , 2)
double neutronfield_Solve(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_mev, double** &Fn_unitless, double** &Gn_unitless, bool stitch, double r_init_fm, double r_final_fm, double fw, double fp, int ncols_meson) {

    // Integration parameters
    double r_left_fm = r_init_fm;
    double r_right_fm = r_final_fm;
    double r_match_fm = 0.6*r0_fm*pow(A,1.0/3.0); // fm
    double FG_unitless[2];
    char direction;

    // convert unitless
    double r_left_unitless = r_left_fm/r0_fm;
    double r_right_unitless = r_right_fm/r0_fm;
    double r_match_unitless = r_match_fm/r0_fm;
    double energy_unitless = en_mev/enscale_mev;

    // create temporary arrays for interior and exterior solutions
    int r_match_row = dm2.findvalue(meson_fields_unitless,nrows_meson,ncols_meson,r_match_unitless,0,1e-2); // fm
    int nrows_int = r_match_row + 1;
    int nrows_ext = nrows_meson - nrows_int + 1;
    double** Fn_int_unitless; double** Fn_ext_unitless;
    double** Gn_int_unitless; double** Gn_ext_unitless;
    dm2.create(Fn_int_unitless,nrows_int,2);
    dm2.create(Fn_ext_unitless,nrows_ext,2);
    dm2.create(Gn_int_unitless,nrows_int,2);
    dm2.create(Gn_ext_unitless,nrows_ext,2);
    
    double F_int_unitless_rmatch; double G_int_unitless_rmatch;
    double F_ext_unitless_rmatch; double G_ext_unitless_rmatch;

    // integrate from origin to r_match
    direction = 'r';
    rk4_n(r_left_unitless, meson_fields_unitless[r_match_row][0], nrows_int, alpha, energy_unitless, FG_unitless, meson_fields_unitless, nrows_meson, direction, Fn_int_unitless, Gn_int_unitless, fw, fp, ncols_meson);
    F_int_unitless_rmatch = FG_unitless[0];
    G_int_unitless_rmatch = FG_unitless[1];

    // integrate from infinity to r_match
    direction = 'l';
    rk4_n(r_right_unitless, meson_fields_unitless[r_match_row][0], nrows_ext, alpha, energy_unitless, FG_unitless, meson_fields_unitless, nrows_meson, direction, Fn_ext_unitless, Gn_ext_unitless, fw, fp, ncols_meson);
    F_ext_unitless_rmatch = FG_unitless[0];
    G_ext_unitless_rmatch = FG_unitless[1];
    
    if (stitch==true) {
        // stitch the interior and exterior solutions together
        dm2.invert_order(Fn_ext_unitless,nrows_ext,2);
        dm2.invert_order(Gn_ext_unitless,nrows_ext,2);
        for (int i=0; i<nrows_int; ++i) {
            Fn_unitless[i][0] = Fn_int_unitless[i][0]; Fn_unitless[i][1] = Fn_int_unitless[i][1];
            Gn_unitless[i][0] = Gn_int_unitless[i][0]; Gn_unitless[i][1] = Gn_int_unitless[i][1];
        }
        for (int i=0; i<(nrows_ext-1); ++i) {
            Fn_unitless[i+nrows_int][0] = Fn_ext_unitless[i+1][0]; Fn_unitless[i+nrows_int][1] = (F_int_unitless_rmatch/F_ext_unitless_rmatch)*Fn_ext_unitless[i+1][1];
            Gn_unitless[i+nrows_int][0] = Gn_ext_unitless[i+1][0]; Gn_unitless[i+nrows_int][1] = (G_int_unitless_rmatch/G_ext_unitless_rmatch)*Gn_ext_unitless[i+1][1];
        }
        dm2.cleanup(Fn_int_unitless,nrows_int);
        dm2.cleanup(Fn_ext_unitless,nrows_ext);
        dm2.cleanup(Gn_int_unitless,nrows_int);
        dm2.cleanup(Gn_ext_unitless,nrows_ext);

        return -G_ext_unitless_rmatch*F_int_unitless_rmatch + F_ext_unitless_rmatch*G_int_unitless_rmatch;
    } else {
        //dm2.print(Gn_int_unitless,npoints_rk4,2,true,"left.txt");
        //dm2.print(Gn_ext_unitless,npoints_rk4,2,true,"right.txt");
        //cout << (F_int_unitless_rmatch/F_ext_unitless_rmatch) - (G_int_unitless_rmatch/G_ext_unitless_rmatch) << endl;
        dm2.cleanup(Fn_int_unitless,nrows_int);
        dm2.cleanup(Fn_ext_unitless,nrows_ext);
        dm2.cleanup(Gn_int_unitless,nrows_int);
        dm2.cleanup(Gn_ext_unitless,nrows_ext);
        return -G_ext_unitless_rmatch*F_int_unitless_rmatch + F_ext_unitless_rmatch*G_int_unitless_rmatch;
    }
}

// Compute wave functions for given energy and alpha (returns the determinant: -G_ext_unitless*F_int_unitless + F_ext_unitless*G_int_unitless) Ap and Bp must have dimensions (npoints_rk4*2 , 2)
double protonfield_Solve(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_mev, double** &Ap_unitless, double** &Bp_unitless, bool stitch, double r_init_fm, double r_final_fm, double fw, double fp, int ncols_meson) {

    // Integration parameters
    double r_left_fm = r_init_fm;
    double r_right_fm = r_final_fm;
    double r_match_fm = 0.6*r0_fm*pow(A,1.0/3.0); // fm
    double AB_unitless[2];
    char direction;

    // convert unitless
    double r_left_unitless = r_left_fm/r0_fm;
    double r_right_unitless = r_right_fm/r0_fm;
    double r_match_unitless = r_match_fm/r0_fm;
    double energy_unitless = en_mev/enscale_mev;

    // create temporary arrays for interior and exterior solutions
    int r_match_row = dm2.findvalue(meson_fields_unitless,nrows_meson,ncols_meson,r_match_unitless,0,1e-2); // fm
    int nrows_int = r_match_row + 1;
    int nrows_ext = nrows_meson - nrows_int + 1;
    double** Ap_int_unitless; double** Ap_ext_unitless;
    double** Bp_int_unitless; double** Bp_ext_unitless;
    dm2.create(Ap_int_unitless,nrows_int,2);
    dm2.create(Ap_ext_unitless,nrows_ext,2);
    dm2.create(Bp_int_unitless,nrows_int,2);
    dm2.create(Bp_ext_unitless,nrows_ext,2);

    double A_int_unitless_rmatch; double B_int_unitless_rmatch;
    double A_ext_unitless_rmatch; double B_ext_unitless_rmatch;

    // integrate from origin to r_match
    direction = 'r';
    rk4_p(r_left_unitless, r_match_unitless, nrows_int, alpha, energy_unitless, AB_unitless, meson_fields_unitless, nrows_meson, direction, Ap_int_unitless, Bp_int_unitless, fw, fp, ncols_meson);
    A_int_unitless_rmatch = AB_unitless[0];
    B_int_unitless_rmatch = AB_unitless[1];

    // integrate from infinity to r_match
    direction = 'l';
    rk4_p(r_right_unitless, r_match_unitless, nrows_ext, alpha, energy_unitless, AB_unitless, meson_fields_unitless, nrows_meson, direction, Ap_ext_unitless, Bp_ext_unitless, fw, fp, ncols_meson);
    A_ext_unitless_rmatch = AB_unitless[0];
    B_ext_unitless_rmatch = AB_unitless[1];

    if (stitch==true) {
        // stitch the interior and exterior solutions together
        dm2.invert_order(Ap_ext_unitless,nrows_ext,2);
        dm2.invert_order(Bp_ext_unitless,nrows_ext,2);
        for (int i=0; i<nrows_int; ++i) {
            Ap_unitless[i][0] = Ap_int_unitless[i][0]; Ap_unitless[i][1] = Ap_int_unitless[i][1];
            Bp_unitless[i][0] = Bp_int_unitless[i][0]; Bp_unitless[i][1] = Bp_int_unitless[i][1];
        }
        for (int i=0; i<(nrows_ext-1); ++i) {
            Ap_unitless[i+nrows_int][0] = Ap_ext_unitless[i+1][0]; Ap_unitless[i+nrows_int][1] = (A_int_unitless_rmatch/A_ext_unitless_rmatch)*Ap_ext_unitless[i+1][1];
            Bp_unitless[i+nrows_int][0] = Bp_ext_unitless[i+1][0]; Bp_unitless[i+nrows_int][1] = (B_int_unitless_rmatch/B_ext_unitless_rmatch)*Bp_ext_unitless[i+1][1];
        }
        dm2.cleanup(Ap_int_unitless,nrows_int);
        dm2.cleanup(Ap_ext_unitless,nrows_ext);
        dm2.cleanup(Bp_int_unitless,nrows_int);
        dm2.cleanup(Bp_ext_unitless,nrows_ext);

        return -B_ext_unitless_rmatch*A_int_unitless_rmatch + A_ext_unitless_rmatch*B_int_unitless_rmatch;
    } else {
        //dm2.print(Gn_int_unitless,npoints_rk4,2,true,"left.txt");
        //dm2.print(Gn_ext_unitless,npoints_rk4,2,true,"right.txt");
        //cout << (F_int_unitless_rmatch/F_ext_unitless_rmatch) - (G_int_unitless_rmatch/G_ext_unitless_rmatch) << endl;
        dm2.cleanup(Ap_int_unitless,nrows_int);
        dm2.cleanup(Ap_ext_unitless,nrows_ext);
        dm2.cleanup(Bp_int_unitless,nrows_int);
        dm2.cleanup(Bp_ext_unitless,nrows_ext);
        return -B_ext_unitless_rmatch*A_int_unitless_rmatch + A_ext_unitless_rmatch*B_int_unitless_rmatch;
    }
}

// normalize the proton or neutron wave functions
void normalize(double** &AF_unitless, double** &BG_unitless, int nrows) {

    // range and step size to integrate
    double r_init_unitless = AF_unitless[0][0];
    double r_final_unitless = AF_unitless[nrows-1][0];
    double h = (r_final_unitless-r_init_unitless)/(nrows-1);

    // initial conditions
    double r_unitless = r_init_unitless;
    double psi2 = 0;
    double integrand,integrand_mid,integrand_next,r_mid,AF_mid_unitless,BG_mid_unitless;

    // integrate psi^2 
    for (int i=0; i<(nrows-1); ++i) {
        r_unitless = AF_unitless[i][0];
        r_mid = r_unitless + h/2.0;
        AF_mid_unitless = dm2.interpolate(nrows,2,AF_unitless,r_mid,0,1,true);
        BG_mid_unitless = dm2.interpolate(nrows,2,BG_unitless,r_mid,0,1,true);
        integrand = pow(AF_unitless[i][1],2.0) + pow(BG_unitless[i][1],2.0);
        integrand_mid = pow(AF_mid_unitless,2.0) + pow(BG_mid_unitless,2.0);
        integrand_next = pow(AF_unitless[i+1][1],2.0) + pow(BG_unitless[i+1][1],2.0);
        psi2 = psi2 + h/6.0*(integrand + 4.0*integrand_mid + integrand_next);
    }
    double norm = sqrt(1.0/psi2); // normalization factor 

    // normalize the wave functions
    for (int i=0; i<nrows; ++i) {
        AF_unitless[i][1] = norm*AF_unitless[i][1];
        BG_unitless[i][1] = norm*BG_unitless[i][1];
    }
}

// returns the energy solution for a given range where the en_determinant changes sign
double energy_bisect_n(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Fn_unitless, double** &Gn_unitless, double fw, double fp, int ncols_meson) {
    int count = 0;
    double midx_mev = 0.0;
    double sol_mev = 0.0;
    double en_determinant;
    bool stitch = false;
    double r_init_fm = meson_fields_unitless[0][0]*r0_fm;
    double r_final_fm = meson_fields_unitless[nrows_meson-1][0]*r0_fm;
    double min_en_error = fabs(neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson));
    double max_en_error = fabs(neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson));
    double thresh = 1e-7;
    double error = thresh*min(min_en_error,max_en_error);
    double midy = error*2;

    // bisect the given range until convergence or too many bisections
    while (fabs(midy)>error) {
        count = count + 1;
        midx_mev = (en_min_mev+en_max_mev)/2.0;
        en_determinant = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        midy = en_determinant;
        en_determinant = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);

        if ((midy*en_determinant)<0) {
            en_max_mev = midx_mev;
        } else {
            en_min_mev = midx_mev;
        }
        sol_mev = midx_mev;

        // Check for divergence
        if (count>75 && fabs(midy)>thresh) {
            cout << "No en found for en range: (" << en_min_mev << "," << en_max_mev << ")" << " error= " << fabs(midy) << endl;
            exit(0);
        } else if (count>75 && fabs(midy)<thresh) {
            stitch=true;
            neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            normalize(Fn_unitless,Gn_unitless,nrows_meson);
            return sol_mev;
        }
    }
    
    stitch=true;
    neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
    normalize(Fn_unitless,Gn_unitless,nrows_meson);
    return sol_mev;
}

// returns the energy solution for a given range where the en_determinant changes sign
double energy_bisect_p(double alpha, double** meson_fields_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Ap_unitless, double** &Bp_unitless, double fw, double fp, int ncols_meson) {
    int count = 0;
    double midx_mev = 0.0;
    double sol_mev = 0.0;
    double en_determinant;
    bool stitch = false;
    double r_init_fm = meson_fields_unitless[0][0]*r0_fm;
    double r_final_fm = meson_fields_unitless[nrows_meson-1][0]*r0_fm;
    double min_en_error = fabs(protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_min_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson));
    double max_en_error = fabs(protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_max_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson));
    double thresh = 1e-7;
    double error = thresh*min(min_en_error,max_en_error);
    double midy = error*2;

    // bisect the given range until convergence or too many bisections
    while (fabs(midy)>error) {
        count = count + 1;
        midx_mev = (en_min_mev+en_max_mev)/2.0;
        en_determinant = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        midy = en_determinant;
        en_determinant = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_min_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);

        if ((midy*en_determinant)<0) {
            en_max_mev = midx_mev;
        } else {
            en_min_mev = midx_mev;
        }
        sol_mev = midx_mev;
        //cout << scientific << setprecision(10) << count << "  " << midx_mev << "  " << midy << endl;
        // Check for divergence
        if (count>75 && fabs(midy)>thresh) {
            cout << "No en found for en range: (" << en_min_mev << "," << en_max_mev << ")" << " error= " << fabs(midy) << endl;
            exit(0);
        } else if (count>75 && fabs(midy)<thresh) {
            stitch=true;
            cout << "flag" << endl;
            protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            normalize(Ap_unitless,Bp_unitless,nrows_meson);
            return sol_mev;
        }
    }
    
    stitch=true;
    protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm,fw,fp,ncols_meson);
    normalize(Ap_unitless,Bp_unitless,nrows_meson);
    return sol_mev;
}

void energy_spectrum_neutron(double** meson_fields_unitless, int nrows_meson, int A, double fw, double fp, int ncols_meson) {
    double en_n_mev;
    double en_n_determinant, en_n_determinant_prev, bound_state_mev;
    double alpha, j;
    double h = 4.0; // energy step
    int node = 0;   // set nodes to zero

    double r_init_fm = meson_fields_unitless[0][0]*r0_fm;
    double r_final_fm = meson_fields_unitless[nrows_meson-1][0]*r0_fm;

    // create the nucleon arrays
    double** Fn_unitless; double** Gn_unitless;
    dm2.create(Fn_unitless,nrows_meson,2);
    dm2.create(Gn_unitless,nrows_meson,2);

    // print out the neutron energy spectrum (energy, node, l, j, alpha)
    ofstream nout("neutron_spectrum.txt");

    // initiate the while loop
    int nstates_l = 1;
    int l=0;

    // find neutron states for different l
    while(nstates_l>0) {
        nstates_l = 0;

        // check the j = l+1/2 states for bound states
        j = l*1.0+0.5;
        alpha = l*1.0+0.5+0.5; // j=l+1/2 state
        en_n_mev = mN_mev-1e-8-h*20.0; // start at a low energy

        // find where the energy determinant is zero and start at the zeroth node
        node = 0;
        en_n_determinant_prev = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        while(en_n_mev<mN_mev) {
            en_n_determinant = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            if (en_n_determinant_prev*en_n_determinant < 0) {
                bound_state_mev = energy_bisect_n(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev-h,en_n_mev,Fn_unitless,Gn_unitless,fw,fp,ncols_meson);
                nout << fixed << setprecision(15) << bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                node = node+1;
            }
            en_n_determinant_prev = en_n_determinant;
            en_n_mev = en_n_mev + h;
        }
        nstates_l = node;   
        // check the j = l-1/2 states for bound states for l>0
        if (l!=0) {
            j = l*1.0-0.5;
            alpha = -l*1.0-0.5+0.5; // j=l-1/2 state
            en_n_mev = mN_mev-1e-8-h*20.0;;

            node = 0;
            en_n_determinant_prev = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            while(en_n_mev<mN_mev) {
                en_n_determinant = neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
                if (en_n_determinant_prev*en_n_determinant < 0) {
                    bound_state_mev = energy_bisect_n(alpha,meson_fields_unitless,nrows_meson,A,en_n_mev-h,en_n_mev,Fn_unitless,Gn_unitless,fw,fp,ncols_meson);
                    nout << fixed << setprecision(15)<< bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                    node = node+1;
                }
                en_n_determinant_prev = en_n_determinant;
                en_n_mev = en_n_mev + h;
            }
            nstates_l = nstates_l + node;
        }

        l = l + 1;
    }

    dm2.cleanup(Fn_unitless,nrows_meson);
    dm2.cleanup(Gn_unitless,nrows_meson);
}

void energy_spectrum_proton(double** meson_fields_unitless,int nrows_meson, int A, double fw, double fp, int ncols_meson) {
    double en_p_mev;
    double en_p_determinant, en_p_determinant_prev, bound_state_mev;
    double alpha, j;
    double h = 4.0; // energy step
    int node = 0;   // set nodes to zero

    double r_init_fm = meson_fields_unitless[0][0]*r0_fm;
    double r_final_fm = meson_fields_unitless[nrows_meson-1][0]*r0_fm;

    // create the nucleon arrays
    double** Ap_unitless; double** Bp_unitless;
    dm2.create(Ap_unitless,nrows_meson,2);
    dm2.create(Bp_unitless,nrows_meson,2);

    // print out the proton energy spectrum (energy, node, l, j, alpha)
    ofstream pout("proton_spectrum.txt");

    // initiate the while loop
    int nstates_l = 1;
    int l=0;

    // find proton states for different l
    while(nstates_l>0) {
        nstates_l = 0;

        // check the j = l+1/2 states for bound states
        j = l*1.0+0.5;
        alpha = l*1.0+0.5+0.5; // j=l+1/2 state
        en_p_mev = mP_mev-1e-8-h*20.0;

        node = 0;
        en_p_determinant_prev = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        while(en_p_mev<mP_mev) {
            en_p_determinant = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            if (en_p_determinant_prev*en_p_determinant < 0) {
                bound_state_mev = energy_bisect_p(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev-h,en_p_mev,Ap_unitless,Bp_unitless,fw,fp,ncols_meson);
                pout << fixed << setprecision(15) << bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                node = node+1;
            }
            en_p_determinant_prev = en_p_determinant;
            en_p_mev = en_p_mev + h;
        }
        nstates_l = node;

        if (l!=0) {
            // check the j = l-1/2 states for bound states
            j = l*1.0-0.5;
            alpha = -l*1.0-0.5+0.5; // j=l-1/2 state
            en_p_mev = mP_mev-1e-8-h*20.0;

            node = 0;
            en_p_determinant_prev = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
            while(en_p_mev<mP_mev) {
                en_p_determinant = protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm,fw,fp,ncols_meson);
                if (en_p_determinant_prev*en_p_determinant < 0) {
                    bound_state_mev = energy_bisect_p(alpha,meson_fields_unitless,nrows_meson,A,en_p_mev-h,en_p_mev,Ap_unitless,Bp_unitless,fw,fp,ncols_meson);
                    pout << fixed << setprecision(15) << bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                    node = node+1;
                }
                en_p_determinant_prev = en_p_determinant;
                en_p_mev = en_p_mev + h;
            }
            nstates_l = nstates_l + node;
        }
        l = l + 1;
    }

    dm2.cleanup(Ap_unitless,nrows_meson);
    dm2.cleanup(Bp_unitless,nrows_meson);
}

string l_map(int l) {
    string l_array[8] = {"s","p","d","f","g","h","i","j"};
    return l_array[l];
}

string j_map(double j) {
    return to_string(int(2*j)) + ";2";
}

// returns exit code -1 if no states are found, 0 if success
int shell_fill(int A, int Z, int en_col, int j_col) {
    int N = A-Z;
    double quantum_j, fill_frac;
    int nstates;

    int nrows_N = dm2.rowcount("neutron_spectrum.txt");
    int ncols_N = dm2.colcount("neutron_spectrum.txt");
    double** neutron_array;
    if (nrows_N == 0) {
        cout << "No states exist." << endl;
        return -1;
    }

    dm2.importdata("neutron_spectrum.txt",neutron_array);
    dm2.sortasc(neutron_array,0,nrows_N,ncols_N);

    ofstream nout("neutron_spectrum.txt");
    int i=0;
    while (N>0) {
        if (i>(nrows_N-1)) {
            cout << "Not enough states to fill all shells: neutron" << endl;
            cout << "States filled: " << A-Z-N << "/" << A-Z << endl;
            dm2.cleanup(neutron_array,nrows_N);
            return -1;
        }
        quantum_j = neutron_array[i][j_col];
        nstates = 2*quantum_j+1;

        if (nstates <= N) {
            fill_frac = 1.0;
            N = N - nstates;
        } else {
            fill_frac = 1.0*N/nstates;
            N = 0;
        }
        
        for (int k=0; k<ncols_N; ++k) {
            nout << fixed << setprecision(15) << neutron_array[i][k] << "  ";
        }
        nout << fixed << setprecision(15) << fill_frac << endl;

        i = i+1;
    }
    
    dm2.cleanup(neutron_array,nrows_N);

    int nrows_P = dm2.rowcount("proton_spectrum.txt");
    int ncols_P = dm2.colcount("proton_spectrum.txt");
    double** proton_array;

    if (nrows_P == 0) {
        cout << "No states exist." << endl;
        return -1;
    }
    dm2.importdata("proton_spectrum.txt",proton_array);
    dm2.sortasc(proton_array,0,nrows_P,ncols_P);

    ofstream pout("proton_spectrum.txt");
    i=0;
    while (Z>0) {
        if (i>(nrows_P-1)) {
            cout << "Not enough states to fill all shells: proton" << endl;
            cout << "States left: " << Z << endl;
            dm2.cleanup(proton_array,nrows_P);
            return -1;
        }
        quantum_j =proton_array[i][j_col];
        nstates = 2*quantum_j+1;

        if (nstates <= Z) {
            fill_frac = 1.0;
            Z = Z - nstates;
        } else {
            fill_frac = 1.0*Z/nstates;
            Z = 0;
        }
        
        for (int k=0; k<ncols_P; ++k) {
            pout << fixed << setprecision(15) << proton_array[i][k] << "  ";
        }
        pout << fixed << setprecision(15) << fill_frac << endl;

        i = i+1;
    }
    
    dm2.cleanup(proton_array,nrows_P);
    return 0;
}

void get_densities(double** meson_fields_unitless, int A, string energy_spectrum_neutron, string energy_spectrum_proton, double** &densities_svtnp_unitless, int nrows_meson, int ncols_dens, double fw, double fp, int ncols_meson) {
    // set up density matrix
    dm2.create(densities_svtnp_unitless,nrows_meson,ncols_dens);
    dm2.zero(densities_svtnp_unitless,nrows_meson,ncols_dens);

    double r_init_fm = meson_fields_unitless[0][0]*r0_fm;
    double r_final_fm = meson_fields_unitless[nrows_meson-1][0]*r0_fm;
    
    // initialize the quantum numbers, energies, and fill fractions
    double quantum_jn, quantum_jp;
    double nfrac, pfrac, en_n, en_p, alpha;
    
    // fill energy matrix to store energies and quantum numbers
    double** energy_array_neutron;
    double** energy_array_proton;
    int n_levels = dm2.rowcount(energy_spectrum_neutron);
    int p_levels = dm2.rowcount(energy_spectrum_proton);
    dm2.importdata(energy_spectrum_neutron,energy_array_neutron);
    dm2.importdata(energy_spectrum_proton,energy_array_proton);
    
    // initialize wave function arrays
    double Fn_r_unitless, Gn_r_unitless, Ap_r_unitless, Bp_r_unitless, r_unitless;
    double*** Fn_unitless_all; double*** Gn_unitless_all; double*** Ap_unitless_all; double*** Bp_unitless_all;
    dm2.create3d(Fn_unitless_all,nrows_meson,2,n_levels);
    dm2.create3d(Gn_unitless_all,nrows_meson,2,n_levels);
    dm2.create3d(Ap_unitless_all,nrows_meson,2,p_levels);
    dm2.create3d(Bp_unitless_all,nrows_meson,2,p_levels);
    
    // fill the wave function arrays for neutron
    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(dynamic) private(alpha,en_n)
    for (int i=0; i<n_levels; ++i) {
        double** Fn_unitless; double** Gn_unitless;
        en_n = energy_array_neutron[i][0];
        alpha = energy_array_neutron[i][4];
        dm2.create(Fn_unitless,nrows_meson,2); 
        dm2.create(Gn_unitless,nrows_meson,2);
        neutronfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_n,Fn_unitless,Gn_unitless,true,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        normalize(Fn_unitless,Gn_unitless,nrows_meson);
        
        for (int j=0; j<nrows_meson; ++j) {
            Fn_unitless_all[j][0][i] = Fn_unitless[j][0];
            Fn_unitless_all[j][1][i] = Fn_unitless[j][1];
            Gn_unitless_all[j][0][i] = Gn_unitless[j][0];
            Gn_unitless_all[j][1][i] = Gn_unitless[j][1];
        }
        dm2.cleanup(Fn_unitless,nrows_meson);
        dm2.cleanup(Gn_unitless,nrows_meson);
    }
    
    // fill the wave function arrays for proton
    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(dynamic) private(alpha,en_p)
    for (int i=0; i<p_levels; ++i) {
        double** Ap_unitless; double** Bp_unitless;
        en_p = energy_array_proton[i][0];
        alpha = energy_array_proton[i][4];
        dm2.create(Ap_unitless,nrows_meson,2); 
        dm2.create(Bp_unitless,nrows_meson,2);
        protonfield_Solve(alpha,meson_fields_unitless,nrows_meson,A,en_p,Ap_unitless,Bp_unitless,true,r_init_fm,r_final_fm,fw,fp,ncols_meson);
        normalize(Ap_unitless,Bp_unitless,nrows_meson);
        for (int j=0; j<nrows_meson; ++j) {
            Ap_unitless_all[j][0][i] = Ap_unitless[j][0];
            Ap_unitless_all[j][1][i] = Ap_unitless[j][1];
            Bp_unitless_all[j][0][i] = Bp_unitless[j][0];
            Bp_unitless_all[j][1][i] = Bp_unitless[j][1];
        }
        dm2.cleanup(Ap_unitless,nrows_meson);
        dm2.cleanup(Bp_unitless,nrows_meson);
    }
    
    
    // compute the densities
    for (int i=0; i<nrows_meson; ++i) {
        r_unitless = Fn_unitless_all[i][0][0];
        densities_svtnp_unitless[i][0] = r_unitless;
        for (int k=0; k<n_levels; ++k) {
            Fn_r_unitless = Fn_unitless_all[i][1][k];
            Gn_r_unitless = Gn_unitless_all[i][1][k];
            quantum_jn = energy_array_neutron[k][3];
            nfrac = energy_array_neutron[k][5];
            densities_svtnp_unitless[i][1] = densities_svtnp_unitless[i][1] + neutron_scalardens(quantum_jn,nfrac,Fn_r_unitless, Gn_r_unitless,r_unitless);
            densities_svtnp_unitless[i][3] = densities_svtnp_unitless[i][3] + neutron_vectordens(quantum_jn,nfrac,Fn_r_unitless,Gn_r_unitless,r_unitless);
            densities_svtnp_unitless[i][5] = densities_svtnp_unitless[i][5] + neutron_tensordens(quantum_jn,nfrac,Fn_r_unitless,Gn_r_unitless,r_unitless);
        }
        for (int k=0; k<p_levels; ++k) {
            Ap_r_unitless = Ap_unitless_all[i][1][k];
            Bp_r_unitless = Bp_unitless_all[i][1][k];
            quantum_jp = energy_array_proton[k][3];
            pfrac = energy_array_proton[k][5];
            densities_svtnp_unitless[i][2] = densities_svtnp_unitless[i][2] + proton_scalardens(quantum_jp,pfrac,Ap_r_unitless, Bp_r_unitless,r_unitless);
            densities_svtnp_unitless[i][4] = densities_svtnp_unitless[i][4] + proton_vectordens(quantum_jp,pfrac,Ap_r_unitless,Bp_r_unitless,r_unitless);
            densities_svtnp_unitless[i][6] = densities_svtnp_unitless[i][6] + proton_tensordens(quantum_jp,pfrac,Ap_r_unitless,Bp_r_unitless,r_unitless);
        }
    }
    
    // get derivatives of the tensor densities
    divergence_der(nrows_meson,densities_svtnp_unitless,5,9);
    divergence_der(nrows_meson,densities_svtnp_unitless,6,10);
    
    //dm2.print(densities_svtnp_unitless,nrows_meson,ncols_dens,true,"dtensor.txt");
    
    
    //Used for analyzing wave functions
    ofstream out("Fn.txt");
    ofstream lout("Gn.txt");
    ofstream pout("Ap.txt");
    ofstream mout("Bp.txt");
    for (int i=0; i<nrows_meson; ++i) {
        out << scientific << setprecision(10) << Fn_unitless_all[i][0][0];
        lout << scientific << setprecision(10) << Gn_unitless_all[i][0][0];
        pout << scientific << setprecision(10)<< Ap_unitless_all[i][0][0];
        mout << scientific << setprecision(10) << Bp_unitless_all[i][0][0];
        for (int j=0; j<n_levels; ++j) {
            out << "  " << Fn_unitless_all[i][1][j];
            lout << "  " << Gn_unitless_all[i][1][j];
        }
        for (int j=0; j<p_levels; ++j) {
            pout << "  " << Ap_unitless_all[i][1][j];
            mout << "  " << Bp_unitless_all[i][1][j];
        }
        out << endl;
        lout << endl;
        pout << endl;
        mout << endl;
    }
    
    // cleanup
    dm2.cleanup(energy_array_neutron,n_levels);
    dm2.cleanup(energy_array_proton,p_levels);

    dm2.cleanup3d(Fn_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Gn_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Ap_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Bp_unitless_all,nrows_meson,2);
}

double greens_meson(double r_unitless, double rp_unitless, double meson_mass_unitless) {
    double res = 0;

    if (r_unitless>rp_unitless) {
        res = 1.0/meson_mass_unitless*rp_unitless/r_unitless*exp(-meson_mass_unitless*r_unitless*conv_r0_en)*sinh(meson_mass_unitless*rp_unitless*conv_r0_en);
    } else {
        res = 1.0/meson_mass_unitless*rp_unitless/r_unitless*exp(-meson_mass_unitless*rp_unitless*conv_r0_en)*sinh(meson_mass_unitless*r_unitless*conv_r0_en);
    }
    return res;
}

// make sure correct
double greens_coulomb(double r_unitless, double rp_unitless) {
    double res;

    if (r_unitless>rp_unitless) {
        res = 1.0/r_unitless*pow(rp_unitless,2.0);
    } else {
        res = rp_unitless;
    }
    return res;
}

double simpsons(double a, double fa, double b, double fb, double** spline, int low_bound, int up_bound) {
    double m = (a+b)/2.0;
    double fm = dm2.splinecalc(spline,low_bound,up_bound,m);
    double res = (b-a)/6.0*(fa + 4.0*fm + fb);
    return res;
}

#define MAX_DEPTH 20  // Maximum depth limit
double adaptive_simpsons(double a, double fa, double b, double fb, double eps, double S, double m, double fm, double** spline, int low_bound, int up_bound, int depth){
    if (depth > MAX_DEPTH) {
        // Reached maximum depth limit, return an error value or handle the case appropriately
        exit(0);  // Return an error value or handle the case as needed
    }

    double lm = (a+m)/2.0;
    double flm = dm2.splinecalc(spline,low_bound,up_bound,lm);
    double S_left = simpsons(a,fa,m,fm,spline,low_bound,up_bound);

    double rm = (m+b)/2.0;
    double frm = dm2.splinecalc(spline,low_bound,up_bound,rm);
    double S_right = simpsons(m,fm,b,fb,spline,low_bound,up_bound);

    double est_err = 0.0667*fabs(S - S_left - S_right);
    if (est_err < eps) {
        return S;
    } else {
        return adaptive_simpsons(a,fa,m,fm,eps,S_left,lm,flm,spline,low_bound,up_bound,depth+1) + adaptive_simpsons(m,fm,b,fb,eps,S_right,rm,frm,spline,low_bound,up_bound,depth+1);
    }
}

// gives integral from a to b 
double rk4_delta(int index_a, int index_b, double** densities, double mDelta_unitless, int nrows, int ncols, double h_rk4, int integrand, double gd2, double lambda_s, double** old_field_mesons) {
    double res = 0.0;
    double k1,k2,k3;
    double sdensn, sdensp;
    double rp_unitless, rp_unitless_n, r_mid, dens_eff_n, dens_eff_mid;
    double dens_eff;
    double oldfield_sigma_unitless, oldfield_delta_unitless;

    if (integrand == 1) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            sdensn = densities[j][1];
            sdensp = densities[j][2];
            oldfield_sigma_unitless = old_field_mesons[j][1];
            oldfield_delta_unitless = old_field_mesons[j][4];
            dens_eff = 0.5*(sdensp-sdensn) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_delta_unitless*pow(oldfield_sigma_unitless,2.0);

            rp_unitless_n = densities[j+1][0];
            sdensn = densities[j+1][1];
            sdensp = densities[j+1][2];
            oldfield_sigma_unitless = old_field_mesons[j+1][1];
            oldfield_delta_unitless = old_field_mesons[j+1][4];
            dens_eff_n = 0.5*(sdensp-sdensn) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_delta_unitless*pow(oldfield_sigma_unitless,2.0);

            r_mid = (rp_unitless + rp_unitless_n)/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;
            
            k1 = rp_unitless*dens_eff*exp(mDelta_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(mDelta_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(mDelta_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }

    } else if (integrand == 2) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            sdensn = densities[j][1];
            sdensp = densities[j][2];
            oldfield_sigma_unitless = old_field_mesons[j][1];
            oldfield_delta_unitless = old_field_mesons[j][4];
            dens_eff = 0.5*(sdensp-sdensn) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_delta_unitless*pow(oldfield_sigma_unitless,2.0);

            rp_unitless_n = densities[j+1][0];
            sdensn = densities[j+1][1];
            sdensp = densities[j+1][2];
            oldfield_sigma_unitless = old_field_mesons[j+1][1];
            oldfield_delta_unitless = old_field_mesons[j+1][4];
            dens_eff_n = 0.5*(sdensp-sdensn) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_delta_unitless*pow(oldfield_sigma_unitless,2.0);

            r_mid = (rp_unitless + rp_unitless_n)/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;
            
            k1 = rp_unitless*dens_eff*exp(-mDelta_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(-mDelta_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(-mDelta_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }
    } else {
        cout << "invalid integral selection" << endl;
        exit(0);
    }

    return res;
}

// gives integral from a to b 
double rk4_sigma(int index_a, int index_b, double** densities, double mSigma_unitless, int nrows, int ncols, double h_rk4, int integrand, double gs2, double kappa, double lambda, double lambda_s, double** old_field_mesons) {
    double res = 0.0;
    double k1,k2,k3;
    double sdensn, sdensp, oldfield_sigma_unitless, oldfield_delta_unitless;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;
    double rp_unitless, rp_unitless_n, r_mid, dens_eff_n, dens_eff_mid;
    double dens_eff;
    double kappa_unitless = kappa/enscale_mev;

    if (integrand == 1) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            sdensn = densities[j][1];
            sdensp = densities[j][2];
            oldfield_sigma_unitless = old_field_mesons[j][1];
            oldfield_delta_unitless = old_field_mesons[j][4];
            dens_eff = sdensp+sdensn - pow(conv_r0_en,3.0)*0.5*kappa_unitless*pow(oldfield_sigma_unitless,2.0) - pow(conv_r0_en,3.0)*1.0/6.0*lambda*pow(oldfield_sigma_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_sigma_unitless*pow(oldfield_delta_unitless,2.0);
            
            rp_unitless_n = densities[j+1][0];
            sdensn = densities[j+1][1];
            sdensp = densities[j+1][2];
            oldfield_sigma_unitless = old_field_mesons[j+1][1];
            oldfield_delta_unitless = old_field_mesons[j+1][4];
            dens_eff_n = sdensp+sdensn - pow(conv_r0_en,3.0)*0.5*kappa_unitless*pow(oldfield_sigma_unitless,2.0) - pow(conv_r0_en,3.0)*1.0/6.0*lambda*pow(oldfield_sigma_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_sigma_unitless*pow(oldfield_delta_unitless,2.0);

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;

            k1 = rp_unitless*dens_eff*exp(mSigma_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(mSigma_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(mSigma_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }

    } else if (integrand == 2) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            sdensn = densities[j][1];
            sdensp = densities[j][2];
            oldfield_sigma_unitless = old_field_mesons[j][1];
            oldfield_delta_unitless = old_field_mesons[j][4];
            dens_eff = sdensp+sdensn - pow(conv_r0_en,3.0)*0.5*kappa_unitless*pow(oldfield_sigma_unitless,2.0) - pow(conv_r0_en,3.0)*1.0/6.0*lambda*pow(oldfield_sigma_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_sigma_unitless*pow(oldfield_delta_unitless,2.0);

            rp_unitless_n = densities[j+1][0];
            sdensn = densities[j+1][1];
            sdensp = densities[j+1][2];
            oldfield_sigma_unitless = old_field_mesons[j+1][1];
            oldfield_delta_unitless = old_field_mesons[j+1][4];
            dens_eff_n = sdensp+sdensn - pow(conv_r0_en,3.0)*0.5*kappa_unitless*pow(oldfield_sigma_unitless,2.0) - pow(conv_r0_en,3.0)*1.0/6.0*lambda*pow(oldfield_sigma_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_s*oldfield_sigma_unitless*pow(oldfield_delta_unitless,2.0);

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;

            k1 = rp_unitless*dens_eff*exp(-mSigma_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(-mSigma_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(-mSigma_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }
    } else {
        cout << "invalid integral selection" << endl;
        exit(0);
    }

    return res;
}

// gives integral from a to b 
double rk4_omega(int index_a, int index_b, double** densities, double mOmega_unitless, int nrows, int ncols, double h_rk4, int integrand, double gw2, double zeta, double lambda_v, double fw, double** old_field_mesons) {
    double res = 0.0;
    double k1,k2,k3;
    double vdensn, vdensp, oldfield_omega_unitless, oldfield_rho_unitless, div_tdensn, div_tdensp;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;
    double rp_unitless, rp_unitless_n, dens_eff_n, r_mid, dens_eff_mid;
    double dens_eff;

    if (integrand == 1) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            vdensn = densities[j][3];
            vdensp = densities[j][4];
            div_tdensn = densities[j][9];
            div_tdensp = densities[j][10];
            oldfield_omega_unitless = old_field_mesons[j][2];
            oldfield_rho_unitless = old_field_mesons[j][3];
            dens_eff = vdensp+vdensn - pow(conv_r0_en,3.0)*zeta/6.0*pow(oldfield_omega_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless
                        - 0.5/mNuc_unitless*fw*(div_tdensn+div_tdensp)/conv_r0_en;

            rp_unitless_n = densities[j+1][0];
            vdensn = densities[j+1][3];
            vdensp = densities[j+1][4];
            div_tdensn = densities[j+1][9];
            div_tdensp = densities[j+1][10];
            oldfield_omega_unitless = old_field_mesons[j+1][2];
            oldfield_rho_unitless = old_field_mesons[j+1][3];
            dens_eff_n = vdensp+vdensn - pow(conv_r0_en,3.0)*zeta/6.0*pow(oldfield_omega_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless
                        - 0.5/mNuc_unitless*fw*(div_tdensn+div_tdensp)/conv_r0_en;

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;

            k1 = rp_unitless*dens_eff*exp(mOmega_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(mOmega_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(mOmega_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
            //cout << index_a << " : " << j << " : " << index_b << "  " << dens_eff << "  " << dens_eff_n << endl;
        }

    } else if (integrand == 2) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            vdensn = densities[j][3];
            vdensp = densities[j][4];
            div_tdensn = densities[j][9];
            div_tdensp = densities[j][10];
            oldfield_omega_unitless = old_field_mesons[j][2];
            oldfield_rho_unitless = old_field_mesons[j][3];
            dens_eff = vdensp+vdensn - pow(conv_r0_en,3.0)*zeta/6.0*pow(oldfield_omega_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless
                        - 0.5/mNuc_unitless*fw*(div_tdensn+div_tdensp)/conv_r0_en;

            rp_unitless_n = densities[j+1][0];
            vdensn = densities[j+1][3];
            vdensp = densities[j+1][4];
            div_tdensn = densities[j+1][9];
            div_tdensp = densities[j+1][10];
            oldfield_omega_unitless = old_field_mesons[j+1][2];
            oldfield_rho_unitless = old_field_mesons[j+1][3];
            dens_eff_n = vdensp+vdensn - pow(conv_r0_en,3.0)*zeta/6.0*pow(oldfield_omega_unitless,3.0) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless
                        - 0.5/mNuc_unitless*fw*(div_tdensn+div_tdensp)/conv_r0_en;

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;

            k1 = rp_unitless*dens_eff*exp(-mOmega_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(-mOmega_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(-mOmega_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }
    } else {
        cout << "invalid integral selection" << endl;
        exit(0);
    }

    return res;
}

// gives integral from a to b 
double rk4_rho(int index_a, int index_b, double** densities, double mRho_unitless, int nrows, int ncols, double h_rk4, int integrand, double gp2, double xi, double lambda_v, double fp, double** old_field_mesons) {
    double res = 0.0;
    double k1,k2,k3;
    double vdensn, vdensp, oldfield_omega_unitless, oldfield_rho_unitless, div_tdensn, div_tdensp;
    double rp_unitless, rp_unitless_n, dens_eff_n, r_mid, dens_eff_mid;
    double dens_eff;

    if (integrand == 1) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            vdensn = densities[j][3];
            vdensp = densities[j][4];
            div_tdensn = densities[j][9];
            div_tdensp = densities[j][10];
            oldfield_omega_unitless = old_field_mesons[j][2];
            oldfield_rho_unitless = old_field_mesons[j][3];
            dens_eff = 0.5*(vdensp-vdensn) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless - pow(conv_r0_en,3.0)*xi/6.0*pow(oldfield_rho_unitless,3.0)
            -0.25/mNuc_unitless*fp*(div_tdensp-div_tdensn)/conv_r0_en;

            rp_unitless_n = densities[j+1][0];
            vdensn = densities[j+1][3];
            vdensp = densities[j+1][4];
            div_tdensn = densities[j+1][9];
            div_tdensp = densities[j+1][10];
            oldfield_omega_unitless = old_field_mesons[j+1][2];
            oldfield_rho_unitless = old_field_mesons[j+1][3];
            dens_eff_n = 0.5*(vdensp-vdensn) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless - pow(conv_r0_en,3.0)*xi/6.0*pow(oldfield_rho_unitless,3.0)
            -0.25/mNuc_unitless*fp*(div_tdensp-div_tdensn)/conv_r0_en;

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;
            
            k1 = rp_unitless*dens_eff*exp(mRho_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(mRho_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(mRho_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }

    } else if (integrand == 2) {
        for (int j=index_a; j<index_b; ++j) {
            rp_unitless = densities[j][0];
            vdensn = densities[j][3];
            vdensp = densities[j][4];
            div_tdensn = densities[j][9];
            div_tdensp = densities[j][10];
            oldfield_omega_unitless = old_field_mesons[j][2];
            oldfield_rho_unitless = old_field_mesons[j][3];
            dens_eff = 0.5*(vdensp-vdensn) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless - pow(conv_r0_en,3.0)*xi/6.0*pow(oldfield_rho_unitless,3.0)
            -0.25/mNuc_unitless*fp*(div_tdensp-div_tdensn)/conv_r0_en;

            rp_unitless_n = densities[j+1][0];
            vdensn = densities[j+1][3];
            vdensp = densities[j+1][4];
            div_tdensn = densities[j+1][9];
            div_tdensp = densities[j+1][10];
            oldfield_omega_unitless = old_field_mesons[j+1][2];
            oldfield_rho_unitless = old_field_mesons[j+1][3];
            dens_eff_n = 0.5*(vdensp-vdensn) - pow(conv_r0_en,3.0)*2.0*lambda_v*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless - pow(conv_r0_en,3.0)*xi/6.0*pow(oldfield_rho_unitless,3.0)
            -0.25/mNuc_unitless*fp*(div_tdensp-div_tdensn)/conv_r0_en;

            r_mid = rp_unitless + h_rk4/2.0;
            dens_eff_mid = (dens_eff_n-dens_eff)/(rp_unitless_n-rp_unitless)*(r_mid-rp_unitless) + dens_eff;

            k1 = rp_unitless*dens_eff*exp(-mRho_unitless*rp_unitless*conv_r0_en);
            k2 = (rp_unitless+h_rk4/2.0)*dens_eff_mid*exp(-mRho_unitless*(rp_unitless+h_rk4/2.0)*conv_r0_en);
            k3 = (rp_unitless+h_rk4)*dens_eff_n*exp(-mRho_unitless*(rp_unitless+h_rk4)*conv_r0_en);
            res = res + h_rk4/6.0*(k1 + 4.0*k2 + k3);
        }
    } else {
        cout << "invalid integral selection" << endl;
        exit(0);
    }

    return res;
}

double rk4_coulomb(double r_unitless, double** densities, int nrows, double h_rk4) {
    double rp_unitless;
    double vdensp_unitless;
    double e2_charge = 1.4399627/(enscale_mev*r0_fm);
    double res = 0.0;

    // set up integral
    double** integrand;
    dm2.create(integrand,nrows,2);
    for (int i=0; i<nrows; ++i) {
        rp_unitless = densities[i][0];
        integrand[i][0] = densities[i][0];
        vdensp_unitless = densities[i][4];

        integrand[i][1] = 4.0*pi*e2_charge*(vdensp_unitless)*greens_coulomb(r_unitless,rp_unitless);
    }

    double** spline;
    dm2.cubicspline(integrand,spline,nrows,0,1);
    double fa, fb, fmid, a, b, mid, S;
    for (int j=0; j<(nrows-1); j=j+1) {
        a = integrand[j][0];
        b = integrand[j+1][0];
        mid = (a+b)/2.0;
        fa = integrand[j][1];
        fb = integrand[j+1][1];
        fmid = dm2.splinecalc(spline,j,j+1,mid);
        S = h_rk4/6.0*(fa + 4.0*fmid + fb);
        res = res + S;
    }

    dm2.cleanup(spline,nrows);
    dm2.cleanup(integrand,nrows);

    return res;
}

void get_nonlinear_meson_fields(double** &meson_fields_unitless, int npoints_meson, int A, double** densities, int nrows_dens, int ncols_dens, int sdens_n_col, int vdens_n_col, int sdens_p_col, int vdens_p_col, double gs2, double gw2, double gp2, double gd2, double kappa, double lambda, double zeta, double xi, double lambda_v, double lambda_s, double fw, double fp, double mSigma_mev, double mOmega_mev, double mRho_mev, double mDelta_mev, int gridsize_meson, int meson_iterations, int ncols_meson) {
    double I1r, I2r, I20;
    double r_unitless;

    double rp_final_unitless = densities[nrows_dens-1][0];
    double rp_init_unitless = densities[0][0];

    int npoints = gridsize_meson;
    double h_rk4 = (rp_final_unitless-rp_init_unitless)/(npoints-1);

    // get coulomb field
    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(static,6) private(r_unitless)
    for (int i=0; i<npoints_meson; ++i) {
        r_unitless = densities[i][0];
        meson_fields_unitless[i][5] = rk4_coulomb(r_unitless,densities,npoints_meson,h_rk4);
    }

    // Nonlinear convergence
    int MAX_ITER = meson_iterations;
    double** old_meson_fields;
    dm2.create(old_meson_fields,npoints_meson,ncols_meson);
    dm2.copy_pointer(meson_fields_unitless,old_meson_fields,npoints_meson,ncols_meson);

    double mOmega_unitless = mOmega_mev/enscale_mev;
    double mRho_unitless = mRho_mev/enscale_mev;
    double mSigma_unitless = mSigma_mev/enscale_mev;
    double mDelta_unitless = mDelta_mev/enscale_mev;
    
    #pragma omp parallel sections private(r_unitless,I1r,I2r,I20)
    {
        #pragma omp section
        {   
            // Convergence for the Sigma/Delta Field
            for (int k=0; k<MAX_ITER; ++k) {
                // Get the Sigma Field
                r_unitless = densities[0][0];
                I1r = 0.0;
                I2r = rk4_sigma(0,npoints-1,densities,mSigma_unitless,nrows_dens,ncols_dens,h_rk4,2,gs2,kappa,lambda,lambda_s,old_meson_fields);
                I20 = rk4_sigma(0,npoints-1,densities,mSigma_unitless,nrows_dens,ncols_dens,h_rk4,2,gs2,kappa,lambda,lambda_s,old_meson_fields);
                meson_fields_unitless[0][0] = r_unitless;
                meson_fields_unitless[0][1] = gs2/(2.0*mSigma_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mSigma_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mSigma_unitless*r_unitless*conv_r0_en)*I2r);
                old_meson_fields[0][1] = meson_fields_unitless[0][1];

                for (int i=1; i<npoints; ++i) {
                    r_unitless = densities[i][0];
                    meson_fields_unitless[i][0] = r_unitless;

                    I1r = I1r + rk4_sigma((i-1),i,densities,mSigma_unitless,nrows_dens,ncols_dens,h_rk4,1,gs2,kappa,lambda,lambda_s,old_meson_fields);
                    I2r = rk4_sigma(i,npoints-1,densities,mSigma_unitless,nrows_dens,ncols_dens,h_rk4,2,gs2,kappa,lambda,lambda_s,old_meson_fields);
                    meson_fields_unitless[i][1] = gs2/(2.0*mSigma_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mSigma_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mSigma_unitless*r_unitless*conv_r0_en)*I2r);
                    old_meson_fields[i][1] = meson_fields_unitless[i][1];
                }

                r_unitless = densities[0][0];
                I1r = 0.0;
                I2r = rk4_delta(0,npoints-1,densities,mDelta_unitless,nrows_dens,ncols_dens,h_rk4,2,gd2,lambda_s,old_meson_fields);
                I20 = rk4_delta(0,npoints-1,densities,mDelta_unitless,nrows_dens,ncols_dens,h_rk4,2,gd2,lambda_s,old_meson_fields);
                meson_fields_unitless[0][0] = r_unitless;
                meson_fields_unitless[0][4] = gd2/(2.0*mDelta_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mDelta_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mDelta_unitless*r_unitless*conv_r0_en)*I2r);
                old_meson_fields[0][4] = meson_fields_unitless[0][4];

                for (int i=1; i<npoints; ++i) {
                    r_unitless = densities[i][0];
                    meson_fields_unitless[i][0] = r_unitless;

                    I1r = I1r + rk4_delta((i-1),i,densities,mDelta_unitless,nrows_dens,ncols_dens,h_rk4,1,gd2,lambda_s,old_meson_fields);
                    I2r = rk4_delta(i,npoints-1,densities,mDelta_unitless,nrows_dens,ncols_dens,h_rk4,2,gd2,lambda_s,old_meson_fields);
                    meson_fields_unitless[i][4] = gd2/(2.0*mDelta_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mDelta_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mDelta_unitless*r_unitless*conv_r0_en)*I2r);
                    old_meson_fields[i][4] = meson_fields_unitless[i][4];
                }
            }
        }
        #pragma omp section
        {
            // Convergence for the Omega/Rho Field (new method)
            for (int k=0; k<MAX_ITER; ++k) {
                // Get the Omega field
                r_unitless = densities[0][0];
                I1r = 0.0;
                I2r = rk4_omega(0,npoints-1,densities,mOmega_unitless,nrows_dens,ncols_dens,h_rk4,2,gw2,zeta,lambda_v,fw,old_meson_fields);
                I20 = rk4_omega(0,npoints-1,densities,mOmega_unitless,nrows_dens,ncols_dens,h_rk4,2,gw2,zeta,lambda_v,fw,old_meson_fields);
                meson_fields_unitless[0][0] = r_unitless;
                meson_fields_unitless[0][2] = gw2/(2.0*mOmega_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mOmega_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mOmega_unitless*r_unitless*conv_r0_en)*I2r);
                old_meson_fields[0][2] = meson_fields_unitless[0][2];

                for (int i=1; i<npoints; ++i) {
                    r_unitless = densities[i][0];
                    meson_fields_unitless[i][0] = r_unitless;

                    I1r = I1r + rk4_omega((i-1),i,densities,mOmega_unitless,nrows_dens,ncols_dens,h_rk4,1,gw2,zeta,lambda_v,fw,old_meson_fields);
                    I2r = rk4_omega(i,npoints-1,densities,mOmega_unitless,nrows_dens,ncols_dens,h_rk4,2,gw2,zeta,lambda_v,fw,old_meson_fields);
                    meson_fields_unitless[i][2] = gw2/(2.0*mOmega_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mOmega_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mOmega_unitless*r_unitless*conv_r0_en)*I2r);
                    old_meson_fields[i][2] = meson_fields_unitless[i][2];
                }

                r_unitless = densities[0][0];
                I1r = 0.0;
                I2r = rk4_rho(0,npoints-1,densities,mRho_unitless,nrows_dens,ncols_dens,h_rk4,2,gp2,xi,lambda_v,fp,old_meson_fields);
                I20 = rk4_rho(0,npoints-1,densities,mRho_unitless,nrows_dens,ncols_dens,h_rk4,2,gp2,xi,lambda_v,fp,old_meson_fields);
                meson_fields_unitless[0][0] = r_unitless;
                meson_fields_unitless[0][3] = gp2/(2.0*mRho_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mRho_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mRho_unitless*r_unitless*conv_r0_en)*I2r);
                meson_fields_unitless[0][3] = (old_meson_fields[0][3] + meson_fields_unitless[0][3])/2.0;
                old_meson_fields[0][3] = meson_fields_unitless[0][3];

                for (int i=1; i<npoints; ++i) {
                    r_unitless = densities[i][0];
                    meson_fields_unitless[i][0] = r_unitless;

                    I1r = I1r + rk4_rho((i-1),i,densities,mRho_unitless,nrows_dens,ncols_dens,h_rk4,1,gp2,xi,lambda_v,fp,old_meson_fields);
                    I2r = rk4_rho(i,npoints-1,densities,mRho_unitless,nrows_dens,ncols_dens,h_rk4,2,gp2,xi,lambda_v,fp,old_meson_fields);
                    meson_fields_unitless[i][3] = gp2/(2.0*mRho_unitless*r_unitless*pow(conv_r0_en,2.0))*(exp(-mRho_unitless*r_unitless*conv_r0_en)*(I1r - I20) + exp(mRho_unitless*r_unitless*conv_r0_en)*I2r);
                    meson_fields_unitless[i][3] = (old_meson_fields[i][3] + meson_fields_unitless[i][3])/2.0;
                    old_meson_fields[i][3] = meson_fields_unitless[i][3];
                }
            }
        }
    }

    // get the derivative for the omega and rho
    meson_der(npoints_meson,meson_fields_unitless,2,6);
    meson_der(npoints_meson,meson_fields_unitless,3,7);

    dm2.cleanup(old_meson_fields,npoints_meson);

}

// return the Binding energy in MeV
double get_BA(double** meson_field_unitless, double** densities_unitless, string n_energies, string p_energies, int npoints_meson, int npoints_densities, int ncols_density, int A, double kappa, double lambda, double zeta, double xi, double lambda_v, double lambda_s, double fw, double fp) {
    double BA_unitless = 0;
    double** n_spectrum;
    double** p_spectrum;
    double en_neutrons = 0;
    double en_protons = 0;

    dm2.importdata(n_energies,n_spectrum);
    dm2.importdata(p_energies,p_spectrum);
    double nrows = dm2.rowcount(n_energies);
    double prows = dm2.rowcount(p_energies);

    for (int i=0; i<nrows; ++i) {
        en_neutrons = en_neutrons + (n_spectrum[i][0])*(2.0*n_spectrum[i][3]+1)*n_spectrum[i][5]/enscale_mev;
        //cout << "neutron spectrum: " << n_spectrum[i][0]/enscale_mev << endl;
    }

    for (int i=0; i<prows; ++i) {
        en_protons = en_protons + (p_spectrum[i][0])*(2.0*p_spectrum[i][3]+1)*p_spectrum[i][5]/enscale_mev;
        //cout << "proton spectrum: " << p_spectrum[i][0]/enscale_mev << endl;
    }
    
    double r_init_unitless = densities_unitless[0][0];
    double r_final_unitless = densities_unitless[npoints_meson-1][0];
    double h_rk4 = (r_final_unitless-r_init_unitless)/(npoints_meson-1);
    double r_unitless;
    double integral = 0.0;
    double en_mesons_integrand, en_coulomb_integrand;
    double** integrand;
    dm2.create(integrand,npoints_meson,2);
    double gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, gd_delta_r_unitless, e_coulomb_r_unitless, ns_density_r_unitless, ps_density_r_unitless, nv_density_r_unitless, pv_density_r_unitless;
    double div_nt_density_r_unitless, div_pt_density_r_unitless; 
    double kappa_unitless = kappa/enscale_mev;
    //ofstream out("dtensor.txt");
    for (int i=0; i<npoints_meson; ++i) {
        r_unitless = meson_field_unitless[i][0];
        gs_sigma_r_unitless = meson_field_unitless[i][1];
        gw_omega_r_unitless = meson_field_unitless[i][2];
        gp_rho_r_unitless = meson_field_unitless[i][3];
        gd_delta_r_unitless = meson_field_unitless[i][4];
        e_coulomb_r_unitless = meson_field_unitless[i][5];
        ns_density_r_unitless = densities_unitless[i][1];
        ps_density_r_unitless = densities_unitless[i][2];
        nv_density_r_unitless = densities_unitless[i][3];
        pv_density_r_unitless = densities_unitless[i][4];
        div_nt_density_r_unitless = densities_unitless[i][9];
        div_pt_density_r_unitless = densities_unitless[i][10];

        
        en_mesons_integrand = gs_sigma_r_unitless*(ns_density_r_unitless+ps_density_r_unitless) + 0.5*gd_delta_r_unitless*(ps_density_r_unitless - ns_density_r_unitless) - gw_omega_r_unitless*(nv_density_r_unitless+pv_density_r_unitless)
                            - 0.5*gp_rho_r_unitless*(pv_density_r_unitless-nv_density_r_unitless) - 1.0/6.0*kappa_unitless*pow(conv_r0_en,3.0)*pow(gs_sigma_r_unitless,3.0) - 1.0/12.0*lambda*pow(conv_r0_en,3.0)*pow(gs_sigma_r_unitless,4.0)
                            - 2.0*lambda_s*pow(conv_r0_en,3.0)*pow(gs_sigma_r_unitless,2.0)*pow(gd_delta_r_unitless,2.0) + 1.0/12.0*zeta*pow(gw_omega_r_unitless,4.0)*pow(conv_r0_en,3.0)
                            + 2.0*lambda_v*pow(gw_omega_r_unitless,2.0)*pow(gp_rho_r_unitless,2.0)*pow(conv_r0_en,3.0) + 1.0/12.0*xi*pow(gp_rho_r_unitless,4.0)*pow(conv_r0_en,3.0)
                            + 0.5/mNuc_unitless*fw*gw_omega_r_unitless*(div_nt_density_r_unitless+div_pt_density_r_unitless)/conv_r0_en
                            + 0.25/mNuc_unitless*fp*gp_rho_r_unitless*(div_pt_density_r_unitless-div_nt_density_r_unitless)/conv_r0_en;
    
        en_coulomb_integrand = - e_coulomb_r_unitless*pv_density_r_unitless;
        integrand[i][0] = r_unitless;
        integrand[i][1] = (en_mesons_integrand+en_coulomb_integrand)*pow(r_unitless,2.0);
        //out << r_unitless*r0_fm << "  " << (div_pt_density_r_unitless+div_nt_density_r_unitless)/pow(r0_fm,4) << endl;
        //cout << gs_sigma_r_unitless << "  " << gw_omega_r_unitless << "  " << gp_rho_r_unitless << "  " << e_coulomb_r_unitless << endl;
    }


    double** spline;
    dm2.cubicspline(integrand,spline,npoints_meson,0,1);
    double fa, fb, fmid, a, b, mid, S;
    for (int j=0; j<(npoints_meson-1); ++j) {
        a = integrand[j][0];
        b = integrand[j+1][0];
        mid = (a+b)/2.0;
        fa = integrand[j][1];
        fb = integrand[j+1][1];
        fmid = dm2.splinecalc(spline,j,j+1,mid);
        S = h_rk4/6.0*(fa + 4.0*fmid + fb);
        integral = integral + S;
    }
    dm2.cleanup(integrand,npoints_meson);
    dm2.cleanup(spline,npoints_meson);

    BA_unitless = en_neutrons + en_protons + 2.0*pi*integral;
    //cout << " " << en_neutrons*enscale_mev/A << "  " << en_protons*enscale_mev/A << "  " << 2.0*pi*integral*enscale_mev/A << endl;

    dm2.cleanup(n_spectrum,nrows);
    dm2.cleanup(p_spectrum,prows);
    return (BA_unitless*enscale_mev - 0.75*41.0*pow(A,-1.0/3.0))/A - mNuc_mev;
    //return (BA_unitless*enscale_mev - 17.2*pow(A,-0.2))/A - mNuc_mev;
}

void get_nucleon_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_N_P_fm2[2]) {
    int nsteps_rk4 = npoints_densities;
    double r_init_unitless = densities_unitless[0][0];
    double r_final_unitless = densities_unitless[npoints_densities-1][0];
    double h = (r_final_unitless-r_init_unitless)/(nsteps_rk4-1);
    double r_unitless;
    double nv_density_r_unitless, pv_density_r_unitless;
    double integralp = 0; double integraln = 0;

    for (int i=0; i<nsteps_rk4; ++i) {
        r_unitless = densities_unitless[i][0];
        nv_density_r_unitless = densities_unitless[i][3];
        pv_density_r_unitless = densities_unitless[i][4];
        integralp = integralp + h*(pow(r_unitless,4.0)*pv_density_r_unitless);
        integraln = integraln + h*(pow(r_unitless,4.0)*nv_density_r_unitless);
    }
    Radii2_N_P_fm2[0] = pow(r0_fm,2.0)*(4.0*pi/(A-Z))*integraln;
    Radii2_N_P_fm2[1] = pow(r0_fm,2.0)*(4.0*pi/Z)*integralp;
}

void subtract_restmass(string proton_spectrum, string neutron_spectrum) {
    double** n_array;
    double** p_array;

    dm2.importdata(neutron_spectrum,n_array);
    dm2.importdata(proton_spectrum,p_array);
    double nrows = dm2.rowcount(neutron_spectrum);
    double prows = dm2.rowcount(proton_spectrum);

    for (int i=0; i<nrows; ++i) {
        //n_array[i][0] = n_array[i][0] - mNuc_mev;
        n_array[i][0] = n_array[i][0]/enscale_mev;
    }

    for (int i=0; i<prows; ++i) {
        //p_array[i][0] = p_array[i][0] - mNuc_mev;
        p_array[i][0] = p_array[i][0]/enscale_mev;
    }

    ofstream nout("neutron_spectrum.txt");
    ofstream pout("proton_spectrum.txt");

    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<6; ++j) {
            nout << n_array[i][j] << "  ";
        }
        nout << to_string(int(n_array[i][1]+1)) + l_map(int(n_array[i][2])) + j_map(n_array[i][3]) << endl;
    }

    for (int i=0; i<prows; ++i) {
        for (int j=0; j<6; ++j) {
            pout << p_array[i][j] << "  ";
        }
        pout << to_string(int(p_array[i][1]+1)) + l_map(int(p_array[i][2])) + j_map(p_array[i][3]) << endl;
    }
    dm2.cleanup(n_array,nrows);
    dm2.cleanup(p_array,prows);
}

void EM_formfactors(double GE_pn[2], double GM_pn[2], double q) {
    double hbar_gevfm = hbar_mevfm*1e-3;
    double Q2 = pow(q*hbar_gevfm,2.0);  // GeV^2
    double mup = 2.79284356;    // proton magnetic moment
    double mun = -1.91304272;    // proton magnetic moment

    double a_GE_p[13] = {0.239163298067,-1.10985857441,1.44438081306,0.479569465603,-2.28689474187,1.12663298498,1.250619843540,
                    -3.63102047159,4.08221702379,0.504097346499,-5.08512046051,3.96774254395,-0.981529071103};
    double a_GM_p[13] = {0.264142994136,-1.09530612212,1.21855378178,0.661136493537,-1.40567892503,-1.35641843888,1.447029155340,
                    4.23566973590,-5.33404565341,-2.916300520960,8.70740306757,-5.70699994375,1.280814375890};
    double a_GE_n[13] = {0.048919981379,-0.064525053912,-0.240825897382,0.392108744873,0.300445258602,-0.661888687179,-0.175639769687,
                    0.624691724461,-0.077684299367,-0.236003975259,0.090401973470,0.000000000000,0.000000000000};
    double a_GM_n[13] = {0.257758326959,-1.079540642058,1.182183812195,0.711015085833,-1.348080936796,-1.662444025208,2.624354426029,
                    1.751234494568,-4.922300878888,3.197892727312,-0.712072389946,0.000000000000,0.000000000000};
    
    double t0 = -0.70000000000;   // GeV^2
    double tcut = 0.0779191396;    // GeV^2
    double z = (sqrt(tcut+Q2)-sqrt(tcut-t0))/(sqrt(tcut+Q2)+sqrt(tcut-t0)); // UNITLESS

    double GE_proton = 0; double GM_proton = 0;
    double GE_neutron = 0; double GM_neutron = 0;

    for (int i=0; i<13; ++i) {
        GE_proton = GE_proton + a_GE_p[i]*pow(z,i*1.0);
        GM_proton = GM_proton + a_GM_p[i]*pow(z,i*1.0);
        GE_neutron = GE_neutron + a_GE_n[i]*pow(z,i*1.0);
        GM_neutron = GM_neutron + a_GM_n[i]*pow(z,i*1.0);
    }
    GM_proton = mup*GM_proton;
    GM_neutron = mun*GM_neutron;

    GE_pn[0] = GE_proton;
    GE_pn[1] = GE_neutron;
    GM_pn[0] = GM_proton;
    GM_pn[1] = GM_neutron;
}

void WEAK_formfactors(double WGE_pn[2], double WGM_pn[2], double q, double qwn) {
    double GE_pn[2];
    double GM_pn[2];
    EM_formfactors(GE_pn, GM_pn, q);
    
    WGE_pn[0] = qwp*GE_pn[0]+qwn*GE_pn[1];  //GE_weak (proton)
    WGE_pn[1] = qwp*GE_pn[1]+qwn*GE_pn[0];  // GE_weak (neutron)
    WGM_pn[0] = qwp*GM_pn[0]+qwn*GM_pn[1];  //GE_weak (proton)
    WGM_pn[1] = qwp*GM_pn[1]+qwn*GM_pn[0];  // GE_weak (neutron)
}

void vector_formfactor(double** densities_svtnp_unitless, double q_unitless, int nrows, double vFF_pn[2]) {

    double r_unitless_a, r_unitless_b, r_unitless_ab2, vectordens_unitless_n_a, vectordens_unitless_n_b, vectordens_unitless_n_ab2, fa, fb, fab2;
    vFF_pn[1] = 0;  // neutron
    for (int j=0; j<(nrows-1); ++j) {
        r_unitless_a = densities_svtnp_unitless[j][0];
        r_unitless_b = densities_svtnp_unitless[j+1][0];
        r_unitless_ab2 = 0.5*(r_unitless_a+r_unitless_b);
        vectordens_unitless_n_a = densities_svtnp_unitless[j][3];
        vectordens_unitless_n_b = densities_svtnp_unitless[j+1][3];
        vectordens_unitless_n_ab2 = dm2.interpolate(nrows,7,densities_svtnp_unitless,r_unitless_ab2,0,3,true);
        fa = 4.0*pi*vectordens_unitless_n_a*sin(q_unitless*r_unitless_a)/(q_unitless*r_unitless_a)*pow(r_unitless_a,2.0);
        fb = 4.0*pi*vectordens_unitless_n_b*sin(q_unitless*r_unitless_b)/(q_unitless*r_unitless_b)*pow(r_unitless_b,2.0);
        fab2 = 4.0*pi*vectordens_unitless_n_ab2*sin(q_unitless*r_unitless_ab2)/(q_unitless*r_unitless_ab2)*pow(r_unitless_ab2,2.0);

        vFF_pn[1] = vFF_pn[1] + (r_unitless_b-r_unitless_a)/6.0*(fa + 4.0*fab2 + fb);
    }

    double vectordens_unitless_p_a, vectordens_unitless_p_b, vectordens_unitless_p_ab2;
    vFF_pn[0] = 0; // proton
    for (int j=0; j<(nrows-1); ++j) {
        r_unitless_a = densities_svtnp_unitless[j][0];
        r_unitless_b = densities_svtnp_unitless[j+1][0];
        r_unitless_ab2 = 0.5*(r_unitless_a+r_unitless_b);
        vectordens_unitless_p_a = densities_svtnp_unitless[j][4];
        vectordens_unitless_p_b = densities_svtnp_unitless[j+1][4];
        fa = 4.0*pi*vectordens_unitless_p_a*sin(q_unitless*r_unitless_a)/(q_unitless*r_unitless_a)*pow(r_unitless_a,2.0);
        fb = 4.0*pi*vectordens_unitless_p_b*sin(q_unitless*r_unitless_b)/(q_unitless*r_unitless_b)*pow(r_unitless_b,2.0);
        vectordens_unitless_p_ab2 = dm2.interpolate(nrows,7,densities_svtnp_unitless,r_unitless_ab2,0,4,true);
        fab2 = 4.0*pi*vectordens_unitless_p_ab2*sin(q_unitless*r_unitless_ab2)/(q_unitless*r_unitless_ab2)*pow(r_unitless_ab2,2.0);

        vFF_pn[0] = vFF_pn[0] + (r_unitless_b-r_unitless_a)/6.0*(fa + 4.0*fab2 + fb);
    }
    //cout << q_unitless << "  " << vFF_pn[0] << "  " << vFF_pn[1] << endl;
}

void tensor_formfactor(double** densities_svtnp_unitless, double q_unitless, int nrows, double tFF_pn[2]) {
    
    double r_unitless_a, r_unitless_b, tensordens_unitless_n_a, tensordens_unitless_n_b, fa, fb, r_unitless_ab2, tensordens_unitless_n_ab2, fab2;
    // neutron
    tFF_pn[1] = 0;
    for (int j=0; j<(nrows-1); ++j) {
        r_unitless_a = densities_svtnp_unitless[j][0];
        r_unitless_b = densities_svtnp_unitless[j+1][0];
        r_unitless_ab2 = 0.5*(r_unitless_a+r_unitless_b);
        tensordens_unitless_n_a = densities_svtnp_unitless[j][5];
        tensordens_unitless_n_b = densities_svtnp_unitless[j+1][5];
        fa = 4.0*pi*tensordens_unitless_n_a*(sin(q_unitless*r_unitless_a)/pow(q_unitless*r_unitless_a,2.0) - cos(q_unitless*r_unitless_a)/(q_unitless*r_unitless_a))*pow(r_unitless_a,2.0);
        fb = 4.0*pi*tensordens_unitless_n_b*(sin(q_unitless*r_unitless_b)/pow(q_unitless*r_unitless_b,2.0) - cos(q_unitless*r_unitless_b)/(q_unitless*r_unitless_b))*pow(r_unitless_b,2.0);
        tensordens_unitless_n_ab2 = dm2.interpolate(nrows,7,densities_svtnp_unitless,r_unitless_ab2,0,5,true);
        fab2 = 4.0*pi*tensordens_unitless_n_ab2*(sin(q_unitless*r_unitless_ab2)/pow(q_unitless*r_unitless_ab2,2.0) - cos(q_unitless*r_unitless_ab2)/(q_unitless*r_unitless_ab2))*pow(r_unitless_ab2,2.0);

        tFF_pn[1] = tFF_pn[1] + (r_unitless_b-r_unitless_a)/6.0*(fa + 4.0*fab2 + fb);
    }

    // proton
    double tensordens_unitless_p_a, tensordens_unitless_p_b, tensordens_unitless_p_ab2;
    tFF_pn[0] = 0;
    for (int j=0; j<(nrows-1); ++j) {
        r_unitless_a = densities_svtnp_unitless[j][0];
        r_unitless_b = densities_svtnp_unitless[j+1][0];
        r_unitless_ab2 = 0.5*(r_unitless_a+r_unitless_b);
        tensordens_unitless_p_a = densities_svtnp_unitless[j][6];
        tensordens_unitless_p_b = densities_svtnp_unitless[j+1][6];
        fa = 4.0*pi*tensordens_unitless_p_a*(sin(q_unitless*r_unitless_a)/pow(q_unitless*r_unitless_a,2.0) - cos(q_unitless*r_unitless_a)/(q_unitless*r_unitless_a))*pow(r_unitless_a,2.0);
        fb = 4.0*pi*tensordens_unitless_p_b*(sin(q_unitless*r_unitless_b)/pow(q_unitless*r_unitless_b,2.0) - cos(q_unitless*r_unitless_b)/(q_unitless*r_unitless_b))*pow(r_unitless_b,2.0);
        tensordens_unitless_p_ab2 = dm2.interpolate(nrows,7,densities_svtnp_unitless,r_unitless_ab2,0,6,true);
        fab2 = 4.0*pi*tensordens_unitless_p_ab2*(sin(q_unitless*r_unitless_ab2)/pow(q_unitless*r_unitless_ab2,2.0) - cos(q_unitless*r_unitless_ab2)/(q_unitless*r_unitless_ab2))*pow(r_unitless_ab2,2.0);

        tFF_pn[0] = tFF_pn[0] + (r_unitless_b-r_unitless_a)/6.0*(fa + 4.0*fab2 + fb);
    }
    
}

double charge_formfactor(double q_unitless, double** densities_svtnp_unitless, int nrows, int Z) {
    double q_ifm = q_unitless/r0_fm;
    double GE_pn[2] = {0, 0};
    double GM_pn[2] = {0, 0};
    double vFF_pn[2] = {0, 0};
    double tFF_pn[2] = {0, 0};
    double cFF_p, cFF_n, cFF;

    double hbar_gevfm = hbar_mevfm*1e-3;
    double mNuc_gev = mNuc_mev*1e-3;
    double Q2 = pow(q_ifm*hbar_gevfm,2.0);  // GeV^2
    double tau = 0.25*Q2/pow(mNuc_gev,2.0);

    EM_formfactors(GE_pn, GM_pn, q_ifm);
    vector_formfactor(densities_svtnp_unitless,q_unitless,nrows,vFF_pn);
    tensor_formfactor(densities_svtnp_unitless,q_unitless,nrows,tFF_pn);
    cFF_p = GE_pn[0]*vFF_pn[0] + (GM_pn[0] - GE_pn[0])/(1.0+tau)*(tau*vFF_pn[0] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[0]);
    cFF_n = GE_pn[1]*vFF_pn[1] + (GM_pn[1] - GE_pn[1])/(1.0+tau)*(tau*vFF_pn[1] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[1]);
    cFF = cFF_p + cFF_n;

    return cFF/Z;
}

double weak_formfactor(double q_unitless, double** densities_svtnp_unitless, int nrows, int A, int Z, double qwn) {
    double q_ifm = q_unitless/r0_fm;
    double WGE_pn[2] = {0, 0};
    double WGM_pn[2] = {0, 0};
    double vFF_pn[2] = {0, 0};
    double tFF_pn[2] = {0, 0};
    double wFF_p, wFF_n, wFF;

    double hbar_gevfm = hbar_mevfm*1e-3;
    double mNuc_gev = mNuc_mev*1e-3;
    double Q2 = pow(q_ifm*hbar_gevfm,2.0);  // GeV^2
    double tau = 0.25*Q2/pow(mNuc_gev,2.0);

    WEAK_formfactors(WGE_pn,WGM_pn,q_ifm,qwn);
    vector_formfactor(densities_svtnp_unitless,q_unitless,nrows,vFF_pn);
    tensor_formfactor(densities_svtnp_unitless,q_unitless,nrows,tFF_pn);
    wFF_p = WGE_pn[0]*vFF_pn[0] + (WGM_pn[0] - WGE_pn[0])/(1.0+tau)*(tau*vFF_pn[0] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[0]);
    wFF_n = WGE_pn[1]*vFF_pn[1] + (WGM_pn[1] - WGE_pn[1])/(1.0+tau)*(tau*vFF_pn[1] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[1]);
    wFF = wFF_p + wFF_n;

    double Qwk = Z*qwp + (A-Z)*qwn;
    return wFF/Qwk;
}

// version using the array (returns Fch-Fwk)
double get_WEAK_CHARGE_densities_v2(double** &densities_svtnp_unitless, int Z, int A, int ncols_dens, int nrows_dens, double qwn){

    double r_unitless, q_unitless, charge_dens_unitless, weak_dens_unitless, q_transfer_exp;
    double Qwk = Z*qwp + (A-Z)*qwn;

    // get charge and weak form factors first
    int npoints = nrows_dens;
    int q_transfer_row = 0;
    double** q_charge_weak_factors;
    dm2.create(q_charge_weak_factors,npoints,4);

    //#pragma omp parallel num_threads(12)
    //#pragma omp for schedule(static,1) private(q_unitless)
    for (int i=0; i<npoints; ++i) {
        q_unitless = 1e-15 + 0.01*i;
        q_charge_weak_factors[i][0] = q_unitless;
        q_charge_weak_factors[i][1] = charge_formfactor(q_unitless,densities_svtnp_unitless,nrows_dens,Z);
        q_charge_weak_factors[i][2] = weak_formfactor(q_unitless,densities_svtnp_unitless,nrows_dens,A,Z,qwn);
        q_charge_weak_factors[i][3] = q_charge_weak_factors[i][1]-q_charge_weak_factors[i][2];
    }

    // print the weak skin
    if (A==48) {
        //dm2.print(q_charge_weak_factors,npoints,4,true,"q_transfer_Fch,Fwk,wkskin.txt");
    }
    // inverse fourier transform
    double q_unitless_a, q_unitless_b, charge_factor_a, charge_factor_b, fa, fb, q_unitless_ab2, charge_factor_ab2, fab2, weak_factor_a, weak_factor_b, weak_factor_ab2, ga, gb, gab2;
    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(q_unitless_a,q_unitless_b,q_unitless_ab2,charge_factor_a,charge_factor_b,charge_factor_ab2,fa,fb,fab2,r_unitless,charge_dens_unitless,weak_factor_a, weak_factor_b, weak_factor_ab2, ga, gb, gab2, weak_dens_unitless)
    for (int i=0; i<nrows_dens; ++i) {
        r_unitless = densities_svtnp_unitless[i][0];
        charge_dens_unitless = 0.0;
        weak_dens_unitless = 0.0;
        for (int j=0; j<(npoints-1); ++j) {
            q_unitless_a = q_charge_weak_factors[j][0];
            q_unitless_b = q_charge_weak_factors[j+1][0];
            q_unitless_ab2 = 0.5*(q_unitless_a+q_unitless_b);
            
            charge_factor_a = q_charge_weak_factors[j][1];
            charge_factor_b = q_charge_weak_factors[j+1][1];
            fa = Z/(2.0*pi*pi)*pow(q_unitless_a,2.0)*charge_factor_a*sin(q_unitless_a*r_unitless)/(q_unitless_a*r_unitless);
            fb = Z/(2.0*pi*pi)*pow(q_unitless_b,2.0)*charge_factor_b*sin(q_unitless_b*r_unitless)/(q_unitless_b*r_unitless);
            charge_factor_ab2 = dm2.interpolate(npoints,3,q_charge_weak_factors,q_unitless_ab2,0,1,true);
            fab2 = Z/(2.0*pi*pi)*pow(q_unitless_ab2,2.0)*charge_factor_ab2*sin(q_unitless_ab2*r_unitless)/(q_unitless_ab2*r_unitless);

            weak_factor_a = q_charge_weak_factors[j][2];
            weak_factor_b = q_charge_weak_factors[j+1][2];
            ga = Qwk/(2.0*pi*pi)*pow(q_unitless_a,2.0)*weak_factor_a*sin(q_unitless_a*r_unitless)/(q_unitless_a*r_unitless);
            gb = Qwk/(2.0*pi*pi)*pow(q_unitless_b,2.0)*weak_factor_b*sin(q_unitless_b*r_unitless)/(q_unitless_b*r_unitless);
            weak_factor_ab2 = dm2.interpolate(npoints,3,q_charge_weak_factors,q_unitless_ab2,0,2,true);
            gab2 = Qwk/(2.0*pi*pi)*pow(q_unitless_ab2,2.0)*weak_factor_ab2*sin(q_unitless_ab2*r_unitless)/(q_unitless_ab2*r_unitless);

            charge_dens_unitless = charge_dens_unitless + (q_unitless_b-q_unitless_a)/6.0*(fa + 4.0*fab2 + fb);
            weak_dens_unitless = weak_dens_unitless + (q_unitless_b-q_unitless_a)/6.0*(ga + 4.0*gab2 + gb);
        }

        densities_svtnp_unitless[i][7] = charge_dens_unitless;
        densities_svtnp_unitless[i][8] = weak_dens_unitless;
    }

    if (Z==20 && A==48) {
        q_transfer_exp = q_transfer_Ca;
        q_transfer_row = dm2.findvalue(q_charge_weak_factors,nrows_dens,4,q_transfer_exp*r0_fm,0,0.1);
        double Fch_Fwk = q_charge_weak_factors[q_transfer_row][3];
        dm2.cleanup(q_charge_weak_factors,npoints);
        //cout << "Fch: " << q_charge_weak_factors[q_transfer_row][1] << endl;
        return Fch_Fwk;
    } else if (Z==82 && A ==208) {
        q_transfer_exp = q_transfer_Pb;
        q_transfer_row = dm2.findvalue(q_charge_weak_factors,nrows_dens,4,q_transfer_exp*r0_fm,0,0.1);
        double Fch_Fwk = q_charge_weak_factors[q_transfer_row][3];
        dm2.cleanup(q_charge_weak_factors,npoints);
        //cout << "Fch: " << q_charge_weak_factors[q_transfer_row][1] << endl;
        return Fch_Fwk;
    } else {
        dm2.cleanup(q_charge_weak_factors,npoints);
        return 0;
    }
}

void get_weak_charge_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_C_W_fm2[2], double qwn) {
    int nsteps_rk4 = npoints_densities;
    double r_init_unitless = densities_unitless[0][0];
    double r_final_unitless = densities_unitless[npoints_densities-1][0];
    double h = (r_final_unitless-r_init_unitless)/(nsteps_rk4-1);
    double r_unitless;
    double charge_density_r_unitless, weak_density_r_unitless;
    double integralcharge = 0; double integralweak = 0;
    double Qwk = Z*qwp + (A-Z)*qwn;

    for (int i=0; i<nsteps_rk4; ++i) {
        r_unitless = densities_unitless[i][0];
        charge_density_r_unitless = densities_unitless[i][7];
        weak_density_r_unitless = densities_unitless[i][8];
        integralcharge = integralcharge + h*(pow(r_unitless,4.0)*charge_density_r_unitless);
        integralweak = integralweak + h*(pow(r_unitless,4.0)*weak_density_r_unitless);
    }
    Radii2_C_W_fm2[0] = pow(r0_fm,2.0)*(4.0*pi/Z)*integralcharge;
    Radii2_C_W_fm2[1] = pow(r0_fm,2.0)*(4.0*pi/Qwk)*integralweak;
}

// returns exit code -1 if no states are found and 0 if success and -2 if no convergence
// densities form of (r,psn,psp,pvn,pvp,ptn,ptp,pch,pwk)
int hartree_method(double fin_couplings[16], int A, int Z, int iterations, int gridsize, int meson_iterations, double Observables[7], double convergence_help, bool print_densities, bool print_meson_fields) {
    double R_fm = convergence_help*1.05*pow(A,1.0/3.0); // Wood saxon parameter
    double r_init_fm = 1e-5;    // initial starting point 
    double r_final_fm = 20.0;   // final point
    double a_fm = 0.65; // wood saxon parameter
    double** densities_svtnp_unitless; double** meson_fields_unitless;
    int ncols_dens = 11; // (r,sn,sp,vn,vp,tn,tp,ch,wk,dtn,dtp)
    int ncols_meson = 8; // (r,sigma,omega,rho,delta,coulomb,dsigma,domega)
    int npoints_meson = gridsize;
    double BA_mev_last = 0;
    int exit_code = 0;
    int interation_limit = 100;
    int iter_MAX = 500; int count = 0;
    string cont;
    double gs2, gw2, gp2, gd2, kappa, lambda, zeta, xi, lambda_v, lambda_s, fw, fp, mSigma_mev, mOmega_mev, mRho_mev, mDelta_mev, qwn, BA_mev, Fch_Fwk; 
    double Radii2_N_P_fm2[2]; double Radii2_C_W_fm2[2];
    if (A==48 && Z==20) {
        qwn = qwn_Ca;
    } else {
        qwn = qwn_Pb;
    }

    gs2 = fin_couplings[0]; gw2 = fin_couplings[1]; gp2 = fin_couplings[2]; gd2 = fin_couplings[3];
    kappa = fin_couplings[4]; lambda = fin_couplings[5]; zeta = fin_couplings[6]; xi = fin_couplings[7]; lambda_v = fin_couplings[8]; lambda_s = fin_couplings[9];
    fw = fin_couplings[10]; fp = fin_couplings[11];
    mSigma_mev = fin_couplings[12]; mOmega_mev = fin_couplings[13]; mRho_mev = fin_couplings[14]; mDelta_mev = fin_couplings[15];
    fw = fw/sqrt(gw2);
    fp = fp/sqrt(gp2);
    
    // Create Initial Meson Fields
    //dm2.importdata("meson_fields_pre.txt",meson_fields_unitless);
    dm2.create(meson_fields_unitless,npoints_meson,ncols_meson);
    init_meson_r(npoints_meson,meson_fields_unitless,r_init_fm,r_final_fm);
    init_meson(npoints_meson,41.0,R_fm,a_fm,meson_fields_unitless,r_init_fm,r_final_fm,gs2,1);
    init_meson(npoints_meson,25.0,R_fm,a_fm,meson_fields_unitless,r_init_fm,r_final_fm,gw2,2);
    init_meson(npoints_meson,-0.5,R_fm,a_fm,meson_fields_unitless,r_init_fm,r_final_fm,gp2,3);
    init_meson(npoints_meson,-0.5,R_fm,a_fm,meson_fields_unitless,r_init_fm,r_final_fm,gd2,4);
    init_coulomb(npoints_meson,R_fm,meson_fields_unitless,r_init_fm,r_final_fm,Z,5);
    meson_der(npoints_meson,meson_fields_unitless,2,6);
    meson_der(npoints_meson,meson_fields_unitless,3,7);
    
    // Get the wave functions and energy spectrum
    energy_spectrum_proton(meson_fields_unitless,npoints_meson,A,fw,fp,ncols_meson);
    energy_spectrum_neutron(meson_fields_unitless,npoints_meson,A,fw,fp,ncols_meson);
    exit_code = shell_fill(A,Z,0,3);
    if (exit_code == -1) {
        dm2.cleanup(meson_fields_unitless,npoints_meson);
        return -1;
    }

    // Start the self consistent hartree method
    for (int i=0; i<iterations; ++i) {
        get_densities(meson_fields_unitless,A,"neutron_spectrum.txt","proton_spectrum.txt",densities_svtnp_unitless,npoints_meson, ncols_dens,fw,fp,ncols_meson);
        get_nonlinear_meson_fields(meson_fields_unitless,npoints_meson,A,densities_svtnp_unitless,npoints_meson,ncols_dens,1,3,2,4,gs2,gw2,gp2,gd2,kappa,lambda,zeta,xi,lambda_v,lambda_s,fw,fp,mSigma_mev,mOmega_mev,mRho_mev,mDelta_mev,gridsize,meson_iterations,ncols_meson);
        BA_mev = get_BA(meson_fields_unitless,densities_svtnp_unitless,"neutron_spectrum.txt","proton_spectrum.txt",npoints_meson,npoints_meson,ncols_dens,A,kappa,lambda,zeta,xi,lambda_v,lambda_s,fw,fp);
    
        if (fabs(BA_mev_last-BA_mev)<1e-6) {
            break;
        } else {
            iterations = iterations + 1;
        }

        if (iterations>iter_MAX) {
            cout << "hartree did not converge after " <<  iterations << " iterations" << endl;
            dm2.cleanup(densities_svtnp_unitless,npoints_meson);
            dm2.cleanup(meson_fields_unitless,npoints_meson);
            return -2;
        }

        // kill hartree if too many iterations
        if (iterations > interation_limit) {
            if (fabs(BA_mev)>10) {              // check if overbound
                cout << "Way overbound" << endl;
                dm2.cleanup(densities_svtnp_unitless,npoints_meson);
                dm2.cleanup(meson_fields_unitless,npoints_meson);
                return -3;
            }
            if (fabs(BA_mev-BA_mev_last) < 0.1/pow(10.0,count)) {   // check if its even converging
                interation_limit = interation_limit + 100;
                count = count + 1;
                cout << "Adding 100 more iterations" << endl;
            } else {
                dm2.cleanup(densities_svtnp_unitless,npoints_meson);
                dm2.cleanup(meson_fields_unitless,npoints_meson);
                cout << "hartree did not converge after " <<  iterations << " iterations" << endl;
                return -2;
            }
        }
        BA_mev_last = BA_mev;

        // calculate wave functions
        #pragma omp parallel sections
        {
            #pragma omp section
                energy_spectrum_neutron(meson_fields_unitless,npoints_meson,A,fw,fp,ncols_meson);
            #pragma omp section
                energy_spectrum_proton(meson_fields_unitless,npoints_meson,A,fw,fp,ncols_meson);
        }

        // fill nucleon shells
        exit_code = shell_fill(A,Z,0,3);
        if (exit_code == -1) {
            dm2.cleanup(meson_fields_unitless,npoints_meson);
            dm2.cleanup(densities_svtnp_unitless,npoints_meson);
            return -1;
        }   
        //cout << setprecision(8) << "iteration " << i+1 << "  " << fabs(BA_mev) << "  " << densities_svtnp_unitless[0][4]/pow(r0_fm,3.0) << endl;
        dm2.cleanup(densities_svtnp_unitless,npoints_meson);
        
    }

    get_densities(meson_fields_unitless,A,"neutron_spectrum.txt","proton_spectrum.txt",densities_svtnp_unitless,npoints_meson,ncols_dens,fw,fp,ncols_meson);
    
    #pragma omp parallel sections
    {
        #pragma omp section
            Fch_Fwk = get_WEAK_CHARGE_densities_v2(densities_svtnp_unitless,Z,A,ncols_dens,npoints_meson,qwn);
        #pragma omp section
            BA_mev = get_BA(meson_fields_unitless,densities_svtnp_unitless,"neutron_spectrum.txt","proton_spectrum.txt",npoints_meson,npoints_meson,ncols_dens,A,kappa,lambda,zeta,xi,lambda_v,lambda_s,fw,fp);
    }

    #pragma omp parallel sections
    {
        #pragma omp section   
            get_nucleon_radii(densities_svtnp_unitless,npoints_meson,ncols_dens,A,Z,Radii2_N_P_fm2);
        #pragma omp section
            get_weak_charge_radii(densities_svtnp_unitless,npoints_meson,ncols_dens,A,Z,Radii2_C_W_fm2,qwn);
    }
    
    cout << A << "  " << Z << endl;
    cout << "BA: " << BA_mev << endl;
    cout << "Neutron Radius: " << sqrt(Radii2_N_P_fm2[0]) << endl;
    cout << "Proton Radius: " << sqrt(Radii2_N_P_fm2[1]) << endl;
    cout << "Charge Radius: " << sqrt(Radii2_C_W_fm2[0]) << endl;
    cout << "Weak Radius: " << sqrt(Radii2_C_W_fm2[1]) << endl;
    cout << "Rn - Rp: " << sqrt(Radii2_N_P_fm2[0]) - sqrt(Radii2_N_P_fm2[1]) << endl;
    cout << "Rwk - Rch: " << sqrt(Radii2_C_W_fm2[1]) - sqrt(Radii2_C_W_fm2[0]) << endl;
    cout << "Fch - Fwk: " << Fch_Fwk << endl;
    cout << "Rch: " << sqrt(Radii2_N_P_fm2[1] + pow(0.84,2.0)) << endl;
    cout << "-----------------------------------" << endl;
    
    Observables[0] = BA_mev; Observables[1] = sqrt(Radii2_N_P_fm2[0]);
    Observables[2] = sqrt(Radii2_N_P_fm2[1]); Observables[3] = sqrt(Radii2_N_P_fm2[1] + pow(0.84,2.0));
    Observables[4] = sqrt(Radii2_C_W_fm2[1]); Observables[5] = Fch_Fwk;
    Observables[6] = densities_svtnp_unitless[0][7]*pow(r0_fm,-3.0);

    // print meson fields
    if (print_meson_fields == true) {
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,0,r0_fm);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,1,enscale_mev);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,2,enscale_mev);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,3,enscale_mev);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,4,enscale_mev);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,5,enscale_mev);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,6,1.0);
        dm2.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,7,1.0);
        dm2.print(meson_fields_unitless,npoints_meson,ncols_meson,true,"meson_fields.txt"); 
    }
    
    // print densities
    if (print_densities == true) {
        dm2.convert_array(densities_svtnp_unitless,npoints_meson,ncols_dens,0,r0_fm);
        for (int i=1; i<ncols_dens; ++i) {
            dm2.convert_array(densities_svtnp_unitless,npoints_meson,ncols_dens,i,pow(r0_fm,-3.0));
        }

        string outfile = "densities" + to_string(A) + "," + to_string(Z) +".txt";
        dm2.print(densities_svtnp_unitless,npoints_meson,ncols_dens,true,outfile);
    }
    
    dm2.cleanup(densities_svtnp_unitless,npoints_meson);
    dm2.cleanup(meson_fields_unitless,npoints_meson);
    subtract_restmass("proton_spectrum.txt", "neutron_spectrum.txt");
    return 0;
}