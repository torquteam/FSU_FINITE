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

// Generic?
//const double qwp = 0.0713;  // weak vector-charge of the proton
//const double qwn = -0.9821; // weak vector-charge of the neutron

// 208Pb
const double qwp = 0.0713;  // weak vector-charge of the proton
const double qwn = -0.9821; // weak vector-charge of the neutron

// 48Ca
//const double qwp = 0.0713;  // weak vector-charge of the proton
//const double qwn = -0.9795; // weak vector-charge of the neutron

data2 dm2;


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
    double vdens_n = 2.0*nfrac*(2.0*jn+1.0)/(4.0*pi*pow(r_unitless,2.0))*(Fn_r_unitless*Gn_r_unitless);
    return vdens_n;
}

// proton tensor density given proton wave functions
double proton_tensordens(double jp, double pfrac, double Ap_r_unitless, double Bp_r_unitless, double r_unitless) {
    double vdens_p = 2.0*pfrac*(2.0*jp+1.0)/(4.0*pi*pow(r_unitless,2))*(Ap_r_unitless*Bp_r_unitless);
    return vdens_p;
}

// obtain a grid of the initial sigma field (r,sigma) (ouput is unitless)
void init_sigma(int npoints, double V_mev, double R_fm, double a_fm, double** &array, double r_init_fm, double r_final_fm, double gs2) {
    double a_unitless = a_fm/r0_fm;
    double R_unitless = R_fm/r0_fm;
    double r_init_unitless = r_init_fm/r0_fm;
    double r_final_unitless = r_final_fm/r0_fm;
    double V_unitless = V_mev/enscale_mev;

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;

    V_unitless = V_unitless*sqrt(gs2);
    // fill grid for sigma field
    for (int i=0; i<npoints; ++i) {
        array[i][0] = r_unitless;
        array[i][1] = V_unitless/(1.0+exp((r_unitless-R_unitless)/a_unitless));
        r_unitless = r_unitless+stp;
    }
}

// obtain a grid of the initial omega field (r,omega) (ouput is unitless)
void init_omega(int npoints, double V_mev, double R_fm, double a_fm, double** &array, double r_init_fm, double r_final_fm, double gw2) {
    double a_unitless = a_fm/r0_fm;
    double R_unitless = R_fm/r0_fm;
    double r_init_unitless = r_init_fm/r0_fm;
    double r_final_unitless = r_final_fm/r0_fm;
    double V_unitless = V_mev/enscale_mev;

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;

    V_unitless = V_unitless*sqrt(gw2);
    // fill grid for omega field
    for (int i=0; i<npoints; ++i) {
        array[i][0] = r_unitless;
        array[i][1] = V_unitless/(1.0+exp((r_unitless-R_unitless)/a_unitless));
        r_unitless = r_unitless+stp;
    }
}

// obtain a grid of the initial rho field (r,omega) (ouput is unitless)
void init_rho(int npoints, double V_mev, double R_fm, double a_fm, double** &array, double r_init_fm, double r_final_fm, double gp2) {
    double a_unitless = a_fm/r0_fm;
    double R_unitless = R_fm/r0_fm;
    double r_init_unitless = r_init_fm/r0_fm;
    double r_final_unitless = r_final_fm/r0_fm;
    double V_unitless = V_mev/enscale_mev;

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;

    V_unitless = V_unitless*sqrt(gp2);
    // fill grid for rho field
    for (int i=0; i<npoints; ++i) {
        array[i][0] = r_unitless;
        array[i][1] = V_unitless/(1.0+exp((r_unitless-R_unitless)/a_unitless));
        r_unitless = r_unitless+stp;
    }
}

// obtain a grid of the initial coulomb field (r,V(r)) (ouput is unitless)
void init_coulomb(int npoints, double R_fm, double** &array, double r_init_fm, double r_final_fm, int Z) {
    double R_unitless = R_fm/r0_fm;   // convert unitless
    double r_init_unitless = r_init_fm/r0_fm;     // convert unitless
    double r_final_unitless = r_final_fm/r0_fm;   // convert unitless

    double stp = (r_final_unitless-r_init_unitless)/(npoints-1);
    double r_unitless = r_init_unitless;
    int i=0;

    // fill grid for coulomb field
    while (r_unitless<R_unitless) {
        array[i][0] = r_unitless;
        array[i][1] = Z/(4.0*pi*e0_unitless)/(2.0*pow(R_unitless,3))*(3.0*pow(R_unitless,2) - pow(r_unitless,2));
        i = i+1;
        r_unitless = r_unitless+stp;
    }
    while (i<npoints) {
        array[i][0] = r_unitless;
        array[i][1] = Z/(4.0*pi*e0_unitless)*1.0/(r_unitless);
        i = i+1;
        r_unitless = r_unitless+stp;
    }
}

// input is unitless (r/r0, en/Econv)
double dGdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double Fn_r_unitless, double Gn_r_unitless) {
    double res;
    double mN_unitless = mN_mev/enscale_mev;
    res = -r0_fm*enscale_mev*fm_to_inversemev*(En_unitless - gw_omega_r_unitless + 0.5*gp_rho_r_unitless - mN_unitless + gs_sigma_r_unitless)*Fn_r_unitless - alpha*Gn_r_unitless/r_unitless;
    return res;
}

// input is unitless (r/r0, en/Econv)
double dFdr(double r_unitless, double alpha, double En_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double Fn_r_unitless, double Gn_r_unitless) {
    double res;
    double mN_unitless = mN_mev/enscale_mev;
    res = -r0_fm*enscale_mev*fm_to_inversemev*(-En_unitless + gw_omega_r_unitless - 0.5*gp_rho_r_unitless - mN_unitless + gs_sigma_r_unitless)*Gn_r_unitless + alpha*Fn_r_unitless/r_unitless;
    return res;
}

double dBdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless) {
    double res;
    double mP_unitless = mP_mev/enscale_mev;
    res = -r0_fm*enscale_mev*fm_to_inversemev*(Ep_unitless - gw_omega_r_unitless - 0.5*gp_rho_r_unitless - e_coulomb_r_unitless - mP_unitless + gs_sigma_r_unitless)*Ap_r_unitless - alpha*Bp_r_unitless/r_unitless;
    return res;
}

// input is unitless (r/r0, en/Econv)
double dAdr(double r_unitless, double alpha, double Ep_unitless, double gs_sigma_r_unitless, double gw_omega_r_unitless, double gp_rho_r_unitless, double e_coulomb_r_unitless, double Ap_r_unitless, double Bp_r_unitless) {
    double res;
    double mP_unitless = mP_mev/enscale_mev;
    res = -r0_fm*enscale_mev*fm_to_inversemev*(-Ep_unitless + gw_omega_r_unitless + 0.5*gp_rho_r_unitless + e_coulomb_r_unitless - mP_unitless + gs_sigma_r_unitless)*Bp_r_unitless + alpha*Ap_r_unitless/r_unitless;
    return res;
}

void rk4_n(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double FG_unitless[2], double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, int nrows, char direction, double** &Fn_unitless, double** &Gn_unitless) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    double gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless;
    double h = abs(r_final_unitless-r_init_unitless)/(nsteps-1);

    // need initial conditions for Fn_r, Gn_r
    double Fn_r_unitless, Gn_r_unitless;
    double r_unitless = r_init_unitless;
    
    if (direction == 'r') {
        if (alpha>0) {
            Fn_r_unitless = pow(r_unitless,alpha);
            Gn_r_unitless = pow(r_unitless,alpha+1);
        } else {
            Fn_r_unitless = pow(r_unitless,1.0-alpha);
            Gn_r_unitless = pow(r_unitless,-alpha);
        }
        Fn_unitless[0][0] = r_unitless; Fn_unitless[0][1] = Fn_r_unitless;
        Gn_unitless[0][0] = r_unitless; Gn_unitless[0][1] = Gn_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            gs_sigma_r_unitless = gs_sigma_unitless[i-1][1];
            gw_omega_r_unitless = gw_omega_unitless[i-1][1];
            gp_rho_r_unitless = gp_rho_unitless[i-1][1];
            k1 = h*dFdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless, Gn_r_unitless);
            l1 = h*dGdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless, Gn_r_unitless);
            k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0);
            l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0);
            k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0);
            l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0);
            k4 = h*dFdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k3, Gn_r_unitless+l3);
            l4 = h*dGdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k3, Gn_r_unitless+l3);
            Fn_r_unitless = Fn_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Gn_r_unitless = Gn_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = gs_sigma_unitless[i][0];
            Fn_unitless[i][0] = r_unitless; Fn_unitless[i][1] = Fn_r_unitless;
            Gn_unitless[i][0] = r_unitless; Gn_unitless[i][1] = Gn_r_unitless;
        }
    } else {
        Fn_r_unitless = exp(-sqrt(mN_mev*mN_mev - en_unitless*enscale_mev*en_unitless*enscale_mev)*r_unitless*r0_fm*fm_to_inversemev);
        Gn_r_unitless = exp(-sqrt(mN_mev*mN_mev - en_unitless*enscale_mev*en_unitless*enscale_mev)*r_unitless*r0_fm*fm_to_inversemev);
        h = -h;
        Fn_unitless[0][0] = r_unitless; Fn_unitless[0][1] = Fn_r_unitless;
        Gn_unitless[0][0] = r_unitless; Gn_unitless[0][1] = Gn_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            gs_sigma_r_unitless = gs_sigma_unitless[nrows-1-i][1];
            gw_omega_r_unitless = gw_omega_unitless[nrows-1-i][1];
            gp_rho_r_unitless = gp_rho_unitless[nrows-1-i][1];
            k1 = h*dFdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless, Gn_r_unitless);
            l1 = h*dGdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless, Gn_r_unitless);
            k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0);
            l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k1/2.0, Gn_r_unitless+l1/2.0);
            k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0);
            l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k2/2.0, Gn_r_unitless+l2/2.0);
            k4 = h*dFdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k3, Gn_r_unitless+l3);
            l4 = h*dGdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, Fn_r_unitless+k3, Gn_r_unitless+l3);
            Fn_r_unitless = Fn_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Gn_r_unitless = Gn_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = gs_sigma_unitless[nrows-1-i][0];
            Fn_unitless[i][0] = r_unitless; Fn_unitless[i][1] = Fn_r_unitless;
            Gn_unitless[i][0] = r_unitless; Gn_unitless[i][1] = Gn_r_unitless;
        }
    }

    FG_unitless[0] = Fn_r_unitless;
    FG_unitless[1] = Gn_r_unitless;
}

void rk4_p(double r_init_unitless, double r_final_unitless, int nsteps, double alpha, double en_unitless, double AB_unitless[2], double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless, int nrows, char direction, double** &Ap_unitless, double** &Bp_unitless) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    double gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless;
    double h = abs(r_final_unitless-r_init_unitless)/(nsteps-1);

    // need initial conditions for Ap_r, Bp_r
    double Ap_r_unitless, Bp_r_unitless;
    double r_unitless = r_init_unitless;
    
    if (direction == 'r') {
        if (alpha>0) {
            Ap_r_unitless = pow(r_unitless,alpha);
            Bp_r_unitless = pow(r_unitless,alpha+1);
        } else {
            Ap_r_unitless = pow(r_unitless,1.0-alpha);
            Bp_r_unitless = pow(r_unitless,-alpha);
        }
        Ap_unitless[0][0] = r_unitless; Ap_unitless[0][1] = Ap_r_unitless;
        Bp_unitless[0][0] = r_unitless; Bp_unitless[0][1] = Bp_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            gs_sigma_r_unitless = gs_sigma_unitless[i-1][1];
            gw_omega_r_unitless = gw_omega_unitless[i-1][1];
            gp_rho_r_unitless = gp_rho_unitless[i-1][1];
            e_coulomb_r_unitless = e_coulomb_unitless[i-1][1];
            k1 = h*dAdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless, Bp_r_unitless);
            l1 = h*dBdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless, Bp_r_unitless);
            k2 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0);
            l2 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0);
            k3 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0);
            l3 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0);
            k4 = h*dAdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k3, Bp_r_unitless+l3);
            l4 = h*dBdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k3, Bp_r_unitless+l3);
            Ap_r_unitless = Ap_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Bp_r_unitless = Bp_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = gs_sigma_unitless[i][0];
            Ap_unitless[i][0] = r_unitless; Ap_unitless[i][1] = Ap_r_unitless;
            Bp_unitless[i][0] = r_unitless; Bp_unitless[i][1] = Bp_r_unitless;
        }
    } else {
        Ap_r_unitless = exp(-sqrt(mP_mev*mP_mev - en_unitless*enscale_mev*en_unitless*enscale_mev)*r_unitless*r0_fm*fm_to_inversemev);
        Bp_r_unitless = exp(-sqrt(mP_mev*mP_mev - en_unitless*enscale_mev*en_unitless*enscale_mev)*r_unitless*r0_fm*fm_to_inversemev);
        h = -h;
        Ap_unitless[0][0] = r_unitless; Ap_unitless[0][1] = Ap_r_unitless;
        Bp_unitless[0][0] = r_unitless; Bp_unitless[0][1] = Bp_r_unitless;
        for (int i=1; i<nsteps; ++i) {
            gs_sigma_r_unitless = gs_sigma_unitless[nrows-1-i][1];
            gw_omega_r_unitless = gw_omega_unitless[nrows-1-i][1];
            gp_rho_r_unitless = gp_rho_unitless[nrows-1-i][1];
            e_coulomb_r_unitless = e_coulomb_unitless[nrows-1-i][1];
            k1 = h*dAdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless, Bp_r_unitless);
            l1 = h*dBdr(r_unitless, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless, Bp_r_unitless);
            k2 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0);
            l2 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k1/2.0, Bp_r_unitless+l1/2.0);
            k3 = h*dAdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0);
            l3 = h*dBdr(r_unitless+h/2.0, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k2/2.0, Bp_r_unitless+l2/2.0);
            k4 = h*dAdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k3, Bp_r_unitless+l3);
            l4 = h*dBdr(r_unitless+h, alpha, en_unitless, gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, Ap_r_unitless+k3, Bp_r_unitless+l3);
            Ap_r_unitless = Ap_r_unitless + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            Bp_r_unitless = Bp_r_unitless + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r_unitless = gs_sigma_unitless[nrows-1-i][0];
            Ap_unitless[i][0] = r_unitless; Ap_unitless[i][1] = Ap_r_unitless;
            Bp_unitless[i][0] = r_unitless; Bp_unitless[i][1] = Bp_r_unitless;
        }
    }

    AB_unitless[0] = Ap_r_unitless;
    AB_unitless[1] = Bp_r_unitless;
}

// Compute wave functions for given energy and alpha (returns the determinant: -G_ext_unitless*F_int_unitless + F_ext_unitless*G_int_unitless) Fn and Gn must have dimensions (npoints_rk4*2 , 2)
double neutronfield_Solve(double alpha, double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, int nrows_meson, int A, double en_mev, double** &Fn_unitless, double** &Gn_unitless, bool stitch, double r_init_fm, double r_final_fm) {

    // Integration parameters
    double r_left_fm = r_init_fm;
    double r_right_fm = r_final_fm;
    double r_match_fm = 2.0*r0_fm*pow(A,1.0/3.0); // fm
    double FG_unitless[2];
    char direction;

    // convert unitless
    double r_left_unitless = r_left_fm/r0_fm;
    double r_right_unitless = r_right_fm/r0_fm;
    double r_match_unitless = r_match_fm/r0_fm;
    double energy_unitless = en_mev/enscale_mev;

    // create temporary arrays for interior and exterior solutions
    int r_match_row = dm2.findvalue(gs_sigma_unitless,nrows_meson,2,r_match_unitless,0,1e-3); // fm
    int nrows_int = r_match_row + 1;
    int nrows_ext = nrows_meson - nrows_int;
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
    rk4_n(r_left_unitless, gs_sigma_unitless[r_match_row][0], nrows_int, alpha, energy_unitless, FG_unitless, gs_sigma_unitless, gw_omega_unitless, gp_rho_unitless, nrows_meson, direction, Fn_int_unitless, Gn_int_unitless);
    F_int_unitless_rmatch = FG_unitless[0];
    G_int_unitless_rmatch = FG_unitless[1];

    // integrate from infinity to r_match
    direction = 'l';
    rk4_n(r_right_unitless, gs_sigma_unitless[r_match_row][0], nrows_ext, alpha, energy_unitless, FG_unitless, gs_sigma_unitless, gw_omega_unitless, gp_rho_unitless, nrows_meson, direction, Fn_ext_unitless, Gn_ext_unitless);
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
        for (int i=0; i<nrows_ext; ++i) {
            Fn_unitless[i+nrows_int][0] = Fn_ext_unitless[i][0]; Fn_unitless[i+nrows_int][1] = (F_int_unitless_rmatch/F_ext_unitless_rmatch)*Fn_ext_unitless[i][1];
            Gn_unitless[i+nrows_int][0] = Gn_ext_unitless[i][0]; Gn_unitless[i+nrows_int][1] = (G_int_unitless_rmatch/G_ext_unitless_rmatch)*Gn_ext_unitless[i][1];
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
double protonfield_Solve(double alpha, double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless, int nrows_meson, int A, double en_mev, double** &Ap_unitless, double** &Bp_unitless, bool stitch, double r_init_fm, double r_final_fm) {

    // Integration parameters
    double r_left_fm = r_init_fm;
    double r_right_fm = r_final_fm;
    double r_match_fm = 2.0*r0_fm*pow(A,1.0/3.0); // fm
    double AB_unitless[2];
    char direction;

    // convert unitless
    double r_left_unitless = r_left_fm/r0_fm;
    double r_right_unitless = r_right_fm/r0_fm;
    double r_match_unitless = r_match_fm/r0_fm;
    double energy_unitless = en_mev/enscale_mev;

    // create temporary arrays for interior and exterior solutions
    int r_match_row = dm2.findvalue(gs_sigma_unitless,nrows_meson,2,r_match_unitless,0,1e-3); // fm
    int nrows_int = r_match_row + 1;
    int nrows_ext = nrows_meson - nrows_int;
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
    rk4_p(r_left_unitless, r_match_unitless, nrows_int, alpha, energy_unitless, AB_unitless, gs_sigma_unitless, gw_omega_unitless, gp_rho_unitless, e_coulomb_unitless, nrows_meson, direction, Ap_int_unitless, Bp_int_unitless);
    A_int_unitless_rmatch = AB_unitless[0];
    B_int_unitless_rmatch = AB_unitless[1];

    // integrate from infinity to r_match
    direction = 'l';
    rk4_p(r_right_unitless, r_match_unitless, nrows_ext, alpha, energy_unitless, AB_unitless, gs_sigma_unitless, gw_omega_unitless, gp_rho_unitless, e_coulomb_unitless, nrows_meson, direction, Ap_ext_unitless, Bp_ext_unitless);
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
        for (int i=0; i<nrows_ext; ++i) {
            Ap_unitless[i+nrows_int][0] = Ap_ext_unitless[i][0]; Ap_unitless[i+nrows_int][1] = (A_int_unitless_rmatch/A_ext_unitless_rmatch)*Ap_ext_unitless[i][1];
            Bp_unitless[i+nrows_int][0] = Bp_ext_unitless[i][0]; Bp_unitless[i+nrows_int][1] = (B_int_unitless_rmatch/B_ext_unitless_rmatch)*Bp_ext_unitless[i][1];
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
    int nsteps_rk4 = 5000;

    // range and step size to integrate
    double r_init_unitless = min(AF_unitless[0][0], BG_unitless[0][0]);
    double r_final_unitless = min(AF_unitless[nrows-1][0], BG_unitless[nrows-1][0]);
    double h = abs(r_final_unitless-r_init_unitless)/nsteps_rk4;

    // initial conditions
    double r_unitless = r_init_unitless;
    double AF_r_unitless = 0; double BG_r_unitless = 0;
    double psi2 = 0;

    // integrate psi^2 
    for (int i=0; i<nsteps_rk4; ++i) {
        AF_r_unitless = dm2.interpolate(nrows,2,AF_unitless,r_unitless,0,1,true);
        BG_r_unitless = dm2.interpolate(nrows,2,BG_unitless,r_unitless,0,1,true);
        psi2 = psi2 + h*(pow(AF_r_unitless,2.0) + pow(BG_r_unitless,2.0));
        r_unitless = r_unitless + h;
    }
    double norm = sqrt(1.0/psi2); // normalization factor 

    // normalize the wave functions
    for (int i=0; i<nrows; ++i) {
        AF_unitless[i][1] = norm*AF_unitless[i][1];
        BG_unitless[i][1] = norm*BG_unitless[i][1];
    }
}

// returns the energy solution for a given range where the en_determinant changes sign
double energy_bisect_n(double alpha, double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Fn_unitless, double** &Gn_unitless) {
    int count = 0;
    double midx_mev = 0.0;
    double sol_mev = 0.0;
    double en_determinant;
    bool stitch = false;
    double r_init_fm = gs_sigma_unitless[0][0]*r0_fm;
    double r_final_fm = gs_sigma_unitless[nrows_meson-1][0]*r0_fm;
    double min_en_error = abs(neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm));
    double max_en_error = abs(neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm));
    double thresh = 1e-7;
    double error = thresh*min(min_en_error,max_en_error);
    double midy = error*2;

    // bisect the given range until convergence or too many bisections
    while (abs(midy)>error) {
        count = count + 1;
        midx_mev = (en_min_mev+en_max_mev)/2.0;
        en_determinant = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm);
        midy = en_determinant;
        en_determinant = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_min_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm);

        if ((midy*en_determinant)<0) {
            en_max_mev = midx_mev;
        } else {
            en_min_mev = midx_mev;
        }
        sol_mev = midx_mev;

        // Check for divergence
        if (count>75 && abs(midy)>thresh) {
            cout << "No en found for en range: (" << en_min_mev << "," << en_max_mev << ")" << " error= " << abs(midy) << endl;
            exit(0);
        } else if (count>75 && abs(midy)<thresh) {
            stitch=true;
            neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm);
            normalize(Fn_unitless,Gn_unitless,nrows_meson);
            return sol_mev;
        }
    }
    
    stitch=true;
    neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,midx_mev,Fn_unitless,Gn_unitless,stitch,r_init_fm,r_final_fm);
    normalize(Fn_unitless,Gn_unitless,nrows_meson);
    return sol_mev;
}

// returns the energy solution for a given range where the en_determinant changes sign
double energy_bisect_p(double alpha, double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless, int nrows_meson, int A, double en_min_mev, double en_max_mev, double** &Ap_unitless, double** &Bp_unitless) {
    int count = 0;
    double midx_mev = 0.0;
    double sol_mev = 0.0;
    double en_determinant;
    bool stitch = false;
    double r_init_fm = gs_sigma_unitless[0][0]*r0_fm;
    double r_final_fm = gs_sigma_unitless[nrows_meson-1][0]*r0_fm;
    double min_en_error = abs(protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_min_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm));
    double max_en_error = abs(protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_max_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm));
    double thresh = 1e-7;
    double error = thresh*min(min_en_error,max_en_error);
    double midy = error*2;

    // bisect the given range until convergence or too many bisections
    while (abs(midy)>error) {
        count = count + 1;
        midx_mev = (en_min_mev+en_max_mev)/2.0;
        en_determinant = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm);
        midy = en_determinant;
        en_determinant = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_min_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm);

        if ((midy*en_determinant)<0) {
            en_max_mev = midx_mev;
        } else {
            en_min_mev = midx_mev;
        }
        sol_mev = midx_mev;
        //cout << scientific << setprecision(10) << count << "  " << midx_mev << "  " << midy << endl;
        // Check for divergence
        if (count>75 && abs(midy)>thresh) {
            cout << "No en found for en range: (" << en_min_mev << "," << en_max_mev << ")" << " error= " << abs(midy) << endl;
            exit(0);
        } else if (count>75 && abs(midy)<thresh) {
            stitch=true;
            cout << "flag" << endl;
            protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm);
            normalize(Ap_unitless,Bp_unitless,nrows_meson);
            return sol_mev;
        }
    }
    
    stitch=true;
    protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,midx_mev,Ap_unitless,Bp_unitless,stitch,r_init_fm,r_final_fm);
    normalize(Ap_unitless,Bp_unitless,nrows_meson);
    return sol_mev;
}

void energy_spectrum_neutron(double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, int nrows_meson, int A) {
    double en_n_mev;
    double en_n_determinant, en_n_determinant_prev, bound_state_mev;
    double alpha, j;
    double h = 4.0; // energy step
    int node = 0;   // set nodes to zero

    double r_init_fm = gs_sigma_unitless[0][0]*r0_fm;
    double r_final_fm = gs_sigma_unitless[nrows_meson-1][0]*r0_fm;

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
        en_n_determinant_prev = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm);
        while(en_n_mev<mN_mev) {
            en_n_determinant = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm);
            if (en_n_determinant_prev*en_n_determinant < 0) {
                bound_state_mev = energy_bisect_n(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev-h,en_n_mev,Fn_unitless,Gn_unitless);
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
            en_n_determinant_prev = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm);
            while(en_n_mev<mN_mev) {
                en_n_determinant = neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev,Fn_unitless,Gn_unitless,false,r_init_fm,r_final_fm);
                if (en_n_determinant_prev*en_n_determinant < 0) {
                    bound_state_mev = energy_bisect_n(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n_mev-h,en_n_mev,Fn_unitless,Gn_unitless);
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

void energy_spectrum_proton(double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless,int nrows_meson, int A) {
    double en_p_mev;
    double en_p_determinant, en_p_determinant_prev, bound_state_mev;
    double alpha, j;
    double h = 4.0; // energy step
    int node = 0;   // set nodes to zero

    double r_init_fm = gs_sigma_unitless[0][0]*r0_fm;
    double r_final_fm = gs_sigma_unitless[nrows_meson-1][0]*r0_fm;

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
        en_p_determinant_prev = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm);
        while(en_p_mev<mP_mev) {
            en_p_determinant = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm);
            if (en_p_determinant_prev*en_p_determinant < 0) {
                bound_state_mev = energy_bisect_p(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev-h,en_p_mev,Ap_unitless,Bp_unitless);
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
            en_p_determinant_prev = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm);
            while(en_p_mev<mP_mev) {
                en_p_determinant = protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev,Ap_unitless,Bp_unitless,false,r_init_fm,r_final_fm);
                if (en_p_determinant_prev*en_p_determinant < 0) {
                    bound_state_mev = energy_bisect_p(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p_mev-h,en_p_mev,Ap_unitless,Bp_unitless);
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

void shell_fill(int A, int Z, int en_col, int j_col) {
    int N = A-Z;
    double quantum_j, fill_frac;
    int nstates;

    int nrows_N = dm2.rowcount("neutron_spectrum.txt");
    int ncols_N = dm2.colcount("neutron_spectrum.txt");
    double** neutron_array;

    dm2.importdata("neutron_spectrum.txt",neutron_array);
    dm2.sortasc(neutron_array,0,nrows_N,ncols_N);

    ofstream nout("neutron_spectrum.txt");
    int i=0;
    while (N>0) {
        if (i>(nrows_N-1)) {
            cout << "Not enough states to fill all shells: neutron" << endl;
            cout << "States filled: " << A-Z-N << "/" << A-Z << endl;
            exit(0);
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
            nout << fixed << setprecision(10) << neutron_array[i][k] << "  ";
        }
        nout << fixed << setprecision(10) << fill_frac << endl;

        i = i+1;
    }
    
    dm2.cleanup(neutron_array,nrows_N);

    int nrows_P = dm2.rowcount("proton_spectrum.txt");
    int ncols_P = dm2.colcount("proton_spectrum.txt");
    double** proton_array;

    dm2.importdata("proton_spectrum.txt",proton_array);
    dm2.sortasc(proton_array,0,nrows_P,ncols_P);

    ofstream pout("proton_spectrum.txt");
    i=0;
    while (Z>0) {
        if (i>(nrows_P-1)) {
            cout << "Not enough states to fill all shells: proton" << endl;
            cout << "States left: " << Z << endl;
            exit(0);
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
            pout << fixed << setprecision(10) << proton_array[i][k] << "  ";
        }
        pout << fixed << setprecision(10) << fill_frac << endl;

        i = i+1;
    }
    
    dm2.cleanup(proton_array,nrows_P);
}

void get_densities(double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless, int A, string energy_spectrum_neutron, string energy_spectrum_proton, double** &densities_svtnp_unitless, int nrows_meson, int ncols_dens) {
    // set up density matrix
    dm2.create(densities_svtnp_unitless,nrows_meson,ncols_dens);
    dm2.zero(densities_svtnp_unitless,nrows_meson,ncols_dens);

    double r_init_fm = gs_sigma_unitless[0][0]*r0_fm;
    double r_final_fm = gs_sigma_unitless[nrows_meson-1][0]*r0_fm;
    
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
        neutronfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,nrows_meson,A,en_n,Fn_unitless,Gn_unitless,true,r_init_fm,r_final_fm);
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
        protonfield_Solve(alpha,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,nrows_meson,A,en_p,Ap_unitless,Bp_unitless,true,r_init_fm,r_final_fm);
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
    /*
    //Used for analyzing wave functions
    ofstream out("Fn.txt");
    ofstream lout("Gn.txt");
    ofstream pout("Ap.txt");
    ofstream mout("Bp.txt");
    for (int i=0; i<nrows_meson; ++i) {
        out << scientific << setprecision(10) << Fn_unitless_all[i][0][0];
        lout << Gn_unitless_all[i][0][0];
        pout << Ap_unitless_all[i][0][0];
        mout << Bp_unitless_all[i][0][0];
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
    */
    // cleanup
    dm2.cleanup(energy_array_neutron,n_levels);
    dm2.cleanup(energy_array_proton,p_levels);

    dm2.cleanup3d(Fn_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Gn_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Ap_unitless_all,nrows_meson,2);
    dm2.cleanup3d(Bp_unitless_all,nrows_meson,2);
}

double greens_sigma(double r_unitless, double rp_unitless, double mSigma_mev) {
    double mSigma_unitless = mSigma_mev/enscale_mev;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;
    double res;

    if (r_unitless>rp_unitless) {
        res = 1.0/mSigma_unitless*rp_unitless/r_unitless*exp(-mSigma_unitless*r_unitless*conv_r0_en)*sinh(mSigma_unitless*rp_unitless*conv_r0_en);
    } else {
        res = 1.0/mSigma_unitless*rp_unitless/r_unitless*exp(-mSigma_unitless*rp_unitless*conv_r0_en)*sinh(mSigma_unitless*r_unitless*conv_r0_en);
    }
    return res;
}

double greens_omega(double r_unitless, double rp_unitless, double mOmega_mev) {
    double mOmega_unitless = mOmega_mev/enscale_mev;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;
    double res;

    if (r_unitless>rp_unitless) {
        res = 1.0/mOmega_unitless*rp_unitless/r_unitless*exp(-mOmega_unitless*r_unitless*conv_r0_en)*sinh(mOmega_unitless*rp_unitless*conv_r0_en);
    } else {
        res = 1.0/mOmega_unitless*rp_unitless/r_unitless*exp(-mOmega_unitless*rp_unitless*conv_r0_en)*sinh(mOmega_unitless*r_unitless*conv_r0_en);
    }
    return res;
}

double greens_rho(double r_unitless, double rp_unitless, double mRho_mev) {
    double mRho_unitless = mRho_mev/enscale_mev;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;
    double res;

    if (r_unitless>rp_unitless) {
        res = 1.0/mRho_unitless*rp_unitless/r_unitless*exp(-mRho_unitless*r_unitless*conv_r0_en)*sinh(mRho_unitless*rp_unitless*conv_r0_en);
    } else {
        res = 1.0/mRho_unitless*rp_unitless/r_unitless*exp(-mRho_unitless*rp_unitless*conv_r0_en)*sinh(mRho_unitless*r_unitless*conv_r0_en);
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

void meson_interpolator(int npoints_meson, int space, double** &gs_sigma_unitless, double** &gw_omega_unitless, double** &gp_rho_unitless, double** &e_coulomb_unitless, double** densities) {
    // Fill in gaps since integration only fills every n points
    double phi, r_unitless;
    int rmax, rmin;
    for (int i=0; i<npoints_meson; ++i) {
        phi = i%space;
        if (phi != 0) {
            rmax = i-phi+space;
            rmin = i-phi;

            if (rmax == npoints_meson) {
                r_unitless = densities[rmin][0] + (densities[rmin][0] - densities[rmin-1][0]);
                gs_sigma_unitless[i][0] = r_unitless;
                gw_omega_unitless[i][0] = r_unitless;
                gp_rho_unitless[i][0] = r_unitless;
                e_coulomb_unitless[i][0] = r_unitless;

                gs_sigma_unitless[i][1] = gs_sigma_unitless[rmin][1] + (gs_sigma_unitless[rmin][1] - gs_sigma_unitless[rmin-1][1]);
                gw_omega_unitless[i][1] = gw_omega_unitless[rmin][1] + (gw_omega_unitless[rmin][1] - gw_omega_unitless[rmin-1][1]);
                gp_rho_unitless[i][1] = gp_rho_unitless[rmin][1] + (gp_rho_unitless[rmin][1] - gp_rho_unitless[rmin-1][1]);
                e_coulomb_unitless[i][1] = e_coulomb_unitless[rmin][1] + (e_coulomb_unitless[rmin][1] - e_coulomb_unitless[rmin-1][1]);
            } else {
                r_unitless = densities[rmin][0] + (densities[rmax][0] - densities[rmin][0])*phi/(space*1.0);
                gs_sigma_unitless[i][0] = r_unitless;
                gw_omega_unitless[i][0] = r_unitless;
                gp_rho_unitless[i][0] = r_unitless;
                e_coulomb_unitless[i][0] = r_unitless;
            
                gs_sigma_unitless[i][1] = gs_sigma_unitless[rmin][1] + (gs_sigma_unitless[rmax][1] - gs_sigma_unitless[rmin][1])*phi/(space*1.0);
                gw_omega_unitless[i][1] = gw_omega_unitless[rmin][1] + (gw_omega_unitless[rmax][1] - gw_omega_unitless[rmin][1])*phi/(space*1.0);
                gp_rho_unitless[i][1] = gp_rho_unitless[rmin][1] + (gp_rho_unitless[rmax][1] - gp_rho_unitless[rmin][1])*phi/(space*1.0);
                e_coulomb_unitless[i][1] = e_coulomb_unitless[rmin][1] + (e_coulomb_unitless[rmax][1] - e_coulomb_unitless[rmin][1])*phi/(space*1.0);
            }
        }
    }
}

void get_meson_fields(double** &gs_sigma_unitless, double** &gw_omega_unitless, double** &gp_rho_unitless, double** &e_coulomb_unitless, int npoints_meson, int A, double** densities, int nrows_dens, int ncols_dens, int sdens_n_col, int vdens_n_col, int sdens_p_col, int vdens_p_col, double gs2, double gw2, double gp2, double mSigma_mev, double mOmega_mev, double mRho_mev) {
    double k1, k2, k3, k4;
    double rp_unitless, r_unitless;
    double scalardens_n_unitless, vectordens_n_unitless, scalardens_p_unitless, vectordens_p_unitless;
    double e2_charge = 1.4399627/(enscale_mev*r0_fm);
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;

    double rp_final_unitless = densities[nrows_dens-1][0];
    double rp_init_unitless = densities[0][0];

    int npoints = 2000;
    double h_rk4 = (rp_final_unitless-rp_init_unitless)/npoints_meson;
    int space = npoints_meson/npoints;
    
    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(dynamic) private(r_unitless,rp_unitless,scalardens_n_unitless,scalardens_p_unitless,vectordens_n_unitless,vectordens_p_unitless,k1,k2,k3,k4)
    for (int i=0; i<npoints; ++i) {
        r_unitless = densities[space*i][0];
        gs_sigma_unitless[space*i][0] = r_unitless;
        gw_omega_unitless[space*i][0] = r_unitless;
        gp_rho_unitless[space*i][0] = r_unitless;
        e_coulomb_unitless[space*i][0] = r_unitless;
        for (int j=0; j<npoints_meson; ++j) {
            rp_unitless = densities[j][0];
            scalardens_n_unitless = densities[j][sdens_n_col];
            scalardens_p_unitless = densities[j][sdens_p_col];
            vectordens_n_unitless = densities[j][vdens_n_col];
            vectordens_p_unitless = densities[j][vdens_p_col];
            k1 = h_rk4*gs2*1.0/pow(conv_r0_en,2.0)*(scalardens_n_unitless+scalardens_p_unitless)*greens_sigma(r_unitless,rp_unitless,mSigma_mev);
            k2 = h_rk4*gs2*1.0/pow(conv_r0_en,2.0)*(scalardens_n_unitless+scalardens_p_unitless)*greens_sigma(r_unitless,rp_unitless+h_rk4/2.0,mSigma_mev);
            k3 = h_rk4*gs2*1.0/pow(conv_r0_en,2.0)*(scalardens_n_unitless+scalardens_p_unitless)*greens_sigma(r_unitless,rp_unitless+h_rk4/2.0,mSigma_mev);
            k4 = h_rk4*gs2*1.0/pow(conv_r0_en,2.0)*(scalardens_n_unitless+scalardens_p_unitless)*greens_sigma(r_unitless,rp_unitless+h_rk4,mSigma_mev);
            gs_sigma_unitless[space*i][1] = gs_sigma_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

            k1 = h_rk4*gw2*1.0/pow(conv_r0_en,2.0)*(vectordens_n_unitless+vectordens_p_unitless)*greens_omega(r_unitless,rp_unitless,mOmega_mev);
            k2 = h_rk4*gw2*1.0/pow(conv_r0_en,2.0)*(vectordens_n_unitless+vectordens_p_unitless)*greens_omega(r_unitless,rp_unitless+h_rk4/2.0,mOmega_mev);
            k3 = h_rk4*gw2*1.0/pow(conv_r0_en,2.0)*(vectordens_n_unitless+vectordens_p_unitless)*greens_omega(r_unitless,rp_unitless+h_rk4/2.0,mOmega_mev);
            k4 = h_rk4*gw2*1.0/pow(conv_r0_en,2.0)*(vectordens_n_unitless+vectordens_p_unitless)*greens_omega(r_unitless,rp_unitless+h_rk4,mOmega_mev);
            gw_omega_unitless[space*i][1] = gw_omega_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

            k1 = h_rk4*0.5*gp2*1.0/pow(conv_r0_en,2.0)*(vectordens_p_unitless-vectordens_n_unitless)*greens_rho(r_unitless,rp_unitless,mRho_mev);
            k2 = h_rk4*0.5*gp2*1.0/pow(conv_r0_en,2.0)*(vectordens_p_unitless-vectordens_n_unitless)*greens_rho(r_unitless,rp_unitless+h_rk4/2.0,mRho_mev);
            k3 = h_rk4*0.5*gp2*1.0/pow(conv_r0_en,2.0)*(vectordens_p_unitless-vectordens_n_unitless)*greens_rho(r_unitless,rp_unitless+h_rk4/2.0,mRho_mev);
            k4 = h_rk4*0.5*gp2*1.0/pow(conv_r0_en,2.0)*(vectordens_p_unitless-vectordens_n_unitless)*greens_rho(r_unitless,rp_unitless+h_rk4,mRho_mev);
            gp_rho_unitless[space*i][1] = gp_rho_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

            k1 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless);
            k2 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4/2.0);
            k3 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4/2.0);
            k4 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4);
            e_coulomb_unitless[space*i][1] = e_coulomb_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }
    }
    meson_interpolator(npoints_meson,space,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,densities);
}

void get_nonlinear_meson_fields(double** &gs_sigma_unitless, double** &gw_omega_unitless, double** &gp_rho_unitless, double** &e_coulomb_unitless, int npoints_meson, int A, double** densities, int nrows_dens, int ncols_dens, int sdens_n_col, int vdens_n_col, int sdens_p_col, int vdens_p_col, double gs2, double gw2, double gp2, double b, double c, double h, double lambda, double mSigma_mev, double mOmega_mev, double mRho_mev, int gridsize_meson, int meson_iterations) {
    double k1, k2, k3, k4;
    double rp_unitless, r_unitless;
    double scalardens_n_unitless, vectordens_n_unitless, scalardens_p_unitless, vectordens_p_unitless;
    double e2_charge = 1.4399627/(enscale_mev*r0_fm);
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;

    double rp_final_unitless = densities[nrows_dens-1][0];
    double rp_init_unitless = densities[0][0];

    int npoints = gridsize_meson;
    double h_rk4 = (rp_final_unitless-rp_init_unitless)/npoints;
    int space = npoints_meson/npoints;

    // get coulomb field
    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(static,6) private(r_unitless,rp_unitless,vectordens_p_unitless,k1,k2,k3,k4)
    for (int i=0; i<npoints; ++i) {
        r_unitless = densities[space*i][0];
        e_coulomb_unitless[space*i][0] = r_unitless;
        e_coulomb_unitless[space*i][1] = 0.0;
        for (int j=0; j<npoints; ++j) {
            rp_unitless = densities[j*space][0];
            vectordens_p_unitless = densities[j*space][vdens_p_col];

            k1 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless);
            k2 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4/2.0);
            k3 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4/2.0);
            k4 = h_rk4*4.0*pi*e2_charge*(vectordens_p_unitless)*greens_coulomb(r_unitless,rp_unitless+h_rk4);
            e_coulomb_unitless[space*i][1] = e_coulomb_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
        }
    }

    // Nonlinear convergence
    int MAX_ITER = meson_iterations;
    double oldfield_sigma_unitless, oldfield_omega_unitless, oldfield_rho_unitless;
    double mNuc_unitless = mNuc_mev/enscale_mev;
    double** oldfield_sigma; double** oldfield_omega; double** oldfield_rho;
    dm2.create(oldfield_sigma,npoints_meson,2);
    dm2.create(oldfield_omega,npoints_meson,2);
    dm2.create(oldfield_rho,npoints_meson,2);
    dm2.copy_pointer(gs_sigma_unitless,oldfield_sigma,npoints_meson,2);
    dm2.copy_pointer(gw_omega_unitless,oldfield_omega,npoints_meson,2);
    dm2.copy_pointer(gp_rho_unitless,oldfield_rho,npoints_meson,2);

    // Convergence for the Sigma Field
    for (int k=0; k<MAX_ITER; ++k) {
        #pragma omp parallel num_threads(12)
        #pragma omp for schedule(static,6) private(r_unitless,rp_unitless,scalardens_n_unitless,scalardens_p_unitless,oldfield_sigma_unitless,k1,k2,k3,k4)
        for (int i=0; i<npoints; ++i) {
            r_unitless = densities[space*i][0];
            gs_sigma_unitless[space*i][0] = r_unitless;
            gs_sigma_unitless[space*i][1] = 0.0;
            for (int j=0; j<npoints; ++j) {
                rp_unitless = densities[j*space][0];
                scalardens_n_unitless = densities[j*space][sdens_n_col];
                scalardens_p_unitless = densities[j*space][sdens_p_col];
                oldfield_sigma_unitless = oldfield_sigma[j*space][1];

                k1 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gs2*(scalardens_n_unitless+scalardens_p_unitless) - conv_r0_en*b*mNuc_unitless*pow(oldfield_sigma_unitless,2.0)*gs2 - conv_r0_en*c*pow(oldfield_sigma_unitless,3.0)*gs2)*greens_sigma(r_unitless,rp_unitless, mSigma_mev);
                k2 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gs2*(scalardens_n_unitless+scalardens_p_unitless) - conv_r0_en*b*mNuc_unitless*pow(oldfield_sigma_unitless,2.0)*gs2 - conv_r0_en*c*pow(oldfield_sigma_unitless,3.0)*gs2)*greens_sigma(r_unitless,rp_unitless+h_rk4/2.0, mSigma_mev);
                k3 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gs2*(scalardens_n_unitless+scalardens_p_unitless) - conv_r0_en*b*mNuc_unitless*pow(oldfield_sigma_unitless,2.0)*gs2 - conv_r0_en*c*pow(oldfield_sigma_unitless,3.0)*gs2)*greens_sigma(r_unitless,rp_unitless+h_rk4/2.0, mSigma_mev);
                k4 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gs2*(scalardens_n_unitless+scalardens_p_unitless) - conv_r0_en*b*mNuc_unitless*pow(oldfield_sigma_unitless,2.0)*gs2 - conv_r0_en*c*pow(oldfield_sigma_unitless,3.0)*gs2)*greens_sigma(r_unitless,rp_unitless+h_rk4, mSigma_mev);
                gs_sigma_unitless[space*i][1] = gs_sigma_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            }
            
            oldfield_sigma[i*space][1] = gs_sigma_unitless[i*space][1];
        }
    }

    // Convergence for the Omega/Rho Field
    for (int k=0; k<MAX_ITER; ++k) {
        #pragma omp parallel num_threads(12)
        #pragma omp for schedule(static,1) private(r_unitless,rp_unitless,vectordens_n_unitless,vectordens_p_unitless,oldfield_omega_unitless,oldfield_rho_unitless,k1,k2,k3,k4)
        for (int i=0; i<npoints; ++i) {
            r_unitless = densities[space*i][0];
            gw_omega_unitless[space*i][0] = r_unitless;
            gw_omega_unitless[space*i][1] = 0.0;
            gp_rho_unitless[space*i][0] = r_unitless;
            gp_rho_unitless[space*i][1] = 0.0;
            for (int j=0; j<npoints; ++j) {
                rp_unitless = densities[j*space][0];
                vectordens_n_unitless = densities[j*space][vdens_n_col];
                vectordens_p_unitless = densities[j*space][vdens_p_col];
                oldfield_omega_unitless = oldfield_omega[j*space][1];
                oldfield_rho_unitless = oldfield_rho[j*space][1];

                k1 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gw2*(vectordens_n_unitless+vectordens_p_unitless) - conv_r0_en*h*pow(oldfield_omega_unitless,3.0)*gw2 - conv_r0_en*lambda*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless*gw2)*greens_omega(r_unitless,rp_unitless,mOmega_mev);
                k2 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gw2*(vectordens_n_unitless+vectordens_p_unitless) - conv_r0_en*h*pow(oldfield_omega_unitless,3.0)*gw2 - conv_r0_en*lambda*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless*gw2)*greens_omega(r_unitless,rp_unitless+h_rk4/2.0,mOmega_mev);
                k3 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gw2*(vectordens_n_unitless+vectordens_p_unitless) - conv_r0_en*h*pow(oldfield_omega_unitless,3.0)*gw2 - conv_r0_en*lambda*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless*gw2)*greens_omega(r_unitless,rp_unitless+h_rk4/2.0,mOmega_mev);
                k4 = h_rk4*(1.0/pow(conv_r0_en,2.0)*gw2*(vectordens_n_unitless+vectordens_p_unitless) - conv_r0_en*h*pow(oldfield_omega_unitless,3.0)*gw2 - conv_r0_en*lambda*pow(oldfield_rho_unitless,2.0)*oldfield_omega_unitless*gw2)*greens_omega(r_unitless,rp_unitless+h_rk4,mOmega_mev);
                gw_omega_unitless[space*i][1] = gw_omega_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            }
            oldfield_omega[i*space][1] = gw_omega_unitless[i*space][1];

            for (int j=0; j<npoints; ++j) {
                rp_unitless = densities[j*space][0];
                vectordens_n_unitless = densities[j*space][vdens_n_col];
                vectordens_p_unitless = densities[j*space][vdens_p_col];
                oldfield_rho_unitless = oldfield_rho[j*space][1];
                oldfield_omega_unitless = oldfield_omega[j*space][1];

                k1 = h_rk4*(1.0/pow(conv_r0_en,2.0)*0.5*gp2*(vectordens_p_unitless-vectordens_n_unitless) - conv_r0_en*lambda*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless*gp2)*greens_rho(r_unitless,rp_unitless,mRho_mev);
                k2 = h_rk4*(1.0/pow(conv_r0_en,2.0)*0.5*gp2*(vectordens_p_unitless-vectordens_n_unitless) - conv_r0_en*lambda*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless*gp2)*greens_rho(r_unitless,rp_unitless+h_rk4/2.0,mRho_mev);
                k3 = h_rk4*(1.0/pow(conv_r0_en,2.0)*0.5*gp2*(vectordens_p_unitless-vectordens_n_unitless) - conv_r0_en*lambda*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless*gp2)*greens_rho(r_unitless,rp_unitless+h_rk4/2.0,mRho_mev);
                k4 = h_rk4*(1.0/pow(conv_r0_en,2.0)*0.5*gp2*(vectordens_p_unitless-vectordens_n_unitless) - conv_r0_en*lambda*pow(oldfield_omega_unitless,2.0)*oldfield_rho_unitless*gp2)*greens_rho(r_unitless,rp_unitless+h_rk4,mRho_mev);
                gp_rho_unitless[space*i][1] = gp_rho_unitless[space*i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
            }
            oldfield_rho[i*space][1] = gp_rho_unitless[i*space][1];
        }
    }
    meson_interpolator(npoints_meson,space,gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,densities);
    
    dm2.cleanup(oldfield_sigma,npoints_meson);
    dm2.cleanup(oldfield_omega,npoints_meson);
    dm2.cleanup(oldfield_rho,npoints_meson);
}

// return the Binding energy in MeV
double get_BA(double** gs_sigma_unitless, double** gw_omega_unitless, double** gp_rho_unitless, double** e_coulomb_unitless, double** densities_unitless, string n_energies, string p_energies, int npoints_meson, int npoints_densities, int ncols_density, int A, double b, double c, double h, double lambda) {
    double BA_unitless = 0;
    double** n_spectrum;
    double** p_spectrum;

    double mNuc_unitless = mNuc_mev/enscale_mev;
    double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;

    dm2.importdata(n_energies,n_spectrum);
    dm2.importdata(p_energies,p_spectrum);
    double nrows = dm2.rowcount(n_energies);
    double prows = dm2.rowcount(p_energies);

    for (int i=0; i<nrows; ++i) {
        BA_unitless = BA_unitless + (n_spectrum[i][0])*(2.0*n_spectrum[i][3]+1)*n_spectrum[i][5]/enscale_mev;
    }

    for (int i=0; i<prows; ++i) {
        BA_unitless = BA_unitless + (p_spectrum[i][0])*(2.0*p_spectrum[i][3]+1)*p_spectrum[i][5]/enscale_mev;
    }
    
    int nsteps_rk4 = npoints_meson;
    double r_init_unitless = min(densities_unitless[0][0], gs_sigma_unitless[0][0]);
    double r_final_unitless = min(densities_unitless[npoints_densities-1][0], gs_sigma_unitless[npoints_meson-1][0]);
    double h_rk4 = abs(r_final_unitless-r_init_unitless)/(nsteps_rk4-1);
    double r_unitless;
    double integral = 0.0;
    double gs_sigma_r_unitless, gw_omega_r_unitless, gp_rho_r_unitless, e_coulomb_r_unitless, ns_density_r_unitless, ps_density_r_unitless, nv_density_r_unitless, pv_density_r_unitless;
    for (int i=0; i<nsteps_rk4; ++i) {
        r_unitless = gs_sigma_unitless[i][0];
        gs_sigma_r_unitless = gs_sigma_unitless[i][1];
        gw_omega_r_unitless = gw_omega_unitless[i][1];
        gp_rho_r_unitless = gp_rho_unitless[i][1];
        e_coulomb_r_unitless = e_coulomb_unitless[i][1];
        ns_density_r_unitless = densities_unitless[i][1];
        ps_density_r_unitless = densities_unitless[i][2];
        nv_density_r_unitless = densities_unitless[i][3];
        pv_density_r_unitless = densities_unitless[i][4];
        integral = integral + h_rk4*(gs_sigma_r_unitless*(ns_density_r_unitless+ps_density_r_unitless) - 1.0/3.0*b*mNuc_unitless*pow(conv_r0_en,3.0)*pow(gs_sigma_r_unitless,3.0) - 0.5*c*pow(conv_r0_en,3.0)*pow(gs_sigma_r_unitless,4.0)
                                    - gw_omega_r_unitless*(nv_density_r_unitless+pv_density_r_unitless) + 0.5*h*pow(gw_omega_r_unitless,4.0)*pow(conv_r0_en,3.0)
                                    - 0.5*gp_rho_r_unitless*(pv_density_r_unitless-nv_density_r_unitless) + lambda*pow(gw_omega_r_unitless,2.0)*pow(gp_rho_r_unitless,2.0)*pow(conv_r0_en,3.0) - e_coulomb_r_unitless*pv_density_r_unitless)*pow(r_unitless,2);
    }
    //cout << BA_unitless << "  " << integral << endl;
    BA_unitless = BA_unitless + 2.0*pi*integral;

    dm2.cleanup(n_spectrum,nrows);
    dm2.cleanup(p_spectrum,prows);
    return (BA_unitless*enscale_mev - 0.75*41.0*pow(A,-1.0/3.0))/A - mNuc_mev;
}

void get_nucleon_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_N_P_fm2[2]) {
    int nsteps_rk4 = npoints_densities;
    double r_init_unitless = densities_unitless[0][0];
    double r_final_unitless = densities_unitless[npoints_densities-1][0];
    double h = abs(r_final_unitless-r_init_unitless)/(nsteps_rk4-1);
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
        n_array[i][0] = n_array[i][0] - mNuc_mev;
    }

    for (int i=0; i<prows; ++i) {
        p_array[i][0] = p_array[i][0] - mNuc_mev;
    }

    dm2.print(n_array,nrows,6,true,"neutron_spectrum.txt");
    dm2.print(p_array,prows,6,true,"proton_spectrum.txt");

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

void WEAK_formfactors(double WGE_pn[2], double WGM_pn[2], double q) {
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

double weak_formfactor(double q_unitless, double** densities_svtnp_unitless, int nrows, int A, int Z) {
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

    WEAK_formfactors(WGE_pn,WGM_pn,q_ifm);
    vector_formfactor(densities_svtnp_unitless,q_unitless,nrows,vFF_pn);
    tensor_formfactor(densities_svtnp_unitless,q_unitless,nrows,tFF_pn);
    wFF_p = WGE_pn[0]*vFF_pn[0] + (WGM_pn[0] - WGE_pn[0])/(1.0+tau)*(tau*vFF_pn[0] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[0]);
    wFF_n = WGE_pn[1]*vFF_pn[1] + (WGM_pn[1] - WGE_pn[1])/(1.0+tau)*(tau*vFF_pn[1] + 0.5*q_ifm*hbar_mevfm/mNuc_mev*tFF_pn[1]);
    wFF = wFF_p + wFF_n;

    double Qwk = Z*qwp + (A-Z)*qwn;

    return wFF/Qwk;
}

// version using the txt file
void get_WEAK_CHARGE_densities_v1(string densities, int Z, int A){
    double** densities_svtnp_unitless;
    dm2.importdata(densities,densities_svtnp_unitless);
    int nrows = dm2.rowcount(densities);
    int ncols = dm2.colcount(densities);

    ofstream out("density_new.txt");

    double r_unitless, q_unitless, charge_dens_unitless, weak_dens_unitless;
    double Qwk = Z*qwp + (A-Z)*qwn;
    cout << setprecision(10) << Qwk << endl;

    // get charge and weak form factors first
    int npoints = nrows;
    double** q_charge_weak_factors;
    dm2.create(q_charge_weak_factors,npoints,3);

    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(static,1) private(q_unitless)
    for (int i=0; i<npoints; ++i) {
        q_unitless = 1e-15 + 0.005*i;
        q_charge_weak_factors[i][0] = q_unitless;
        q_charge_weak_factors[i][1] = charge_formfactor(q_unitless,densities_svtnp_unitless,nrows,Z);
        q_charge_weak_factors[i][2] = weak_formfactor(q_unitless,densities_svtnp_unitless,nrows,A,Z);
        
        //cout << setprecision(10) << pow(q_unitless/r0_fm,2.0) << "  " << q_charge_weak_factors[i][1] << "  " << q_charge_weak_factors[i][2] << endl;
    }

    // inverse fourier transform
    double q_unitless_a, q_unitless_b, charge_factor_a, charge_factor_b, fa, fb, q_unitless_ab2, charge_factor_ab2, fab2, weak_factor_a, weak_factor_b, weak_factor_ab2, ga, gb, gab2;
    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(q_unitless_a,q_unitless_b,q_unitless_ab2,charge_factor_a,charge_factor_b,charge_factor_ab2,fa,fb,fab2,r_unitless,charge_dens_unitless,weak_factor_a, weak_factor_b, weak_factor_ab2, ga, gb, gab2, weak_dens_unitless)
    for (int i=0; i<nrows; ++i) {
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

        #pragma omp ordered
        for (int k=0; k<ncols; ++k) {
            out << setprecision(10) << densities_svtnp_unitless[i][k] << "  ";
        }
        out << charge_dens_unitless << "  " << weak_dens_unitless << endl;
    }
    
    dm2.cleanup(densities_svtnp_unitless,nrows);
    dm2.cleanup(q_charge_weak_factors,npoints);
}

// version using the array
void get_WEAK_CHARGE_densities_v2(double** &densities_svtnp_unitless, int Z, int A, int ncols_dens, int nrows_dens){

    double r_unitless, q_unitless, charge_dens_unitless, weak_dens_unitless;
    double Qwk = Z*qwp + (A-Z)*qwn;

    // get charge and weak form factors first
    int npoints = nrows_dens;
    double** q_charge_weak_factors;
    dm2.create(q_charge_weak_factors,npoints,3);

    #pragma omp parallel num_threads(12)
    #pragma omp for schedule(static,1) private(q_unitless)
    for (int i=0; i<npoints; ++i) {
        q_unitless = 1e-15 + 0.01*i;
        q_charge_weak_factors[i][0] = q_unitless;
        q_charge_weak_factors[i][1] = charge_formfactor(q_unitless,densities_svtnp_unitless,nrows_dens,Z);
        q_charge_weak_factors[i][2] = weak_formfactor(q_unitless,densities_svtnp_unitless,nrows_dens,A,Z);
        
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
    
    dm2.cleanup(q_charge_weak_factors,npoints);
}

void get_weak_charge_radii(double** densities_unitless, int npoints_densities, int ncols_density, int A, int Z, double Radii2_C_W_fm2[2]) {
    int nsteps_rk4 = npoints_densities;
    double r_init_unitless = densities_unitless[0][0];
    double r_final_unitless = densities_unitless[npoints_densities-1][0];
    double h = abs(r_final_unitless-r_init_unitless)/(nsteps_rk4-1);
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

void hartee_method(double gs2, double gw2, double gp2, double b, double c, double h, double lambda, double mSigma_mev, double mOmega_mev, double mRho_mev, int A, int Z, int iterations, int gridsize, int gridsize_meson, int meson_iterations, double Observables[5]) {
    double R_fm = 1.2*pow(A,1.0/3.0);
    double r_init_fm = 1e-5;
    double r_final_fm = 5.0*R_fm;
    double a_fm = 0.65;
    double** gs_sigma_unitless; double** gw_omega_unitless; double** gp_rho_unitless; double** e_coulomb_unitless;
    double** densities_svtnp_unitless;
    int ncols_dens = 9;
    int npoints_meson = gridsize;
    double BA_mev;

    // Create Initial Meson Fields
    dm2.create(gs_sigma_unitless,npoints_meson,2);
    dm2.create(gw_omega_unitless,npoints_meson,2);
    dm2.create(gp_rho_unitless,npoints_meson,2);
    dm2.create(e_coulomb_unitless,npoints_meson,2); 
    init_sigma(npoints_meson,41,R_fm,a_fm,gs_sigma_unitless,r_init_fm,r_final_fm,gs2);
    init_omega(npoints_meson,25,R_fm,a_fm,gw_omega_unitless,r_init_fm,r_final_fm,gw2);
    init_rho(npoints_meson,-0.5,R_fm,a_fm,gp_rho_unitless,r_init_fm,r_final_fm,gp2);
    init_coulomb(npoints_meson,R_fm,e_coulomb_unitless,r_init_fm,r_final_fm,Z);

    // Get the wave functions and energy spectrum
    energy_spectrum_proton(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,npoints_meson,A);
    energy_spectrum_neutron(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,npoints_meson,A);
    shell_fill(A,Z,0,3);

    // Start the self consistent hartree method
    for (int i=0; i<iterations; ++i) {
        get_densities(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,A,"neutron_spectrum.txt","proton_spectrum.txt",densities_svtnp_unitless,npoints_meson, ncols_dens);
        get_nonlinear_meson_fields(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,npoints_meson,A,densities_svtnp_unitless,npoints_meson,ncols_dens,1,3,2,4,gs2,gw2,gp2,b,c,h,lambda,mSigma_mev,mOmega_mev,mRho_mev,gridsize_meson,meson_iterations);
        BA_mev = get_BA(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,densities_svtnp_unitless,"neutron_spectrum.txt","proton_spectrum.txt",npoints_meson,npoints_meson,ncols_dens,A,b,c,h,lambda);
        #pragma omp parallel sections
        {
            #pragma omp section
                energy_spectrum_neutron(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,npoints_meson,A);
            #pragma omp section
                energy_spectrum_proton(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,npoints_meson,A);
        }
        shell_fill(A,Z,0,3);
        dm2.cleanup(densities_svtnp_unitless,npoints_meson);
        cout << "iteration " << i+1 << "  " << abs(BA_mev) << endl;
    }

    get_densities(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,A,"neutron_spectrum.txt","proton_spectrum.txt",densities_svtnp_unitless,npoints_meson,ncols_dens);
    get_WEAK_CHARGE_densities_v2(densities_svtnp_unitless,Z,A,ncols_dens,npoints_meson);
    //dm2.print(densities_svtnp_unitless,npoints_meson,ncols_dens,true,"density.txt");
    BA_mev = get_BA(gs_sigma_unitless,gw_omega_unitless,gp_rho_unitless,e_coulomb_unitless,densities_svtnp_unitless,"neutron_spectrum.txt","proton_spectrum.txt",npoints_meson,npoints_meson,ncols_dens,A,b,c,h,lambda);
    double Radii2_N_P_fm2[2]; double Radii2_C_W_fm2[2];
    get_nucleon_radii(densities_svtnp_unitless,npoints_meson,ncols_dens,A,Z,Radii2_N_P_fm2);
    get_weak_charge_radii(densities_svtnp_unitless,npoints_meson,ncols_dens,A,Z,Radii2_C_W_fm2);
    
    cout << "BA: " << BA_mev << endl;
    cout << "Neutron Radius: " << sqrt(Radii2_N_P_fm2[0]) << endl;
    cout << "Proton Radius: " << sqrt(Radii2_N_P_fm2[1]) << endl;
    cout << "Charge Radius: " << sqrt(Radii2_C_W_fm2[0]) << endl;
    cout << "Weak Radius: " << sqrt(Radii2_C_W_fm2[1]) << endl;
    cout << "Rn - Rp: " << sqrt(Radii2_N_P_fm2[0]) - sqrt(Radii2_N_P_fm2[1]) << endl;
    cout << "Rch - Rwk: " << sqrt(Radii2_C_W_fm2[0]) - sqrt(Radii2_C_W_fm2[1]) << endl;

    Observables[0] = BA_mev; Observables[1] = sqrt(Radii2_N_P_fm2[0]);
    Observables[2] = sqrt(Radii2_N_P_fm2[1]); Observables[3] = sqrt(Radii2_C_W_fm2[0]);
    Observables[4] = sqrt(Radii2_C_W_fm2[1]);
    dm2.cleanup(densities_svtnp_unitless,npoints_meson);
    dm2.cleanup(gs_sigma_unitless,npoints_meson);
    dm2.cleanup(gw_omega_unitless,npoints_meson);
    dm2.cleanup(gp_rho_unitless,npoints_meson);
    dm2.cleanup(e_coulomb_unitless,npoints_meson);
    subtract_restmass("proton_spectrum.txt", "neutron_spectrum.txt");
}