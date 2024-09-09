#include "NumMethods.hpp"
#include "infinitematter.hpp"
#include "minpack.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;
data2 dmi;
nummeth nmi;
tools tool;
bulks bulk;
equationofstate eos;

// coupling array [gsoms2,gwomw2,gpomp2,gdomd2,b,c,h,lambda]

// Constants
const double pi = 4.0*atan(1.0);
const double mP = 939.0; //938.27231;    // Mass of proton (MeV)
const double mN = 939.0; //939.56542052;
const double mNuc = (mP+mN)/2.0;
const double mE = 0.511;
const double mMU = 105.7;

//###########################################################
// Tools class
// ##########################################################
double tools :: scalardens(double k, double mstar) {
    double en = sqrt(pow(k,2.0)+pow(mstar,2.0));
    double integral = mstar/2.0*( k*en - pow(mstar,2.0)*log( abs((k+en)/mstar) ) );
    return integral/pow(pi,2.0);
}

double tools :: endensity_integral(double k, double mstar) {
    double x = k/mstar;
    double y = sqrt(1.0 + pow(x,2.0));
    double integral = pow(mstar,4.0)/(8.0*pow(pi,2.0))*(pow(x,3.0)*y + pow(y,3.0)*x  - log(x+y));
    return integral;
}

double tools :: pr_integral(double k, double mstar) {
    double x = k/mstar;
    double y = sqrt(1.0 + pow(x,2.0));
    double integral = pow(mstar,4.0)/(8.0*pow(pi,2.0))*(2.0/3.0*pow(x,3.0)*y - x*y + log(x+y));
    return integral;
}

// integral of k^4/rp^3 dk (comes from derivative of the scalar dens)
double tools :: common_integral(double k, double mstar) {
    double en = sqrt(pow(k,2.0)+pow(mstar,2.0));
    double integral = (0.5*pow(k,3.0) + 3.0/2.0*k*pow(mstar,2.0))/en + 3.0*pow(mstar,2.0)/2.0*log(abs(mstar/(k+en)));
    return integral;
}

// Field Equation for g_sigma sigma (used for bisection method to solve for g_sigma*sigma)
double tools :: gss_FE(double kf, double gss, double gdd, double couplings[10], double t) {
    double res,mstarp,mstarn,kfn,kfp,sdensn,sdensp;
    double gsoms2,kappa,lambda,lambda_s;

    gsoms2 = couplings[0];
    kappa = couplings[4];
    lambda = couplings[5];
    lambda_s = couplings[9];
    
    mstarp = mP-gss-0.5*gdd; mstarn = mN-gss+0.5*gdd;     // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // fermi momenta of neutron and proton
    sdensp = scalardens(kfp,mstarp);     // scalar dens of proton
    sdensn = scalardens(kfn,mstarn);    // scalar dens of neutron
    
    res = gsoms2*(sdensp + sdensn - 0.5*kappa*pow(gss,2.0) - 1.0/6.0*lambda*pow(gss,3.0) - 2.0*lambda_s*gss*pow(gdd,2.0)) - gss;
    return res;
}

// Field Equation for g_delta delta (used for bisection method to solve for g_delta*delta)
double tools :: gdd_FE(double kf, double gss, double gdd, double couplings[10], double t) {
    double res,mstarp,mstarn,kfn,kfp,sdensn,sdensp;
    double gdomd2, lambda_s;

    gdomd2 = couplings[3];
    lambda_s = couplings[9];

    mstarp = mP-gss-0.5*gdd; mstarn = mN-gss+0.5*gdd;     // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // fermi momenta of neutron and proton
    sdensp = scalardens(kfp,mstarp);     // scalar dens of proton
    sdensn = scalardens(kfn,mstarn);    // scalar dens of neutron

    res = gdomd2*(0.5*sdensp - 0.5*sdensn - 2.0*lambda_s*gdd*pow(gss,2.0)) - gdd;
    return res;
}

//########################################################################
//########################################################################
//########################################################################

// Field Equation for g_omega omega (used for bisection method to solve for g_omega*omega)
double gww_FE(double kf, double gww, double gpp, double couplings[10], double t) {
    double res,kfn,kfp,vdensn,vdensp;
    double gwomw2, lambda_v, zeta;

    gwomw2 = couplings[1];
    zeta = couplings[6];
    lambda_v = couplings[8];

    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // fermi momenta of neutron and proton
    vdensp = 1.0/(3.0*pow(pi,2.0))*pow(kfp,3.0);     // scalar dens of proton
    vdensn = 1.0/(3.0*pow(pi,2.0))*pow(kfn,3.0);    // scalar dens of neutron

    res = gwomw2*(vdensp + vdensn - 2.0*lambda_v*gww*pow(gpp,2.0) - 1.0/6.0*zeta*pow(gww,3.0)) - gww;
    return res;
}

// Field Equation for g_omega omega (used for bisection method to solve for g_omega*omega)
double gpp_FE(double kf, double gww, double gpp, double couplings[10], double t) {
    double res,kfn,kfp,vdensn,vdensp;
    double gpomp2, xi, lambda_v;

    gpomp2 = couplings[2];
    xi = couplings[7];
    lambda_v = couplings[8];

    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // fermi momenta of neutron and proton
    vdensp = 1.0/(3.0*pow(pi,2.0))*pow(kfp,3.0);     // scalar dens of proton
    vdensn = 1.0/(3.0*pow(pi,2.0))*pow(kfn,3.0);    // scalar dens of neutron

    res = gpomp2*(0.5*vdensp - 0.5*vdensn - 2.0*lambda_v*gpp*pow(gww,2.0) - 1.0/6.0*xi*pow(gpp,3.0)) - gpp;
    return res;
}

// p = {kf,t,gwomw2,gpomp2,zeta,xi,lambda_v}
void vectorfields_func(int n, double x[], double fvec[], int &iflag, double p[]) {
    double kf = p[0];
    double t = p[1];
    double gwomw2 = p[2];
    double gpomp2 = p[3];
    double zeta = p[4];
    double xi = p[5];
    double lambda_v = p[6];
    double couplings[10] = {0,gwomw2,gpomp2,0,0,0,zeta,xi,lambda_v,0};


    fvec[0] = gww_FE(kf,x[0],x[1],couplings,t);
    fvec[1] = gpp_FE(kf,x[0],x[1],couplings,t);
    return;
}

double vectorfield_2D_NR(double kf, double couplings[10], double t, double eps, double fields[2]) {

    // Initial guesses for the fields
    double gww = 900.0;
    double gpp = -100.0;

    // couplings for vector sector
    double gwomw2 = couplings[1];
    double gpomp2 = couplings[2];
    double zeta = couplings[6];
    double xi = couplings[7];
    double lambda_v = couplings[8];

    double partial_gww_gww, partial_gww_gpp, partial_gpp_gpp, partial_gpp_gww, partial_F_W2, partial_F_B2, partial_F_WB;
    double partial_f_W2, partial_f_B2, partial_f_WB, partial_g_W2, partial_g_B2, partial_g_WB;
    double det, update_gww, update_gpp;

    double y = 2.0*eps;
    int in = 0;
    while (abs(y)>eps) {

        // check if solution satisfies the equation
        y = gww_FE(kf,gww,gpp,couplings,t);

        // calculate Jacobian
        partial_gww_gww = gwomw2*(- 2.0*lambda_v*pow(gpp,2.0) - 0.5*zeta*pow(gww,2.0)) - 1.0;
        partial_gww_gpp = gwomw2*(- 4.0*lambda_v*gww*gpp);
        partial_gpp_gww = gpomp2*(- 4.0*lambda_v*gpp*gww);
        partial_gpp_gpp = gpomp2*(- 2.0*lambda_v*pow(gww,2.0) - 0.5*xi*pow(gpp,2.0)) - 1.0;

        // calculate hessian
        partial_f_W2 = gwomw2*(-zeta*gww);
        partial_f_B2 = gwomw2*(-4.0*lambda_v*gww);
        partial_f_WB = gwomw2*(-4.0*lambda_v*gpp);
        partial_g_W2 = gpomp2*(-4.0*lambda_v*gpp);
        partial_g_B2 = gpomp2*(-xi*gpp);
        partial_g_WB = gpomp2*(-4.0*lambda_v*gww);

        partial_F_W2 = pow(partial_gww_gww,2.0) + gww_FE(kf,gww,gpp,couplings,t)*partial_f_W2 + pow(partial_gpp_gww,2.0) + gpp_FE(kf,gww,gpp,couplings,t)*partial_g_W2;
        partial_F_B2 = pow(partial_gww_gpp,2.0) + gww_FE(kf,gww,gpp,couplings,t)*partial_f_B2 + pow(partial_gpp_gpp,2.0) + gpp_FE(kf,gww,gpp,couplings,t)*partial_g_B2;
        partial_F_WB = partial_gww_gww*partial_gww_gpp + gww_FE(kf,gww,gpp,couplings,t)*partial_f_WB + partial_gpp_gpp*partial_gpp_gww + gpp_FE(kf,gww,gpp,couplings,t)*partial_g_WB;

        if (in>300) {
            // use hessian inverse to calculate newtons update
            det = partial_F_W2*partial_F_B2 - pow(partial_F_WB,2.0);
            update_gww = 1.0/det*(partial_F_B2*(gww_FE(kf,gww,gpp,couplings,t)*partial_gww_gww + gpp_FE(kf,gww,gpp,couplings,t)*partial_gpp_gww) - partial_F_WB*(gww_FE(kf,gww,gpp,couplings,t)*partial_gww_gpp + gpp_FE(kf,gww,gpp,couplings,t)*partial_gpp_gpp));
            update_gpp = 1.0/det*(partial_F_W2*(gww_FE(kf,gww,gpp,couplings,t)*partial_gww_gpp + gpp_FE(kf,gww,gpp,couplings,t)*partial_gpp_gpp) - partial_F_WB*(gww_FE(kf,gww,gpp,couplings,t)*partial_gww_gww + gpp_FE(kf,gww,gpp,couplings,t)*partial_gpp_gww));
        } else {
             // use Jacobian inverse to calculate Newtons update
            det = partial_gww_gww*partial_gpp_gpp - partial_gww_gpp*partial_gpp_gww;
            update_gww = 1.0/det*(partial_gpp_gpp*gww_FE(kf,gww,gpp,couplings,t) - partial_gww_gpp*gpp_FE(kf,gww,gpp,couplings,t));
            update_gpp = 1.0/det*(-partial_gpp_gww*gww_FE(kf,gww,gpp,couplings,t) + partial_gww_gww*gpp_FE(kf,gww,gpp,couplings,t));
        }

        // update the guess
        gww = gww - update_gww;
        gpp = gpp - update_gpp;

        fields[0] = gww;
        fields[1] = gpp;
        //cout << in << "  " << gww << "  " << gpp << "  " << y << endl;
        in = in + 1;
        if (in == 200) {
            gww = 1500;
            gpp = -10.0;
        }
        if (in > 300) {
            cout << "failed to converge on solution for 2Dvector NR at kf: " << kf << endl;
            cout << couplings[0] << "  " << couplings[1] << "  " << couplings[2] << "  " << couplings[3] << "  " << couplings[4] << "  " << couplings[5] << "  " << couplings[6] << "  " << couplings[7] << "  " << couplings[8] << "  " << couplings[9] << endl;
            exit(0);
        }
        //cout << mstarp << "  " << mstarn << endl;
    }
    //cout << gww << "  " <<gpp << "  " << t << endl;
    return 0;

}
//########################################################################
//########################################################################
//########################################################################

// derivative of the sigma field equation (for use in newtons method)
double tools :: gssFE_derivative(double kf, double gss, double gdd, double couplings[10], double t) {
    double res, mstrn,mstrp, kfn,kfp, rn,rp, dsdensp, dsdensn;
    double gsoms2,kappa,lambda,lambda_s;

    gsoms2 = couplings[0];
    kappa = couplings[4];
    lambda = couplings[5];
    lambda_s = couplings[9];

    mstrp = mNuc-gss-0.5*gdd;
    mstrn = mNuc-gss+0.5*gdd;
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // get fermi momenta
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));   // covenient quantities to define
    dsdensp = 1.0/pow(pi,2.0)*(0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrp,2.0))/rp + 3.0*pow(mstrp,2.0)/2.0*log(abs(mstrp/(kf+rp)));
    dsdensn = 1.0/pow(pi,2.0)*(0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrn,2.0))/rn + 3.0*pow(mstrn,2.0)/2.0*log(abs(mstrn/(kf+rn)));
    
    res = gsoms2*(dsdensp + dsdensn - kappa*gss - 0.5*lambda*pow(gss,2.0) - 2.0*lambda_s*pow(gdd,2.0)) - 1.0;
    return res;
}

// derivative of the delta field equation (for use in newtons method) (d/dgdd)
double tools :: gddFE_derivative(double kf, double gss, double gdd, double couplings[10], double t) {
    double res, mstrn,mstrp, kfn,kfp, rn,rp, dsdensp, dsdensn;
    double gdomd2, lambda_s;

    gdomd2 = couplings[3];
    lambda_s = couplings[9];

    mstrp = mNuc-gss-0.5*gdd;
    mstrn = mNuc-gss+0.5*gdd;
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // get fermi momenta
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));   // covenient quantities to define
    dsdensp = -1.0/pow(pi,2.0)*((0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrp,2.0))/rp + 3.0*pow (mstrp,2.0)/2.0*log(abs(mstrp/(kf+rp))));
    dsdensn = 1.0/pow(pi,2.0)*(0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrn,2.0))/rn + 3.0*pow(mstrn,2.0)/2.0*log(abs(mstrn/(kf+rn)));
    
    res = gdomd2*(0.5*dsdensp - 0.5*dsdensn - 2.0*lambda_s*pow(gss,2.0)) - 1.0;
    return res;
}

double tools :: scalarfield_2D_NR(double kf, double couplings[10], double t, double eps, double fields[2]) {
    // get individual nucleon momenta
    double kp = kf*pow(1.0-t,1.0/3.0); 
    double kn = kf*pow(1.0+t,1.0/3.0);

    // Initial guesses for the effective masses
    double mstarp = 900.0;
    double mstarn = 900.0;

    // couplings for scalar sector
    double gsoms2 = couplings[0];
    double gdomd2 = couplings[3];
    double kappa = couplings[4];
    double lambda = couplings[5];
    double lambda_s = couplings[9];
    
    double gss, gdd, ep, en, dsdenspdmstarp, dsdensndmstarn;
    double partial_gss_msp, partial_gss_msn, partial_gdd_msp, partial_gdd_msn;
    double det, update_msp, update_msn;

    double y = 2.0*eps;
    int in = 0;
    while (abs(y)>eps) {
        //cout << mstarn << "  " << mstarp << "  " << y << endl;
        // use effective masses to obtain the fields
        gss = mNuc - 0.5*(mstarn+mstarp);
        gdd = mstarn - mstarp;

        // check if solution satisfies the equation
        y = gss_FE(kf,gss,gdd,couplings,t);

        // calculate Jacobian
        ep = sqrt(pow(kp,2.0) + pow(mstarp,2.0));
        en = sqrt(pow(kn,2.0) + pow(mstarn,2.0));
        dsdenspdmstarp = 0.5/pow(pi,2.0)*(kp*(pow(kp,2.0) + 3.0*pow(mstarp,2.0))*(kp+ep) - 3.0*pow(mstarp,2.0)*(pow(ep,2.0) + kp*ep)*log(abs((kp+ep)/mstarp)))/(pow(ep,2.0) + kp*ep);
        dsdensndmstarn = 0.5/pow(pi,2.0)*(kn*(pow(kn,2.0) + 3.0*pow(mstarn,2.0))*(kn+en) - 3.0*pow(mstarn,2.0)*(pow(en,2.0) + kn*en)*log(abs((kn+en)/mstarn)))/(pow(en,2.0) + kn*en);
        partial_gss_msp = gsoms2*(dsdenspdmstarp + 0.5*kappa*gss + 0.25*lambda*pow(gss,2.0) + lambda_s*pow(gdd,2.0) + 4.0*lambda_s*gss*gdd) + 0.5;
        partial_gss_msn = gsoms2*(dsdensndmstarn + 0.5*kappa*gss + 0.25*lambda*pow(gss,2.0) + lambda_s*pow(gdd,2.0) - 4.0*lambda_s*gss*gdd) + 0.5;
        partial_gdd_msp = gdomd2*(0.5*dsdenspdmstarp + 2.0*lambda_s*pow(gss,2.0) + 2.0*lambda_s*gdd*gss) + 1.0;
        partial_gdd_msn = gdomd2*(-0.5*dsdensndmstarn - 2.0*lambda_s*pow(gss,2.0) + 2.0*lambda_s*gdd*gss) - 1.0;
        
        // use Jacobian inverse to calculare Newtons update
        det = partial_gss_msp*partial_gdd_msn - partial_gss_msn*partial_gdd_msp;
        update_msp = 1.0/det*(partial_gdd_msn*gss_FE(kf,gss,gdd,couplings,t) - partial_gss_msn*gdd_FE(kf,gss,gdd,couplings,t));
        update_msn = 1.0/det*(-partial_gdd_msp*gss_FE(kf,gss,gdd,couplings,t) + partial_gss_msp*gdd_FE(kf,gss,gdd,couplings,t));

        // update the guess
        mstarp = mstarp - update_msp;
        mstarn = mstarn - update_msn;

        if (mstarn < 0) {
            mstarn = (mstarn+mNuc)/2.0;
        }

        if (mstarp < 0) {
            mstarp = (mstarp+mNuc)/2.0;
        }

        fields[0] = gss;
        fields[1] = gdd;

        in = in + 1;
        //cout << gss << "  " << gdd << "  " << y << endl;
        if (in == 50) {
            mstarp = 100;
            mstarn = 100;
        }
        if (in > 100) {
            cout << "failed to converge on solution for 2D NR at: " << kf << endl;
            cout << couplings[0] << "  " << couplings[1] << "  " << couplings[2] << "  " << couplings[3] << "  " << couplings[4] << "  " << couplings[5] << "  " << couplings[6] << "  " << couplings[7] << "  " << couplings[8] << "  " << couplings[9] << endl;
            exit(0);
        }
        //cout << mstarp << "  " << mstarn << endl;
    }
    //cout << mstarp << "  " << mstarn << "  " << t << endl;
    return 0;

}

// bisection for obtaining the sigma field
double tools :: FSUgssnewton(double dens, double gdd, double couplings[10], double t, double eps) {
    double ymin,error,y,min,max,sol,k;
    int ib, in ,MAXIT;
    
    error = eps;           // set min error
    y = error*2;           // initialize the bisection
    min = 1e-5; max = mNuc+1e-5;; // set bounds on solution
    sol = (min+max)/2.0;    // solution starts as midpoint
    k = pow(3.0/2.0*pow(pi,2)*dens,1.0/3.0);    // convert density to fermi momentum
    
    ib=0;   // counter for bisection
    in = 0; // counter for newtons
    MAXIT = 200;    // maximum number of iterations for newtons allowed
    double bis;
    while (abs(y)>error) {
        y = gss_FE(k,sol,gdd,couplings,t);

        ymin = gss_FE(k,min,gdd,couplings,t);
        if (y*ymin<0) {
            max = sol;
        } else {
            min = sol;
        }
        bis = (min+max)/2.0;
        ib = ib +1;

        sol = bis;
        //cout << ib << "  " << in << "  " << ir << "  " << sol << "  " << y << endl;
        if (ib > MAXIT) {
            cout << "bisection failed gss: " << endl;
            cout << "bisect: " << ib << "  " << in << "  " << y << "  " << sol << " for t: " << t << endl;
            exit(0);
        }
    }
    return sol;
}

// newtons method for obtaining the delta field
double tools :: FSUgddnewton(double dens, double couplings[10], double t, double eps) {
    double dy, newt, ymin,error,y,min,max,sol,k,gss,ymax;
    int ib, in;
    int MAXIT = 500;
    
    error = eps;           // set min error
    y = error*2;           // initialize the bisection
    min = -mNuc+2.0; max = mNuc-1.0; // set bounds on solution
    sol = (min+max)/2.0;    // solution starts as midpoint
    k = pow(3.0/2.0*pow(pi,2)*dens,1.0/3.0);    // convert density to fermi momentum
    
    ib=0;   // counter for bisection
    in = MAXIT+1; // counter for newtons

    if (couplings[3] == 0) {
        return 0;
    }

    while (abs(y)>error) {
        gss = FSUgssnewton(dens,sol,couplings,t,eps);
        y = gdd_FE(k,gss,sol,couplings,t);
        //cout << in << "  " << ib << endl;
        dy = gddFE_derivative(k,gss,sol,couplings,t);
        newt = sol - y/dy;
        //cout << sol << "  " << y << endl;
        // if newtons method yields a value outside the solution bounds or its reached its maximum iterations allowed then use bisection
        if (newt > max || newt<min || in > MAXIT) {
            gss = FSUgssnewton(dens,min,couplings,t,eps);
            ymin = gdd_FE(k,gss,min,couplings,t);
            ymax = gdd_FE(k,gss,max,couplings,t);
            if (y*ymin<0) {
                max = sol;
            } else {
                min = sol;
            }
            sol = (ymax*min - ymin*max)/(ymax- ymin);    // new midpoint
            ib = ib+1;              // bisection count
        } else {
            sol = newt;             // newtons solution
            in = in+1;              // newtons count
        }
        
        if (ib > MAXIT && y>1.0e-4) {
            cout << "bisection failed gdd: " << endl;
            cout << "bisect: " << ib << "  " << y << "  " << sol << " for t = " << t << endl;
            exit(0);
        } else if (ib > MAXIT) {
            return sol;
        }
    }
    //cout << y << endl;
    return sol;
}

// newton and bisection method to get omega field
double tools :: FSUgwwnewton(double couplings[10], double dens, double t) {
    double dy, gpp, newt, ymin, error,y,min,max,x;
    int MAXIT, in,ib;
    double gpomp2,lambda_v,gwomw2,zeta;

    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    zeta = couplings[6];
    lambda_v = couplings[8];
    
    error = 1e-11;           // set min error
    y = error*2;            // initialize bisection
    min = 1e-5; max = 2000; // set bounds on solution
    MAXIT = 500;   // maximum number of iterations for newtons method
    in = 0; // count for newtons
    ib = 0; // count for bisection
    x = 400.0;   // solution starts at midpoint
    
    while (abs(y)>error) {
        // newtons method
        gpp = -0.5*gpomp2*dens*t/(1.0+2.0*lambda_v*gpomp2*pow(x,2.0));
        y = x + gwomw2*(2.0*lambda_v*pow(gpp,2.0)*x + 1.0/6.0*zeta*pow(x,3.0) - dens);
        dy = 1.0 + gwomw2*(2.0*lambda_v*pow(gpp,2.0) + 0.5*zeta*pow(x,2.0));
        newt = x - y/dy;
        
        // if newtons method yields a value outside the solution bounds or its reached its maximum iterations allowed then use bisection
        if (newt > max || newt<min || in > MAXIT) {
            gpp = -0.5*gpomp2*dens*t/(1.0+2.0*lambda_v*gpomp2*pow(min,2.0));
            ymin = min + gwomw2*(2.0*lambda_v*pow(gpp,2.0)*min + 1.0/6.0*zeta*pow(min,3.0) - dens);
            if (y*ymin<0) {
                max = x;
            } else {
                min = x;
            }
            x = (min+max)/2.0;  // new midpoint
            ib = ib+1;  // bisection count
        } else {
            x = newt;   // newtons solution
            in = in+1;  // newton count
        }

        // if bisection also fails then exit
        if (ib > MAXIT) {
            cout << "bisection failed gww: " << endl;
            cout << "bisect: " << ib << "  " << y << "  " << x << "  " << gpp << "  " << dens << "  " << t << endl;
            exit(0);
        }
    }
    return x;
}

// Calculate (g_rho/m_rho)^2 Analytically
double tools :: get_gpomp2(double kf, double J, double L, double gss, double gww, double gsoms2, double gwomw2, double gdomd2, double kappa, double lambda, double zeta, double lambda_s) {

    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    double mstar = mP - gss;
    double en = sqrt(pow(kf,2) + pow(mstar,2));
    double integral = tool.common_integral(kf,mstar);
    double alpha_s = 1.0 + gsoms2*(2.0/pow(pi,2.0)*integral + kappa*gss + 0.5*lambda*pow(gss,2.0));
    double dgssdp = gsoms2*(mstar/en)*pow(alpha_s,-1.0);  // good
    double dmstardp = -dgssdp;   // good
    double dendp = 1.0/en*(0.5*pow(pi,2.0)/kf + mstar*dmstardp); // good
    double dgdddt = -pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(mstar/en)*pow(1.0 + gdomd2*0.5/pow(pi,2.0)*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0),-1.0);
    double alpha_d = 1.0 + gdomd2*1.0/(2.0*pow(pi,2.0))*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0);
    double dIdp = 0.5*pow(kf*pi,2.0)/pow(en,3.0) - 3.0*dmstardp*mstar*(-1.0/3.0*pow(kf/en,3.0) - kf/en + log(abs((kf+en)/mstar)));
    double dgdddpdt =  - 0.5*gdomd2*(mstar/en)*pow(alpha_d,-1) - pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(dmstardp/en - mstar/pow(en,2.0)*dendp)*pow(alpha_d,-1.0)
                + pow(kf,3.0)/(3.0*pow(pi,2.0))*pow(gdomd2,2.0)*(mstar/en)*(1.0/(2.0*pow(pi,2.0))*dIdp + 4.0*lambda_s*gss*dgssdp)*pow(alpha_d,-2.0);
    double phi = pow(pi,2.0)/(6.0*kf*en) - pow(kf,2.0)/(6.0*pow(en,2.0))*dendp + 0.25/en*(dmstardp - mstar/en*dendp)*dgdddt + 0.25*mstar/en*dgdddpdt;
    double chi = 12.0*pow(pi,2.0)/pow(kf,3.0)*J - 2.0*pow(pi,2.0)/kf*1.0/en + pow(mstar/en,2.0)*gdomd2*pow(alpha_d,-1.0);
    double gpomp2 = 2.0*gwomw2*pow(1.0+0.5*zeta*gwomw2*pow(gww,2.0),-1.0)*pow(chi,2.0)/(2.0*gwomw2*pow(1.0+0.5*zeta*gwomw2*pow(gww,2.0),-1.0)*chi - (4.0*phi - 4.0*L/(3.0*p0) + 0.5*chi)*3.0*pow(pi,2.0)/pow(kf,3.0)*gww);
    return gpomp2;
}

// need to change
// get gpomp2 along with the two other isovector parameters (gpomp2, lambda) via bisection on Ksym
double tools :: get_gdomd2(double kf, double J, double L, double Ksym, double gss, double gww, double gsoms2, double gwomw2, double kappa, double lambda, double zeta, double lambda_s, int sol) {
    double ymin, error,y,min,max,x;
    int MAXIT,ib;
    double gpomp2, lambda_v;
    
    error = 1e-6;           // set min error
    y = error*2;            // initialize bisection
    min = 1.0e-15; max = 0.01; // set bounds on solution
    MAXIT = 100;   // maximum number of iterations for newtons method
    ib = 0; // count for bisection
    x = min;   // solution starts at midpoint

    // fill couplings array
    double couplings[10];
    couplings[0] = gsoms2; couplings[1] = gwomw2; couplings[4] = kappa; couplings[5] = lambda; couplings[6] = zeta; couplings[9] = lambda_s; 
    
    double** gd_array;
    dmi.create(gd_array,5000,2);
    for (int i=0; i<5000; ++i) {
        gpomp2 = get_gpomp2(kf,J,L,gss,gww,gsoms2,gwomw2,x,kappa,lambda,zeta,lambda_s);
        lambda_v = get_lambda_v(kf,J,gss,gww,x,gpomp2,lambda_s);
        couplings[2] = gpomp2;
        couplings[3] = x;
        couplings[8] = lambda_v;
        y = bulk.get_Ksym(kf,gss,gww,couplings) - Ksym;
        gd_array[i][0] = x;
        gd_array[i][1] = y;
        x = x + 5e-6;
    }
    dmi.print(gd_array,5000,2,true,"gd_array.txt");

    int i = 0;
    double sgn = 1.0;
    if (sol == 1) {
        while (sgn > 0) {
            sgn = gd_array[4999-i][1]*gd_array[4999-i-1][1];
            min = gd_array[4999-i-1][0];
            max = gd_array[4999-i][0];
            i = i+1;
            if (i == 4999) {
                cout << " no zero exists gd2 for J: " << J << endl;
                dmi.print(gd_array,5000,2,true,"gd_array.txt");
                dmi.cleanup(gd_array,5000);
                return -1;
            }
        }
    } else if(sol == -1) {
        while (sgn > 0) {
            sgn = gd_array[i][1]*gd_array[i+1][1];
            min = gd_array[i][0];
            max = gd_array[i+1][0];
            i = i+1;
            if (i == 4999) {
                cout << " no zero exists gd2 for J: " << J << endl;
                dmi.print(gd_array,5000,2,true,"gd_array.txt");
                dmi.cleanup(gd_array,5000);
                exit(0);
            }
        }
    } else {
        cout << "invalid solution selection: 1 for upper and -1 for lower" << endl;
        dmi.cleanup(gd_array,5000);
        exit(0);
    }
    //dmi.print(gd_array,5000,2,true,"gd_array.txt");
    dmi.cleanup(gd_array,5000);

    x = (min+max/2.0);
    //cout << min << "  " << max << endl;
    while (abs(y)>error) {
        gpomp2 = get_gpomp2(kf,J,L,gss,gww,gsoms2,gwomw2,x,kappa,lambda,zeta,lambda_s);
        lambda_v = get_lambda_v(kf,J,gss,gww,x,gpomp2,lambda_s);
        couplings[2] = gpomp2;
        couplings[3] = x;
        couplings[8] = lambda_v;
        y = bulk.get_Ksym(kf,gss,gww,couplings) - Ksym;

        gpomp2 = get_gpomp2(kf,J,L,gss,gww,gsoms2,gwomw2,min,kappa,lambda,zeta,lambda_s);
        lambda_v = get_lambda_v(kf,J,gss,gww,min,gpomp2,lambda_s);
        couplings[2] = gpomp2;
        couplings[3] = min;
        couplings[8] = lambda_v;
        ymin = bulk.get_Ksym(kf,gss,gww,couplings) - Ksym;

        if (y*ymin>0) {
            min = x;
        } else {
            max = x;
        }
        //x = (min*ymax - max*ymin)/(ymax-ymin);  // new midpoint
        x = (min+max)/2.0;
        ib = ib+1;  // bisection count

        // if bisection also fails then exit
        if (ib > MAXIT) {
            cout << "bisection failed gdomd2: " << endl;
            cout << "bisect: " << ib << "  " << y << "  " << x << endl;
            exit(0);
        }
        
    }
    return x;
}

// need to change
double tools :: get_lambda_v(double kf, double J, double gss, double gww, double gdomd2, double gpomp2, double lambda_s) {
    double mstar = mNuc - gss;
    double en = sqrt(pow(kf,2.0) + pow(mstar,2.0));

    double integral = tool.common_integral(kf,mstar);
    double alpha_d = 1.0 + gdomd2*1.0/(2.0*pow(pi,2.0))*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0);
    double lambda_v = 0.5*(pow(12.0*pow(pi,2.0)/pow(kf,3.0)*J - 2.0*pow(pi,2.0)/(kf*en) + pow(mstar/en,2.0)*gdomd2*pow(alpha_d,-1.0),-1.0)*gpomp2 - 1.0)/(pow(gww,2.0)*gpomp2);
    return lambda_v;
}

void tools :: convert_to_inf_couplings(double fin_couplings[16], double inf_couplings[10]) {
    inf_couplings[0] = fin_couplings[0]/pow(fin_couplings[12],2.0);
    inf_couplings[1] = fin_couplings[1]/pow(fin_couplings[13],2.0);
    inf_couplings[2] = fin_couplings[2]/pow(fin_couplings[14],2.0);
    inf_couplings[3] = fin_couplings[3]/pow(fin_couplings[15],2.0);
    inf_couplings[4] = fin_couplings[4];
    inf_couplings[5] = fin_couplings[5];
    inf_couplings[6] = fin_couplings[6];
    inf_couplings[7] = fin_couplings[7];
    inf_couplings[8] = fin_couplings[8];
    inf_couplings[9] = fin_couplings[9];
}

double chneutral(double kf, double couplings[10], double t, double fields_v[2]) {
    double kp, kn, gss, gpp, gdd, munmmup, y, np,ne,nmu, mstarn,mstarp;
    
    kp = kf*pow(1.0-t,1.0/3.0); kn = kf*pow(1.0+t,1.0/3.0); // get nucleon fermi momenta
    np = pow(kp,3.0)/(3.0*pow(pi,2));   // proton density

    double fields_s[2];
    tool.scalarfield_2D_NR(kf,couplings,t,1e-10,fields_s);
    gdd = fields_s[1];
    gss = fields_s[0];
    mstarp = mNuc-gss-0.5*gdd;   mstarn = mNuc-gss+0.5*gdd; // effective mass
    double fvec[2]; int lwa = (2*(3*2+13))/2;  double wa[lwa];
    double params[7] = {kf,t,couplings[1],couplings[2],couplings[6],couplings[7],couplings[8]};
    hybrd1(vectorfields_func,2,fields_v,fvec,1e-6,wa,lwa,params);
    //vectorfield_2D_NR(kf,couplings,t,1e-7,fields);
    gpp = fields_v[1];
    munmmup = sqrt(pow(kn,2) + pow(mstarn,2)) - sqrt(pow(kp,2) + pow(mstarp,2)) - gpp;    // mun - mup = mue
    ne = eos.qknb(munmmup,mE,2.0);
    nmu = eos.qknb(munmmup,mMU,2.0);
    y = np - ne - nmu;
    return y;
}

// perform newtons method and bisection to get t for a given fermi momentum kf
double tools :: get_t_betaeq(double kf, double couplings[10], double t_avg, double fields_v[2]) {
    double min, max;
    min = t_avg - 0.040;
    max = t_avg + 0.040;
    if (max>1.0) {
        max = 0.99999;
    }
    if (min<0.0) {
        min = 1e-8;
    }
    return nmi.zbrent(chneutral,min,max,1e-7,kf,couplings,fields_v);
}

// need to change
double get_J_given_Jtilde(double Jtilde, double kf, double BA, double mstar, double gsoms2, double gwomw2, double kappa, double lambda, double zeta, double xi, double lambda_s, double L, double Ksym, int gd_sol_type, bool delta_coupling) {
    double subkf,subp0,gdomd2,gpomp2,lambda_v,gss,gww,x,ymin;
    double inf_couplings[10];
    double **arrSNM1, **arrSNM2;
    int ncols = 6; int npoints = 1000; int MAXIT,ib;
    double error,y,min,max;
    subkf = 1.15*197.32698;
    subp0 = 2.0/(3.0*pow(pi,2.0))*pow(subkf,3.0);
    
    error = 1e-4;           // set min error
    y = error*2;            // initialize bisection
    min = 25.0; max = 45.0; // set bounds on solution
    MAXIT = 500;   // maximum number of iterations for newtons method
    ib = 0; // count for bisection
    x = (min+max)/2.0;   // solution starts at midpoint

    // bisect in stable region
    gss = mNuc - mstar; // given mstar get gss
    gww = mNuc + BA - sqrt(pow(kf,2) + pow(mstar,2)); 
    if (delta_coupling == true) {
            gdomd2 = tool.get_gdomd2(kf,x,L,Ksym,gss,gww,gsoms2,gwomw2,kappa,lambda,zeta,lambda_s,gd_sol_type);
    } else {
            gdomd2 = 0.0;
    }
    while (gdomd2 < 0) {
        min = min + 1;
        gdomd2 = tool.get_gdomd2(kf,min,L,Ksym,gss,gww,gsoms2,gwomw2,kappa,lambda,zeta,lambda_s,gd_sol_type);
        if (min > 45) {
            cout << " no gdom2 found" << endl;
            exit(0);
        }
    }

    while (abs(y)>error) {
        // newtons method
        gss = mNuc - mstar; // given mstar get gss
        gww = mNuc + BA - sqrt(pow(kf,2) + pow(mstar,2));    // get gww at saturation
        if (delta_coupling == true) {
            gdomd2 = tool.get_gdomd2(kf,x,L,Ksym,gss,gww,gsoms2,gwomw2,kappa,lambda,zeta,lambda_s,gd_sol_type);
        } else {
            gdomd2 = 0.0;
        }
        gpomp2 = tool.get_gpomp2(kf,x,L,gss,gww,gsoms2,gwomw2,gdomd2,kappa,lambda,zeta,lambda_s);
        lambda_v = tool.get_lambda_v(kf,x,gss,gww,gdomd2,gpomp2,lambda_s);
        inf_couplings[0] = gsoms2; inf_couplings[1] = gwomw2; inf_couplings[2] = gpomp2;
        inf_couplings[3] = gdomd2; inf_couplings[4] = kappa; inf_couplings[5] = lambda; inf_couplings[6] = zeta;
        inf_couplings[7] = xi; inf_couplings[8] = lambda_v; inf_couplings[9] = lambda_s;
        eos.get_SNMEOS(inf_couplings,arrSNM1,npoints);
        y = arrSNM1[dmi.findvalue(arrSNM1,npoints,ncols,subp0,0,0.01)][4] - Jtilde;
        //cout << "y: " << x << "  " << y << "  " << gdomd2*pow(980,2.0) << endl;
        dmi.cleanup(arrSNM1,npoints);  // clear the SNM array
            
            
        gss = mNuc - mstar; // given mstar get gss
        gww = mNuc + BA - sqrt(pow(kf,2) + pow(mstar,2)); 
        if (delta_coupling == true) {
            gdomd2 = tool.get_gdomd2(kf,min,L,Ksym,gss,gww,gsoms2,gwomw2,kappa,lambda,zeta,lambda_s,gd_sol_type);
        } else {
            gdomd2 = 0.0;
        }
        gpomp2 = tool.get_gpomp2(kf,min,L,gss,gww,gsoms2,gwomw2,gdomd2,kappa,lambda,zeta,lambda_s);
        lambda_v = tool.get_lambda_v(kf,min,gss,gww,gdomd2,gpomp2,lambda_s);
        inf_couplings[0] = gsoms2; inf_couplings[1] = gwomw2; inf_couplings[2] = gpomp2;
        inf_couplings[3] = gdomd2; inf_couplings[4] = kappa; inf_couplings[5] = lambda; inf_couplings[6] = zeta;
        inf_couplings[7] = xi; inf_couplings[8] = lambda_v; inf_couplings[9] = lambda_s;
        eos.get_SNMEOS(inf_couplings,arrSNM2,npoints);
        ymin = arrSNM2[dmi.findvalue(arrSNM2,npoints,ncols,subp0,0,0.01)][4] - Jtilde;
        //cout << "min: " << min << "  " << ymin << "  " << gdomd2*pow(980,2.0) << endl;
        dmi.cleanup(arrSNM2,npoints);  // clear the SNM array
        
        if (y*ymin<0) {
            max = x;
        } else {
            min = x;
        }
        x = (min+max)/2.0;  // new midpoint
        ib = ib+1;  // bisection count
        // if bisection also fails then exit
        if (ib > MAXIT) {
            cout << "bisection failed Jtilde: " << endl;
            cout << "bisect: " << ib << "  " << y << "  " << x << "  " << endl;
            exit(0);
        }
    }
    return x;
}

//###########################################################
// EOS class
// ##########################################################

// Get the energy density for the FSU+delta model for infinite nuclear matter
double equationofstate :: get_en(double kf, double t, double gss, double gww, double gpp, double gdd, double couplings[10]) {
    double res,mss2,mww2,mpp2,mdd2,mstarp,mstarn,kfn,kfp;
    double gsoms2,gwomw2,gpomp2,gdomd2,kappa,lambda,zeta,xi,lambda_v,lambda_s;
    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];
    xi = couplings[7];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    mss2 = pow(gss,2.0)/gsoms2;  // (m_sigma*sigma)^2
    mww2 = pow(gww,2.0)/gwomw2;  // (m_omgea*omega)^2
    mpp2 = pow(gpp,2.0)/gpomp2;  // (m_rho*rho)^2
    mdd2 = pow(gdd,2.0)/gdomd2;  // (m_delta*delta)^2
    if (gdomd2 == 0) {
        mdd2 = 0;
    }

    mstarp = mP-gss-0.5*gdd; mstarn = mN-gss+0.5*gdd;   // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);     // define the fermi momentum for the proton and neutron given t
    
    // add free contributions
    res = 0.5*mss2 + 0.5*mww2 + 0.5*mpp2 + 0.5*mdd2 + 1.0/6.0*kappa*pow(gss,3.0) + 1.0/24.0*lambda*pow(gss,4.0) + 1.0/8.0*zeta*pow(gww,4.0) + 1.0/8.0*xi*pow(gpp,4.0)
        + 3.0*lambda_v*pow(gpp,2.0)*pow(gww,2.0) + lambda_s*pow(gdd,2.0)*pow(gss,2.0) + tool.endensity_integral(kfp,mstarp) + tool.endensity_integral(kfn,mstarn);
    return res;
}


// Get the pressure for the FSU model
double equationofstate :: get_pr(double kf, double t, double gss, double gww, double gpp, double gdd, double couplings[10]) {
    double res,mss2,mww2,mpp2,mdd2,mstarp,mstarn,kfn,kfp;
    double gsoms2,gwomw2,gpomp2,gdomd2,kappa,lambda,zeta,xi,lambda_v,lambda_s;
    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];
    xi = couplings[7];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    mss2 = pow(gss,2.0)/gsoms2;  // (m_sigma*sigma)^2
    mww2 = pow(gww,2.0)/gwomw2;  // (m_omgea*omega)^2
    mpp2 = pow(gpp,2.0)/gpomp2;  // (m_rho*rho)^2
    mdd2 = pow(gdd,2.0)/gdomd2;  // (m_delta*delta)^2
    if (gdomd2 == 0) {
        mdd2 = 0;
    }

    mstarp = mP-gss-0.5*gdd; mstarn = mN-gss+0.5*gdd;   // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);     // define the fermi momentum for the proton and neutron given t
    
    // add free contributions
    res = -0.5*mss2 - 0.5*mdd2 + 0.5*mww2 + 0.5*mpp2 - 1.0/6.0*kappa*pow(gss,3.0) - 1.0/24.0*lambda*pow(gss,4.0) + 1.0/24.0*zeta*pow(gww,4.0) + 1.0/24.0*xi*pow(gpp,4.0)
            + lambda_v*pow(gpp,2.0)*pow(gww,2.0) - lambda_s*pow(gdd,2.0)*pow(gss,2.0) + tool.pr_integral(kfn,mstarn) + tool.pr_integral(kfp,mstarp);
    return res;
}


double equationofstate :: qken(double cp, double mass, double deg) {
    double z, en;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        en = deg*pow(cp,4)/(8.0*pow(pi,2))*( sqrt(1.0-pow(z,2))*(1.0-0.5*pow(z,2)) - 0.5*pow(z,4)*log((1.0+sqrt(1.0-pow(z,2)))/z) );
        return en;
    }
}


double equationofstate :: qkpr(double cp, double mass, double deg) {
    double z, pr;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        pr = deg*pow(cp,4)/(24.0*pow(pi,2))*( sqrt(1.0-pow(z,2))*(1.0-5.0/2.0*pow(z,2)) + 1.5*pow(z,4)*log((1.0+sqrt(1.0-pow(z,2)))/z) );
        return pr;
    }
}

double equationofstate :: qknb(double cp, double mass, double deg) {
    double z, nb;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        nb = deg*pow(cp,3)/(6.0*pow(pi,2))*pow((1.0-pow(z,2)),3.0/2.0);
        return nb;
    }
}

// get the EOS for a given number of points and set of couplings
// will output an array of the form [nb fm-3, en MeV/fm3, pr MeV/fm3, dpde, mub MeV, dP/drho, Keff MeV]
int equationofstate :: get_SNMEOS(double couplings[10], double** &eos, int npoints) {
    double k,en,dens,gww,gss,gpp,gdd,mstar,t,p0f,ssize;
    k = 40.0;   // initial fermi momentum
    p0f = 1.5; // N TIMES SATURATION DENSITY

    // step size that increases with increasing density
    double conv_mev4 = 1.0/pow(197.32698,3); // mev4 to mev/fm3
    double kf = pow(3.0*pow(pi,2.0)/2.0*p0f*0.15/conv_mev4,1.0/3.0);
    ssize = (kf-k)/npoints;

    // create EOS array
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[6];
    }
    
    // return the EOS for a specified number of points
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density at given fermi momentum
        t = 0;
        gdd = 0.0;
        gss = tool.FSUgssnewton(dens,gdd,couplings,t,1e-7);     // get the sigma field
        gww = tool.FSUgwwnewton(couplings,dens,t);  // get omega field
        gpp = 0.0;
        mstar = mNuc - gss;

        // get the EOS
        en = get_en(k,t,gss,gww,gpp,gdd,couplings);
        eos[i][0] = dens;
        eos[i][1] = en;
        eos[i][2] = mstar/mNuc;
        eos[i][3] = en/dens - mNuc;
        eos[i][4] = bulk.get_J(k,gss,gww,couplings);
        eos[i][5] = bulk.get_L(k,gss,gww,couplings);

        k = k + ssize;
    }
    return 0;
}

int equationofstate :: get_PNMEOS(double couplings[10], double** &eos, int npoints) {
    double k,en,dens,gww,gss,gpp,gdd,t,p0f,ssize;
    k = 40.0;   // initial fermi momentum
    p0f = 2.0; // N TIMES SATURATION DENSITY

    // step size that increases with increasing density
    double conv_mev4 = 1.0/pow(197.32698,3); // mev4 to mev/fm3
    double kf = pow(3.0*pow(pi,2.0)/2.0*p0f*0.15/conv_mev4,1.0/3.0);
    ssize = (kf-k)/npoints;

    // create EOS array
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[2];
    }
    
    double fields[2];
    // return the EOS for a specified number of points
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density at given fermi momentum
        t = 1.0;
        tool.scalarfield_2D_NR(k,couplings,t,1e-7,fields);
        gdd = fields[1];
        gss = fields[0];
        vectorfield_2D_NR(k,couplings,t,1e-7,fields);
        gpp = fields[1];
        gww = fields[0];

        // get the EOS
        en = get_en(k,t,gss,gww,gpp,gdd,couplings);
        eos[i][0] = dens*conv_mev4;
        eos[i][1] = en/dens - mNuc;
        cout << gdd << "  " << gss << "  " << gpp << "  " << gww << "  " << eos[i][1] << endl;

        k = k + ssize;
    }
    return 0;
}

int equationofstate :: get_SymmetryEnergy(double couplings[10], double** &eos, int npoints) {
    double k,en_PNM,en_SNM,dens,gww,gss,gpp,gdd,t,p0f,ssize;
    k = 40.0;   // initial fermi momentum
    p0f = 2.0; // N TIMES SATURATION DENSITY

    // step size that increases with increasing density
    double conv_mev4 = 1.0/pow(197.32698,3); // mev4 to mev/fm3
    double kf = pow(3.0*pow(pi,2.0)/2.0*p0f*0.15/conv_mev4,1.0/3.0);
    ssize = (kf-k)/npoints;

    // create EOS array
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[3];
    }
    
    double fields[2];
    // return the EOS for a specified number of points
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density at given fermi momentum
        t = 1.0;
        tool.scalarfield_2D_NR(k,couplings,t,1e-7,fields);
        gdd = fields[1];
        gss = fields[0];
        vectorfield_2D_NR(k,couplings,t,1e-7,fields);
        gpp = fields[1];
        gww = fields[0];

        // get the EOS
        en_PNM = get_en(k,t,gss,gww,gpp,gdd,couplings)/dens - mNuc;

        t = 0.0;
        tool.scalarfield_2D_NR(k,couplings,t,1e-7,fields);
        gdd = fields[1];
        gss = fields[0];
        vectorfield_2D_NR(k,couplings,t,1e-7,fields);
        gpp = fields[1];
        gww = fields[0];

        en_SNM = get_en(k,t,gss,gww,gpp,gdd,couplings)/dens - mNuc;

        eos[i][0] = dens*conv_mev4;
        eos[i][1] = en_PNM-en_SNM;
        eos[i][2] = bulk.get_J(k,gss,gww,couplings);
        k = k + ssize;
    }
    return 0;
}


// get the EOS for a given number of points and set of couplings
// coupling array of the form [ (gs/ms)^2 , (gw/mw)^2 , (gp/mp)^2, (gd/md)^2, kappa , lambda , zeta, xi, lambda_v , lambda_s ]
// will output an array of the form [nb fm-3, en MeV/fm3, pr MeV/fm3, dpde, mub MeV, dP/drho, Keff MeV]
int equationofstate :: get_EOS_NSM(double couplings[10], double** &eos, int npoints, bool print, bool unstable) {
    double k,en,pr,dens,mue,gww,gss,gpp,gdd,mstarp,mstarn,kp,kn,mun,mup,t,check,conv_mev4,p0f,ssize;
    double dydx1,dydx2;
    double Yp,Ze, Yp_Urca;
    conv_mev4 = 1.0/pow(197.32698,3); // mev4 to mev/fm3
    k = 40.0;   // initial fermi momentum
    p0f = 6.5; // N TIMES SATURATION DENSITY
    double kf = pow(3.0*pow(pi,2.0)/2.0*p0f*0.15/conv_mev4,1.0/3.0);
    // step size that increases with increasing density
    ssize = (kf-k)/npoints;
    double fields_v[2];
    double fields_s[2];

    // create EOS array
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[9];
    }
    
    // return the EOS for a specified number of points
    double t_avg = 1.0; // initial guess for t_avg
    fields_v[0] = 10.0; fields_v[1] = -10.0; // inital guesses for vector fields
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density at given fermi momentum
        t = tool.get_t_betaeq(k,couplings,t_avg,fields_v); // get t from charge neutrality and beta equil
        t_avg = t;
        kp = k*pow(1.0-t,1.0/3.0); kn = k*pow(1.0+t,1.0/3.0);   // get proton and neutron fermi momenta
        tool.scalarfield_2D_NR(k,couplings,t,1e-7,fields_s);
        gdd = fields_s[1];
        gss = fields_s[0];
        double fvec[2]; int lwa = (2*(3*2+13))/2;  double wa[lwa];
        double params[7] = {k,t,couplings[1],couplings[2],couplings[6],couplings[7],couplings[8]};
        hybrd1(vectorfields_func,2,fields_v,fvec,1e-6,wa,lwa,params);

        //vectorfield_2D_NR(k,couplings,t,1e-7,fields);
        gpp = fields_v[1];
        gww = fields_v[0];
        
        //gww = tool.FSUgwwnewton(couplings,dens,t);
        //gpp = -0.5*gpomp2*dens*t/(1.0+2.0*lambda_v*gpomp2*pow(gww,2.0));  // get rho field
        mstarp = mP - gss - 0.5*gdd; mstarn = mN - gss + 0.5*gdd; // get efffecitve masses
        //cout << setprecision(10) << dens*conv_mev4 << "  " << t << "  " << gdd << "  " << gss << "  " << gww << "  " << gpp << endl;
        //cout << "eff rho: " << sqrt(pow(763.0,2.0) + 2.0*couplings[2]*pow(763,2.0)*couplings[8]*pow(gww,2.0)) << endl;
        //double eff_w_mass = 1.0/couplings[1] + 1.0/6.0*couplings[6]*pow(gww,2.0);
        //cout << "check t: " << chneutral(k,couplings,t) << endl;
        //cout << "---------------------------------" << endl;
        
        // get the chemical potentials
        mup = sqrt(pow(kp,2) + pow(mstarp,2)) + gww + 0.5*gpp;  
        mun = sqrt(pow(kn,2) + pow(mstarn,2)) + gww - 0.5*gpp;
        mue = mun - mup;

        Ze = pow(pow(mue,2.0) - pow(mE,2.0),3.0/2.0)/pow(kp,3.0);
        Yp = pow(kp,3.0)/(pow(kn,3.0) + pow(kp,3.0));
        Yp_Urca = 1.0/(1.0 + pow(1.0 + pow(Ze,1.0/3.0),3.0));

        // get the EOS
        en = get_en(k,t,gss,gww,gpp,gdd,couplings) + qken(mue,mE,2.0) + qken(mue,mMU,2.0);
        pr = get_pr(k,t,gss,gww,gpp,gdd,couplings) + qkpr(mue,mE,2.0) + qkpr(mue,mMU,2.0);
        eos[i][0] = dens*conv_mev4;
        eos[i][1] = conv_mev4*en;
        eos[i][2] = conv_mev4*pr;
        eos[i][4] = mun;
        eos[i][5] = 0.5*sqrt(pow(kp,2.0)+pow(mstarp,2.0))*pow(kp/k,2.0)*pow(1.0-t,1.0/3.0) + 0.5*sqrt(pow(k,2.0)+pow(mstarp,2.0))*pow(kn/k,2.0)*pow(1.0+t,1.0/3.0) + gww;   // de/drho
        eos[i][6] = tool.effK2(k,couplings,gss,gww,gpp,gdd,t);
        eos[i][7] = Yp;
        eos[i][8] = Yp_Urca;
        //cout << gss << "  " << gww << "  " << gpp << "  " << gdd << "  " << t << endl;
        //double ntotal = 0.5*dens*(1.0+t) + 0.5*dens*(1.0-t) + qknb(mue,mE,2.0) +  qknb(mue,mMU,2.0);
        //cout << 0.5*dens*(1.0+t)/ntotal << "  " << 0.5*dens*(1.0-t)/ntotal << "  " << qknb(mue,mE,2.0)/ntotal << "  " << qknb(mue,mMU,2.0)/ntotal << endl;
    
        // check for thermodynamic stability
        check = mun*dens - en - pr;
        //cout << "ch: " << -qknb(mue,mE,2.0) - qknb(mue,mMU,2.0) + pow(kp,3.0)/(3.0*pow(pi,2));
        //cout << "  gss: " << gss << "  gww: " << gww << "  gpp: " << gpp << "  gdd: " << gdd << "  t: " << t << endl;
        if (abs(check)>10000.0) {
            cout << "consistency fail for k: " << k << " with value: " << check << endl;
            cout << pr*conv_mev4 << "  " << (mun*dens - en)*conv_mev4 << endl;
            for (int j=0; j<1000; ++j) {
                double w = j*1.0/(1000-1);
                tool.scalarfield_2D_NR(k,couplings,w,1e-10,fields_s);
                gss = fields_s[0];
                vectorfield_2D_NR(k,couplings,w,1e-7,fields_v);
                gpp = fields_v[1];
            }
            exit(0);
        }

        k = k + ssize;
    }

    // get the speed of sound by taking derivatives and dP/drho
    eos[0][3] = (eos[1][2] - eos[0][2])/(eos[1][1] - eos[0][1]);
    eos[1][3] = (eos[2][2] - eos[0][2])/(eos[2][1] - eos[0][1]);
    eos[0][5] = eos[0][3]*eos[0][5];
    eos[1][5] = eos[1][3]*eos[1][5];
    
    for (int i=2; i<(npoints-2); ++i) {
        dydx1 = (eos[i+1][2] - eos[i-1][2])/(eos[i+1][1] - eos[i-1][1]);
        dydx2 = (eos[i+2][2] - eos[i-2][2])/(eos[i+2][1] - eos[i-2][1]);
        eos[i][3] = 0.5*(dydx1 + dydx2);
        eos[i][5] = eos[i][3]*eos[i][5];
        if (eos[i][4] < 0) {
            unstable = true;
        }
    }
    
    eos[npoints-2][3] = (eos[npoints-1][2] - eos[npoints-3][2])/(eos[npoints-1][1] - eos[npoints-3][1]);
    eos[npoints-1][3] = (eos[npoints-1][2] - eos[npoints-2][2])/(eos[npoints-1][1] - eos[npoints-2][1]);
    eos[npoints-2][5] = eos[npoints-2][3]*eos[npoints-2][5];
    eos[npoints-1][5] = eos[npoints-1][3]*eos[npoints-1][5];

    // print in case needed
    if (print == true) {
        dmi.print(eos,npoints,9,true,"FSUEOS.txt");
    }
    return 0;
}

//###########################################################
// Bulk Properties class
// ##########################################################

// Calculate the compressibility for FSU+delta model
double bulks :: get_K(double kf, double gss, double gww, double couplings[10]) {
    double res, mstar,en,integral;
    double gsoms2,gwomw2,kappa,lambda,zeta;
    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];

    mstar = mNuc-gss;   // effective masses
    integral = tool.common_integral(kf,mstar);
    en = sqrt(pow(kf,2.0) + pow(mstar,2.0));
    
    res = 3.0*pow(kf,2.0)/en - 6.0/pow(pi,2.0)*pow(kf,3.0)*pow(mstar/en,2.0)*gsoms2*pow(1.0 + gsoms2*(2.0/pow(pi,2.0)*integral + kappa*gss + 0.5*lambda*pow(gss,2.0)),-1.0) 
            + 6.0/pow(pi,2.0)*pow(kf,3.0)*gwomw2*pow(1.0 + 0.5*zeta*gwomw2*pow(gww,2.0),-1.0);
    return res;
}

// Calculate the symmetry energy at saturation (J) for the FSU+delta model
double bulks :: get_J(double kf, double gss, double gww, double couplings[10]) {
    double res,mstar,en,integral;
    double gpomp2,gdomd2,lambda_v,lambda_s;
    
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    mstar = mNuc-gss;
    en = sqrt(pow(kf,2) + pow(mstar,2));   // convenient quanitity to define since used a lot
    integral = tool.common_integral(kf,mstar);
    
    res = pow(kf,2.0)/(6.0*en) - 1.0/(12.0*pow(pi,2.0))*pow(kf,3.0)*pow(mstar/en,2.0)*gdomd2*pow(1.0 + 0.5/pow(pi,2.0)*gdomd2*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0),-1.0)
            + 1.0/(12.0*pow(pi,2.0))*pow(kf,3.0)*gpomp2/(1.0+gpomp2*2.0*lambda_v*pow(gww,2.0));
    
    return res;
}

// Calculate the drivative of the symmetry energy at saturation (L) for the FSU+delta model
double bulks :: get_L(double kf, double gss, double gww, double couplings[10]) {
    double res_n, mstar,en,integral,dSdp,dmstardp,dgdddt,dgssdp,dendp;
    double dIdp, dgdddpdt,dgwwdp, alpha_d, alpha_p, alpha_s, dgppdpdt;
    double gsoms2,gwomw2,gpomp2,gdomd2,kappa,lambda,zeta,lambda_v,lambda_s;

    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    mstar = mNuc - gss;
    en = sqrt(pow(kf,2) + pow(mstar,2));
    integral = tool.common_integral(kf,mstar);
    alpha_s = 1.0 + gsoms2*(2.0/pow(pi,2.0)*integral + kappa*gss + 0.5*lambda*pow(gss,2.0));
    dgssdp = gsoms2*(mstar/en)*pow(alpha_s,-1.0);  // good
    dmstardp = -dgssdp;   // good
    dgdddt = -pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(mstar/en)*pow(1.0 + gdomd2*0.5/pow(pi,2.0)*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0),-1.0);   // good
    dendp = 1.0/en*(0.5*pow(pi,2.0)/kf + mstar*dmstardp); // good
    dgwwdp = gwomw2*pow(1.0 + 0.5*zeta*gwomw2*pow(gww,2.0),-1.0);  // good
    dIdp = 0.5*pow(kf*pi,2.0)/pow(en,3.0) - 3.0*dmstardp*mstar*(-1.0/3.0*pow(kf/en,3.0) - kf/en + log(abs((kf+en)/mstar)));

    alpha_d = 1.0 + gdomd2*1.0/(2.0*pow(pi,2.0))*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0);
    alpha_p = 1.0 + 2.0*lambda_v*gpomp2*pow(gww,2.0);
    
    // looks good
    dgdddpdt =  - 0.5*gdomd2*(mstar/en)*pow(alpha_d,-1) - pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(dmstardp/en - mstar/pow(en,2.0)*dendp)*pow(alpha_d,-1.0)
                + pow(kf,3.0)/(3.0*pow(pi,2.0))*pow(gdomd2,2.0)*(mstar/en)*(1.0/(2.0*pow(pi,2.0))*dIdp + 4.0*lambda_s*gss*dgssdp)*pow(alpha_d,-2.0);

    dgppdpdt = -0.5*gpomp2*pow(alpha_p,-1.0) + pow(kf,3.0)/(3.0*pow(pi,2.0))*pow(gpomp2,2.0)*(4.0*lambda_v*gww*dgwwdp)*pow(alpha_p,-2.0);

    dSdp = pow(pi,2.0)/(6.0*kf)*1.0/en - pow(kf,2.0)/6.0*1.0/pow(en,2.0)*dendp + 0.25*(dmstardp/en - mstar/pow(en,2.0)*dendp)*dgdddt + 0.25*mstar/en*dgdddpdt - 0.25*dgppdpdt;

    res_n = 2.0/pow(pi,2.0)*pow(kf,3.0)*dSdp;
    return res_n;
}

double bulks :: get_Ksym(double kf, double gss, double gww, double couplings[10]) {
    double mstar,en,integral,dmstardp,dgdddt,dgssdp,dendp,d2gwwdp2,var,varprime,d2Idp2;
    double dIdp, dgdddpdt,dgwwdp,d2mstardp2,d2gssdp2,d2endp2,d3gdddp2dt;
    double gsoms2,gwomw2,gpomp2,gdomd2,kappa,lambda,zeta,lambda_v,lambda_s;
    double alpha_s, alpha_d, dalpha_ddp, d2alpha_ddp2, alpha_p, dalpha_pdp, d2alpha_pdp2, d3gppdp2dt;

    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    mstar = mP - gss;   //good
    en = sqrt(pow(kf,2) + pow(mstar,2));    //good
    integral = tool.common_integral(kf,mstar);  //good
    alpha_s = 1.0 + gsoms2*(2.0/pow(pi,2.0)*integral + kappa*gss + 0.5*lambda*pow(gss,2.0)); //good
    dgssdp = gsoms2*(mstar/en)*pow(alpha_s,-1.0);  // good
    dmstardp = -dgssdp;   // good
    dgwwdp = gwomw2*pow(1.0 + 0.5*zeta*gwomw2*pow(gww,2.0),-1.0);  // good
    dgdddt = -pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(mstar/en)*pow(1.0 + gdomd2*0.5/pow(pi,2.0)*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0),-1.0); //good
    dendp = 1.0/en*(0.5*pow(pi,2.0)/kf + mstar*dmstardp); //good
    dIdp = 0.5*pow(kf*pi,2.0)/pow(en,3.0) - 3.0*dmstardp*mstar*(-1.0/3.0*pow(kf/en,3.0) - kf/en + log(abs((kf+en)/mstar))); //good

    d2gssdp2 = gsoms2*1.0/en*(dmstardp - mstar/en*dendp)*pow(alpha_s,-1.0) - pow(gsoms2,2.0)*mstar/en*(2.0/pow(pi,2.0)*dIdp + kappa*dgssdp + lambda*gss*dgssdp)*pow(alpha_s,-2.0); //good
    d2mstardp2 = -d2gssdp2; //good
    d2endp2 = -1.0/pow(en,2.0)*dendp*(pow(pi,2.0)*0.5/kf + mstar*dmstardp) + 1.0/en*(-0.25*pow(pi/kf,4.0) + pow(dmstardp,2.0) + mstar*d2mstardp2); //good
    d2gwwdp2 = -pow(gwomw2,2.0)*zeta*gww*dgwwdp*pow(1.0 + 0.5*zeta*gwomw2*pow(gww,2.0),-2.0); //good

    var = -1.0/3.0*pow(kf/en,3.0) - kf/en + log(abs((kf+en)/mstar)); 
    varprime = -0.5*pow(pi,2.0)/pow(en,3.0) - 0.5*pow(pi/kf,2.0)/en + pow(kf,3.0)*dendp/pow(en,4.0) + kf*dendp/pow(en,2.0)
                -(kf+en)*dmstardp/(mstar*kf+mstar*en) + mstar*(0.5*pow(pi/kf,2.0) + dendp)/(mstar*kf+mstar*en);
    d2Idp2 = pow(pi,4.0)/(2.0*kf*pow(en,3.0)) - 3.0*pow(kf*pi,2.0)*0.5/pow(en,4.0)*dendp - 3.0*d2mstardp2*mstar*var - 3.0*pow(dmstardp,2.0)*var
                -3.0*dmstardp*mstar*varprime;

    alpha_d = 1.0 + gdomd2*1.0/(2.0*pow(pi,2.0))*integral + 2.0*lambda_s*gdomd2*pow(gss,2.0);
    dalpha_ddp = gdomd2*(1.0/(2.0*pow(pi,2.0))*dIdp + 4.0*lambda_s*gss*dgssdp);
    d2alpha_ddp2 = gdomd2*(1.0/(2.0*pow(pi,2.0))*d2Idp2 + 4.0*lambda_s*pow(dgssdp,2.0) + 4.0*lambda_s*gss*d2gssdp2);

    d3gdddp2dt = -pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*1.0/en*(d2mstardp2 - 2.0*dmstardp/en*dendp + 2.0*mstar/pow(en,2.0)*pow(dendp,2.0) - mstar/en*d2endp2)*pow(alpha_d,-1.0)
                -pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*mstar/en*(2.0*pow(alpha_d,-1.0)*pow(dalpha_ddp,2.0) - d2alpha_ddp2)*pow(alpha_d,-2.0)
                + 2.0*pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2/en*(dmstardp - mstar/en*dendp)*dalpha_ddp*pow(alpha_d,-2.0) + gdomd2*mstar/en*dalpha_ddp*pow(alpha_d,-2.0)
                - gdomd2/en*(dmstardp - mstar/en*dendp)*pow(alpha_d,-1.0);

    alpha_p = 1.0 + 2.0*lambda_v*gpomp2*pow(gww,2.0);
    dalpha_pdp = 4.0*lambda_v*gpomp2*gww*dgwwdp;
    d2alpha_pdp2 = 4.0*lambda_v*gpomp2*(pow(dgwwdp,2.0) + gww*d2gwwdp2);
    d3gppdp2dt = -pow(kf,3.0)/(3.0*pow(pi,2.0))*gpomp2*(2.0*pow(alpha_p,-1.0)*pow(dalpha_pdp,2.0) - d2alpha_pdp2)*pow(alpha_p,-2.0) + gpomp2*dalpha_pdp*pow(alpha_p,-2.0);

    dgdddpdt =  - 0.5*gdomd2*(mstar/en)*pow(alpha_d,-1) - pow(kf,3.0)/(3.0*pow(pi,2.0))*gdomd2*(dmstardp/en - mstar/pow(en,2.0)*dendp)*pow(alpha_d,-1.0)
                + pow(kf,3.0)/(3.0*pow(pi,2.0))*pow(gdomd2,2.0)*(mstar/en)*(1.0/(2.0*pow(pi,2.0))*dIdp + 4.0*lambda_s*gss*dgssdp)*pow(alpha_d,-2.0);

    double d2Sdp2 = 1.0/(6.0*en)*(-0.5*pow(pi/kf,4.0) - 2.0*pow(pi,2.0)/(kf*en)*dendp + 2.0*pow(kf/en,2.0)*pow(dendp,2.0) - pow(kf,2.0)/en*d2endp2) 
                    + 0.25/en*( d2mstardp2*dgdddt + 2.0/pow(en,2.0)*pow(dendp,2.0)*mstar*dgdddt -mstar/en*d2endp2*dgdddt + mstar*d3gdddp2dt 
                    - 2.0*mstar/en*dendp*dgdddpdt + 2.0*dmstardp*dgdddpdt - 2.0/en*dendp*dmstardp*dgdddt  ) - 0.25*d3gppdp2dt;

    double dens = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    double Ksym = 9.0*pow(dens,2.0)*d2Sdp2;

    return Ksym;
}

// get bulk properties from the couplings (needs the EOS for SNM)
void bulks :: get_bulkproperties(double couplings[10]) {
    double **arrSNM;
    int ncols = 5; int npoints = 5000;
    double BA, p0, mstar, K, gss, J, Jtilde, kf, L, gww, Ksym;

    // import SNM EOS
    eos.get_SNMEOS(couplings,arrSNM,npoints);
    //dmi.print(arrSNM,npoints,ncols,true,"snm.txt");

    // get saturation properties from array
    BA = dmi.findmin(arrSNM,3,npoints,ncols);
    p0 = arrSNM[dmi.findvalue(arrSNM,npoints,ncols,BA,3,0.01)][0];
    kf = pow(3.0/2.0*pow(pi,2)*p0,1.0/3.0);
    mstar = arrSNM[dmi.findvalue(arrSNM,npoints,ncols,BA,3,0.01)][2]*mNuc;
    cout << "BA: " << BA << endl;
    cout << "p0: " << p0/pow(197.32698,3) << endl;
    cout << "kf: " << kf/197.32698 << endl;
    cout << "mstr/m: " << mstar/mNuc << endl;

    // get properties from analytic expressions
    gss = mNuc - mstar;
    //cout << setprecision(10) << gdd << endl;
    gww = mNuc + BA - sqrt(pow(kf,2.0) + pow(mstar,2.0));
    K = get_K(kf,gss,gww,couplings);
    J = get_J(kf,gss,gww,couplings);
    L = get_L(kf,gss,gww,couplings);
    Ksym = get_Ksym(kf,gss,gww,couplings);
    
    double subp0 = 2.0/(3.0*pow(pi,2.0))*pow(1.15*197.32698,3.0);
    Jtilde = arrSNM[dmi.findvalue(arrSNM,npoints,ncols,subp0,0,0.01)][4];

    cout << "K: " << K << endl;
    cout << "J: " << J << endl;
    cout << "Jtilde: " << Jtilde << endl;
    cout << "L: " << L << endl;
    cout << "Ksym: " << Ksym << endl;
  
    dmi.cleanup(arrSNM,npoints);  // clear the SNM array
}

// returns finite nuclei parameters (15 parameters)
// (gs2,gw2,gp2,gd2,kappa,lambda,zeta,lambda_v,lambda_s,fw,fp,ms,mw,mp,md)
extern "C" {
int get_parameters(double BA, double p0, double Jtilde, double mstar, double K, double L, double Ksym, double zeta, double xi, double lambda_s, double fw, double fp, double masses[4], double fin_couplings[16], bool flag, int gd_sol_type, bool delta_coupling) {
    double a1,a2,a3,b1,c1,c2,c3,g1,integral,tau;
    double kf,gss,gww,gwomw2,en,sdensn,sdensp,sdens,gsoms2,kappa,lambda,gpomp2,lambda_v,gdomd2;

    p0 = p0*pow(197.32698,3); // convert density from 1/fm3 to MeV^3;
    kf = pow(3.0/2.0*pow(pi,2)*p0,1.0/3.0); // get fermi momentum

    gss = mNuc - mstar; // given mstar get gss
    gww = mNuc + BA - sqrt(pow(kf,2) + pow(mstar,2)); // get gww at saturation
    en = sqrt(pow(kf,2) + pow(mstar,2));
    
    // check to see if the coupling constants are realistic 
    if (p0/pow(gww,3.0) < zeta/6.0) {
        cout << "unrealistic: limit is: " << p0/pow(gww,3.0) << " input is: " << zeta << "  " << gww << endl;
        if (flag == true) {
            return -1;
        }
    }

    gwomw2 = gww/(p0-1.0/6.0*zeta*pow(gww,3.0));   // get (gw/mw)^2

    // alphas, betas, gammas
    a1 = gss;
    a2 = 0.5*pow(gss,2.0);
    a3 = 1.0/6.0*pow(gss,3.0);
    b1 = 1.0/24.0*pow(gss,4.0);
    g1 = 1.0;

    // scalar densities
    sdensp = tool.scalardens(kf,mstar);
    sdensn = tool.scalardens(kf,mstar);
    sdens = sdensp + sdensn;

    // algebraically convenient expressions continued
    c1 = sdens;
    c2 = p0*(mNuc+BA) - tool.endensity_integral(kf,mstar) - tool.endensity_integral(kf,mstar) - 0.5*pow(gww,2.0)/gwomw2 - 1.0/8.0*zeta*pow(gww,4.0);
    integral = tool.common_integral(kf,mstar);
    tau = pow(pi,2.0)/(2.0*kf)*en/pow(mstar,2.0) + gwomw2*pow(en/mstar,2.0)*pow(1.0+0.5*zeta*gwomw2*pow(gww,2.0),-1.0) - pow(pi,2.0)/(6.0*pow(kf,3.0))*pow(en/mstar,2.0)*K;
    c3 = pow(tau,-1.0) - 2.0/pow(pi,2.0)*integral;

    // algebraic expressions for the sigma coupling constants
    gsoms2 = -(a2*a2*a2 - 2.0*a1*a2*a3 + a1*a1*b1 + a3*a3*g1 - a2*b1*g1)/(a2*a3*c1 - a1*b1*c1 - a2*a2*c2 + a1*a3*c2 - a3*a3*c3 + a2*b1*c3);
    kappa = -(-a2*a2*c1 + a1*a2*c2 + a2*a3*c3 - a1*b1*c3 + b1*c1*g1 - a3*c2*g1)/(a2*a2*a2 - 2.0*a1*a2*a3 + a1*a1*b1 + a3*a3*g1 - a2*b1*g1);
    lambda = -(a1*a2*c1 - a1*a1*c2 - a2*a2*c3 + a1*a3*c3 - a3*c1*g1 + a2*c2*g1)/(a2*a2*a2 - 2.0*a1*a2*a3 + a1*a1*b1 + a3*a3*g1 - a2*b1*g1);
    
    //double J = get_J_given_Jtilde(Jtilde,kf,BA,mstar,gsoms2,gwomw2,kappa,lambda,zeta,xi,lambda_s,L,Ksym,gd_sol_type,delta_coupling);
    //double J = 50.0;
    double J = Jtilde;
    //cout << "J: " << J << endl;
    
    if (delta_coupling == true) {
        gdomd2 = tool.get_gdomd2(kf,J,L,Ksym,gss,gww,gsoms2,gwomw2,kappa,lambda,zeta,lambda_s,gd_sol_type);
    } else {
        gdomd2 = 0.0;
    }
    
    gpomp2 = tool.get_gpomp2(kf,J,L,gss,gww,gsoms2,gwomw2,gdomd2,kappa,lambda,zeta,lambda_s);
    lambda_v = tool.get_lambda_v(kf,J,gss,gww,gdomd2,gpomp2,lambda_s);
    
    fin_couplings[0] = gsoms2*pow(masses[0],2.0);
    fin_couplings[1] = gwomw2*pow(masses[1],2.0);
    fin_couplings[2] = gpomp2*pow(masses[2],2.0);
    fin_couplings[3] = gdomd2*pow(masses[3],2.0);
    fin_couplings[4] = kappa;
    fin_couplings[5] = lambda;
    fin_couplings[6] = zeta;
    fin_couplings[7] = xi;
    fin_couplings[8] = lambda_v;
    fin_couplings[9] = lambda_s;
    fin_couplings[10] = fw;
    fin_couplings[11] = fp;
    fin_couplings[12] = masses[0];
    fin_couplings[13] = masses[1];
    fin_couplings[14] = masses[2];
    fin_couplings[15] = masses[3];

    if (gpomp2 < 0 || gsoms2 < 0 || gwomw2 < 0 || zeta < 0 || lambda_v < 0) {
        return -1;
    }

    return 0;
}
}

// get the effective compressibility for determining the crust core transition
double tools :: effK2(double kf, double couplings[10], double gss, double gww, double gpp, double gdd, double t) {
    double res, mstarn, mstarp, kp,kn,rp,rn,dens,Yp,vec_termNn,vec_termPp,vec_termPn,dgssdpn,dgssdpp,dgdddpn,dgdddpp,dmunpn,dmuppp,dmuppn;
    double gsoms2, gwomw2, gpomp2, gdomd2, kappa, lambda, zeta, xi, lambda_v,lambda_s,drndpn,drpdpn,drpdpp;
    mstarn = mNuc - gss + 0.5*gdd; // get the effective mass
    mstarp = mNuc - gss - 0.5*gdd;
    kp = kf*pow(1.0-t,1.0/3.0); kn = kf*pow(1.0+t,1.0/3.0); // get the fermi momenta
    rp = sqrt(pow(kp,2) + pow(mstarp,2)); rn = sqrt(pow(kn,2) + pow(mstarn,2));  // convenient quantities to define
    
    gsoms2 = couplings[0];
    gwomw2 = couplings[1];
    gpomp2 = couplings[2];
    gdomd2 = couplings[3];
    kappa = couplings[4];
    lambda = couplings[5];
    zeta = couplings[6];
    xi = couplings[7];
    lambda_v = couplings[8];
    lambda_s = couplings[9];

    // stops errors for 1/0
    if (kp == 0) {
        kp = 1e-15;
    }

    dens = 1.0/(3.0*pow(pi,2.0))*pow(kn,3.0) + 1.0/(3.0*pow(pi,2.0))*pow(kp,3.0);   // total baryon density
    Yp = 0.5*(1.0-t);   // proton fraction

    // algebraic terms

    double vec_denominator = gpomp2*gwomw2*lambda_v*xi*pow(gpp,2.0) + (1.0 + 0.5*zeta*gwomw2*pow(gww,2.0))*(1.0+2.0*gpomp2*lambda_v*pow(gww,2.0)) + pow(gpp,2.0)*(2.0*gwomw2*lambda_v + 0.5*gpomp2*xi + gpomp2*gwomw2*(-12.0*pow(lambda_v,2.0) + 0.25*zeta*xi)*pow(gww,2.0)) ;
    
    vec_termNn = 0.25*(4.0*gwomw2 + gpomp2 + 2.0*(lambda_v+xi)*gwomw2*gpomp2*pow(gpp,2.0) + 16.0*gwomw2*gpomp2*lambda_v*gpp*gww + (0.5*zeta+8.0*lambda_v)*gpomp2*gwomw2*pow(gww,2.0))/vec_denominator;
    
    vec_termPp = 0.25*(4.0*gwomw2 + gpomp2 + 2.0*(lambda_v+xi)*gwomw2*gpomp2*pow(gpp,2.0) - 16.0*gwomw2*gpomp2*lambda_v*gpp*gww + (0.5*zeta+8.0*lambda_v)*gpomp2*gwomw2*pow(gww,2.0))/vec_denominator;
    
    vec_termPn = 0.25*(4.0*gwomw2 - gpomp2 + 2.0*(-lambda_v+xi)*gwomw2*gpomp2*pow(gpp,2.0) + (-0.5*zeta+8.0*lambda_v)*gpomp2*gwomw2*pow(gww,2.0))/vec_denominator;
    
    double Ip = common_integral(kp,mstarp);   // good
    double In = common_integral(kn,mstarn);   // good
    double alpha = Ip/pow(pi,2.0) + In/pow(pi,2.0) + kappa*gss + 0.5*lambda*pow(gss,2.0) + 2.0*lambda_s*pow(gdd,2.0);
    double beta = Ip/(2.0*pow(pi,2.0)) - In/(2.0*pow(pi,2.0)) + 4.0*lambda_s*gss*gdd;
    double gamma = Ip/(4.0*pow(pi,2.0)) + In/(4.0*pow(pi,2.0)) + 2.0*lambda_s*pow(gss,2.0);

    dgssdpn = mstarn/rn*(2.0*gsoms2 + gsoms2*gdomd2*(beta+2.0*gamma))/( 2.0*(1.0 + gsoms2*(alpha-pow(beta,2.0)*gdomd2) + gamma*gdomd2*(1.0+alpha*gsoms2)) );
    dgdddpn = -mstarn/rn*(gdomd2 + gdomd2*gsoms2*(alpha+2.0*beta))/( 2.0*(1.0 + gsoms2*(alpha-pow(beta,2.0)*gdomd2) + gamma*gdomd2*(1.0+alpha*gsoms2)) );
    dgssdpp = mstarp/rp*(2.0*gsoms2 + gsoms2*gdomd2*(2.0*gamma-beta))/( 2.0*(1.0 + gsoms2*(alpha-pow(beta,2.0)*gdomd2) + gamma*gdomd2*(1.0+alpha*gsoms2)) );
    dgdddpp = mstarp/rp*(gdomd2 + gdomd2*gsoms2*(alpha-2.0*beta))/( 2.0*(1.0 + gsoms2*(alpha-pow(beta,2.0)*gdomd2) + gamma*gdomd2*(1.0+alpha*gsoms2)) );
    
    // derivatives for energy
    drndpn = 1.0/rn*(pow(pi,2.0)/kn - mstarn*dgssdpn + 0.5*mstarn*dgdddpn);
    drpdpn = -mstarp/rp*(dgssdpn + 0.5*dgdddpn);
    drpdpp = 1.0/rp*(pow(pi,2.0)/kp - mstarp*dgssdpp - 0.5*mstarp*dgdddpp);

    // derivatives of chemical potentials
    dmunpn = drndpn + vec_termNn; // dmu_n/drho_n
    dmuppp = drpdpp + vec_termPp; // dmu_p/drho_p
    dmuppn = drpdpn + vec_termPn;   // dmu_p/drho_n


    res = dens/4.0*( (dmunpn + 2.0*dmuppn + dmuppp) + 2.0*(1.0-2.0*Yp)*(dmunpn - dmuppp) + pow(1.0-2.0*Yp,2.0)*(dmunpn - 2.0*dmuppn + dmuppp)
        - pow(dmunpn - dmuppp + (1.0-2.0*Yp)*(dmunpn - 2.0*dmuppn + dmuppp),2.0)/(dmunpn - 2.0*dmuppn + dmuppp) );
    return res;
}


// Use the thermodynamic stability method to find the crust core transition point and glue the two together
// has optional ability to print the new EOS
// crust EOS is of the form (nb, en, p, dpde, mub)
int tools :: ThermalCrust(double** crust, double** core, double** &neweos, int nrowscore, int nrowscrust, bool print, int nbcol, int prcol, int Keffcol) {
    int nrows;
    int ncols_core = 9;
    double trdens = 0;
    double x = core[nrowscore-1][Keffcol];
    double xp;

    // find where the compressibility becomes negative
    for (int i=(nrowscore-1);i>0;i=i-1) {
        xp = core[i][Keffcol];  // store next effK
        if (x>0) {          // check if effK starts stable
            if (x*xp<0) {       // check if effK becomes unstable
                trdens = core[i][0];    // mark tr dens
                break;
            }
        }
        x = xp;     // store old effK
    }

    cout << "Transition at: " << trdens << " 1/fm3" << endl; 
    // Get the row number for the crust and core transitions
    int crustrow = dmi.findvalue(crust, nrowscrust,5,trdens,nbcol,1.0);
    int corerow = dmi.findvalue(core,nrowscore,ncols_core,trdens,nbcol,1.0);
    while (core[corerow][nbcol] < crust[crustrow][nbcol]) {     // make sure the density is monotonic
        corerow = corerow + 1;
        if (corerow > (nrowscore-1)) {
            cout << "Thermal crust invalid nb for core" << endl;
            exit(0);
        }
    }

    while (core[corerow][prcol] < crust[crustrow][prcol]) {     // make sure the pressure is monotonic
        corerow = corerow + 1;
        if (corerow > (nrowscore-1)) {
            cout << "Thermal crust invalid nb for core" << endl;
            exit(0);
        }
    }

    // stitch the EOS
    nrows = crustrow + 1 + (nrowscore-corerow-1);
    neweos = new double*[nrows];
    for (int i=0; i<nrows; i++) {
        neweos[i] = new double[5];
    }

    for (int i=0; i<(crustrow+1); ++i) {
        neweos[i][0] = crust[i][0];
        neweos[i][1] = crust[i][1];
        neweos[i][2] = crust[i][2];
        neweos[i][3] = crust[i][3];
        neweos[i][4] = crust[i][4];
    }
        
    for (int i=(corerow+1); i<nrowscore; ++i) {
        neweos[crustrow+1+i-corerow-1][0] = core[i][0];
        neweos[crustrow+1+i-corerow-1][1] = core[i][1];
        neweos[crustrow+1+i-corerow-1][2] = core[i][2];
        neweos[crustrow+1+i-corerow-1][3] = core[i][3];
        neweos[crustrow+1+i-corerow-1][4] = core[i][4];
    }
    
    if (print == true) {
        dmi.print(neweos,nrows,5,true,"FSUEOSC.txt");
    }
    return nrows;
}