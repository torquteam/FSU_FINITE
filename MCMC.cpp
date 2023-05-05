#include "NumMethods.hpp"
#include "finitenuclei.hpp"
#include "MCMC.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;
data2 dm1;
nummeth nm1;

const double pi = 4.0*atan(1.0);
const double mP = 939; //938.27231;    // Mass of proton (MeV)
const double mN = 939; //939.56542052;
const double mNuc = (mP+mN)/2.0;
const double mOmega = 782.5;
const double mRho = 763.0;

const int field_grid = 2000;
const int meson_grid = 1000; // field_grid/meson_grid must be an integer

double rand_uniform(double cntr, double w) {
    double x,y;
    x = rand()*1.0/RAND_MAX;
    y = cntr + 2.0*w*x - w;
    return y;
}

// Calculate (g_rho/m_rho)^2 Analytically
double get_gpomp2(double kf, double asym, double L, double gss, double gsoms2, double gww, double gwomw2, double h, double b, double c) {
    double gpomp2,z1,z2, mstarp,mstarn,rn,rp,intp,intn,integral,gsdsdk,gwdwdk,alpha,beta,ap,an,gamma;
    mstarp = mP - gss; mstarn = mN - gss;   // effective masses

    // convenient quantities to define (makes expressions easier to read)
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));
    z2 = pow(kf,3.0)/(12.0*pow(pi,2.0)); z1 = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn);

    // integrals for proton and neutron
    intp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));
    intn = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));
    integral = intp + intn;

    // g_sigma*dsigma/dk and g_omega*d_omega/dk
    gsdsdk = 2.0/pow(pi,2.0)*pow(kf,2.0)* gsoms2/2.0*(mstarp/rp + mstarn/rn)*pow(1.0+ gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    gwdwdk = gwomw2*2.0*pow(kf/pi,2.0)/(1.0 + 3.0*gwomw2*h*pow(gww,2.0));

    // more convenient quantities to define (no physical meaning)
    alpha = 1.0/z2*pow(rn*rp,3.0)*gwdwdk*pow(kf,4.0)*pow(asym-z1,2.0);
    beta = asym*pow(rp*rn*kf,3.0)*(-1.5*gww+gwdwdk*kf);
    ap = pow(rp,3.0)*gww*z2*pow(pi*kf,2.0)*(-pow(rn,2.0) + 0.5*pow(kf,2.0) - 0.5*kf*mstarn*gsdsdk);
    an = pow(rn,3.0)*gww*z2*pow(pi*kf,2.0)*(-pow(rp,2.0) + 0.5*pow(kf,2.0) - 0.5*kf*mstarp*gsdsdk);
    gamma = pow(rp*rn,3.0)*(1.5*gww*pow(kf,3.0)*z1 - gwdwdk*pow(kf,4.0)*z1 + 6.0*pow(pi,2.0)*gww*L*z2);
    
    gpomp2 = alpha/(beta+ap+an+gamma);
    return gpomp2;
}

// Get the energy density for the FSU model
double get_en(double kf, double t, double gss, double gsoms2, double gww, double gwomw2, double gpp, double gpomp2, double b, double c, double h, double lambda) {
    double integralp = 0; double integraln = 0; double res,mss2,mww2,mpp2,mstrp,mstrn,kfn,kfp,rn,rp;

    mss2 = pow(gss,2.0)/gsoms2;  // (m_sigma*sigma)^2
    mww2 = pow(gww,2.0)/gwomw2;  // (m_omgea*omega)^2
    mpp2 = pow(gpp,2.0)/gpomp2;  // (m_rho*rho)^2
    mstrp = mP-gss; mstrn = mN-gss;   // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);     // define the fermi momentum for the proton and neutron given t
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));       // convenient quanitity to define since used a lot

    // integrals for the energy density
    integralp = kfp*pow(rp,3)/4.0 - pow(mstrp,2)*kfp*rp/8.0 - pow(mstrp,4)/8.0*log(kfp+rp) + pow(mstrp,4)/8.0*log(abs(mstrp));
    integraln = kfn*pow(rn,3)/4.0 - pow(mstrn,2)*kfn*rn/8.0 - pow(mstrn,4)/8.0*log(kfn+rn) + pow(mstrn,4)/8.0*log(abs(mstrn));
    
    // add free contributions
    res = 0.5*mss2 + 0.5*mww2 + +0.5*mpp2 + 1.0/3.0*b*mNuc*pow(gss,3.0) + 0.25*c*pow(gss,4.0) + 0.75*h*pow(gww,4.0) + 1.5*lambda*pow(gpp,2.0)*pow(gww,2.0)
            + 1.0/pow(pi,2)*integralp + 1.0/pow(pi,2)*integraln;
    return res;
}

// Calculate the compressibility 
double get_K(double kf, double gss, double gwomw2, double gsoms2, double b, double c, double h, double gww) {
    double integralp = 0; double integraln = 0; double integral = 0; double res, mstrp,mstrn,rn,rp;

    mstrp = mP-gss; mstrn = mN-gss;   // effective masses
    rn = sqrt(pow(kf,2) + pow(mstrn,2)); rp = sqrt(pow(kf,2) + pow(mstrp,2));     // convenient quanitity to define since used a lot

    integralp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrp,2.0))/rp + 3.0*pow(mstrp,2.0)/2.0*log(abs(mstrp/(kf+rp)));
    integraln = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrn,2.0))/rn + 3.0*pow(mstrn,2.0)/2.0*log(abs(mstrn/(kf+rn)));
    integral = integralp + integraln;
    res = gwomw2/(1.0 + 3.0*gwomw2*h*pow(gww,2.0))*6.0*pow(kf,3)/pow(pi,2) + 3.0/2.0*pow(kf,2.0)*(1.0/rp + 1.0/rn) 
            - 3.0/2.0*gsoms2*pow(kf,3.0)/pow(pi,2.0)*pow(mstrp/rp + mstrn/rn,2.0)*pow(1.0+gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    return res;
}

// Calculate the symmetry energy at saturation (J)
double get_J(double kf, double gss, double gpomp2, double gww, double lambda) {
    double res,mstarp,mstarn,rn,rp;
    
    mstarp = mP - gss; mstarn = mN - gss; // effective masses
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));   // convenient quanitity to define since used a lot
    res = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn) + 1.0/(12.0*pow(pi,2.0))*gpomp2/(1.0+gpomp2*pow(gww,2.0)*lambda)*pow(kf,3.0);
    return res;
}

// Calculate the drivative of the symmetry energy at saturatiob (L)
double get_L(double kf, double gss, double gsoms2, double gww, double gwomw2, double gpomp2, double h, double lambda, double b, double c) {
    double res, mstarp,mstarn,rn,rp,intp,intn,integral,gsdsdk,gwdwdk;

    mstarp = mP - gss; mstarn = mN - gss; // effective masses
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));   // convenient quanitity to define since used a lot
    intp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));    // proton integral
    intn = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));    // neutron integral
    integral = intp + intn;

    // g_sigma*dsigma/dk
    gsdsdk = 2.0/pow(pi,2.0)*pow(kf,2.0)* gsoms2/2.0*(mstarp/rp + mstarn/rn)*pow(1.0+ gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    
    // g_omega*domega/dk
    gwdwdk = gwomw2*2.0*pow(kf/pi,2.0)/(1.0 + 3.0*gwomw2*h*pow(gww,2.0));
    
    res = pow(kf,2.0)/6.0*(1.0/rp+1.0/rn) + pow(kf,3.0)/12.0*( (mstarp*gsdsdk - kf)/pow(rp,3.0) + (mstarn*gsdsdk - kf)/pow(rn,3.0) ) 
        + pow(kf,3.0)/(4.0*pow(pi,2.0))*gpomp2/(1.0+gpomp2*pow(gww,2.0)*lambda) - pow(kf,4.0)/(6.0*pow(pi,2.0))*pow(gpomp2,2.0)*gww*lambda*gwdwdk/pow(1.0 + gpomp2*lambda*pow(gww,2.0),2.0);
    return res;
}

int get_parameters(double BA, double p0, double J, double mstar, double K, double L, double h, double mSigma, double fp, double params[9], bool flag) {
    double a1,a2,a3,b1,b2,b3,c1,c2,c3,g1,g2,g3,z1,z2,m1,m2,m3,m4,n1,n2,n3,n4;
    double kf,gss,mstarp,mstarn,gww,gwomw2,rn,rp,sdensn,sdensp,sdens,gsoms2,b,c,gpomp2,gpp,lambda,pintegralK,nintegralK,integralK;

    p0 = p0*pow(197.32698,3); // convert density from 1/fm3 to MeV^3;
    kf = pow(3.0/2.0*pow(pi,2)*p0,1.0/3.0); // get fermi momentum

    gss = mNuc - mstar; // given mstar get gss
    mstarp = mP - gss;  mstarn = mN - gss;  // get effective masses (technically redundant since mN = mP)
    gww = mNuc + BA - 0.5*sqrt(pow(kf,2) + pow(mstarp,2)) - 0.5*sqrt(pow(kf,2) + pow(mstarn,2));    // get gww at saturation
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2)); // convenient quantities to define
    
    // check to see if the coupling constants are realistic 
    if (p0/pow(gww,3.0) < h) {
        cout << "unrealistic: limit is: " << p0/pow(gww,3.0) << " input is: " << h << "  " << gww << endl;
        if (flag == true) {
            exit(0);
        }
    }

    gwomw2 = gww/(p0-h*pow(gww,3.0));   // get (gw/mw)^2

    // proton and neutron integrals
    pintegralK = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));
    nintegralK = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));
    integralK = pintegralK + nintegralK;

    // algebraically convenient expressions
    n1 = 6.0*pow(kf,3.0)/pow(pi,2.0)*gwomw2/(1.0+ 3.0*gwomw2*h*pow(gww,2.0));
    n2 = 3.0/2.0*pow(kf,2.0)*(1.0/rp + 1.0/rn);
    n3 = 3.0/2.0*pow(kf,3.0)/pow(pi,2.0)*pow(mstarp/rp + mstarn/rn,2.0);    
    n4 = 1.0/pow(pi,2.0)*integralK;
    a1 = gss;
    a2 = mNuc*pow(gss,2.0);
    a3 = pow(gss,3.0);
    b1 = K - n1 - n2;
    b2 = 2.0*mNuc*gss*(K-n1-n2);
    b3 = 3.0*pow(gss,2.0)*(K-n1-n2);
    g1 = 0.5*pow(gss,2.0);
    g2 = 1.0/3.0*mNuc*pow(gss,3.0);
    g3 = 0.25*pow(gss,4.0);
    m1 = 1.0/pow(pi,2.0)* (kf*pow(rp,3)/4.0 - pow(mstarp,2)*kf*rp/8.0 - pow(mstarp,4)/8.0*log(kf+rp) + pow(mstarp,4)/8.0*log(abs(mstarp)));
    m2 = 1.0/pow(pi,2.0)* (kf*pow(rn,3)/4.0 - pow(mstarn,2)*kf*rn/8.0 - pow(mstarn,4)/8.0*log(kf+rn) + pow(mstarn,4)/8.0*log(abs(mstarn)));
    m3 = 0.5*pow(gwomw2,-1.0)*pow(gww,2.0);
    m4 = 0.75*h*pow(gww,4.0);

    // scalar densities
    sdensp = 1.0/pow(pi,2)*mstarp*(kf*rp/2.0 - pow(mstarp,2)/2.0*log(abs(kf+rp)) + pow(mstarp,2)/2.0*log(abs(mstarp)) );
    sdensn = 1.0/pow(pi,2)*mstarn*(kf*rn/2.0 - pow(mstarn,2)/2.0*log(abs(kf+rn)) + pow(mstarn,2)/2.0*log(abs(mstarn)) );
    sdens = sdensp + sdensn;

    // algebraically convenient expressions continued
    c1 = sdens;
    c2 = -n3 - n4*(K-n1-n2);
    c3 = p0*(mNuc+BA)-m1-m2-m3-m4;

    // algebraic expressions for the sigma coupling constants
    gsoms2 = (a3*b2*g1-a2*b3*g1-a3*b1*g2+a1*b3*g2+a2*b1*g3-a1*b2*g3)/(c3*a3*b2-c3*a2*b3-c2*a3*g2+c1*b3*g2+c2*a2*g3-c1*b2*g3);
    b = (c3*a3*b1-c3*a1*b3-c2*a3*g1+c1*b3*g1+c2*a1*g3-c1*b1*g3)/(-a3*b2*g1+a2*b3*g1+a3*b1*g2-a1*b3*g2-a2*b1*g3+a1*b2*g3);
    c = (c3*a2*b1-c3*a1*b2-c2*a2*g1+c1*b2*g1+c2*a1*g2-c1*b1*g2)/(a3*b2*g1-a2*b3*g1-a3*b1*g2+a1*b3*g2+a2*b1*g3-a1*b2*g3);
    
    // get(gp/mp)^2
    gpomp2 = get_gpomp2(kf,J,L,gss,gsoms2,gww,gwomw2,h,b,c);
    gpp = 0.0;   // valid at saturation

    z2 = pow(kf,3.0)/(12.0*pow(pi,2.0));
    z1 = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn);
    lambda = 1.0/(gpomp2*pow(gww,2.0))*(z2*gpomp2/(J-z1) - 1.0);  // calculate last remaining coupling

    params[0] = gsoms2*pow(mSigma,2.0);
    params[1] = gwomw2*pow(mOmega,2.0);
    params[2] = gpomp2*pow(mRho,2.0);
    params[3] = b;
    params[4] = c;
    params[5] = h;
    params[6] = lambda;
    params[7] = mSigma;
    params[8] = fp;

    return 0;
}

// target liklihood function for the calibration of models to finite nuclei observables
// takes in the proposed bulk properties and the old bulk properties to compute a chi square value
double targetNuclei(double* props, double* olds, double** cov, double** DATA, bool flag) {
    double chisqp = 0; double chisq0 = 0;   // chisq of prior
    double vec1[7]; double vec2[7];         // temporary arrays for matrix mulitplication
    double BA, p0, kf, J, mstar, K, L, h, mSigma, fp, res;   // bulk parameters
    double gs2, gw2, gp2, b, c, lambda;
    double oldparams[9]; double newparams[9];       // arrays for the calculated parameters
    double lkl0, lklp, max0, maxp;
    int A, Z;
    double Observables[5]; double RnRp48o, RnRp208o, RnRp48p, RnRp208p;

    //Make sure proposed changes are physical
    BA = props[0]*DATA[0][3]; kf = props[1]*DATA[1][3]; mstar = props[2]*DATA[2][3]*mNuc; 
    K = props[3]*DATA[3][3]; J = props[4]*DATA[4][3]; L = props[5]*DATA[5][3]; h = props[6]*DATA[6][3];
    mSigma = props[7]*DATA[7][3]; fp = props[8]*DATA[8][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    get_parameters(BA,p0,J,mstar,K,L,h,mSigma,fp,newparams,flag);     // solve for parameters given bulk properties
    // Check if lambda is negative (eff rho mass is imaginary)
    if (newparams[6] < 0 || newparams[5] < 0) {
        return -1.0;
    }

    // initialize the temp arrays and get the prior distribution X^2 = (x-mu)^T Cov (x-mu)
    vec1[0] = 0; vec1[1] = 0; vec1[2] = 0; vec1[3] = 0; vec1[4] = 0; vec1[5] = 0; vec1[6] = 0;
    vec2[0] = 0; vec2[1] = 0; vec2[2] = 0; vec2[3] = 0; vec2[4] = 0; vec2[5] = 0; vec2[6] = 0;
    for (int i=0; i<7; ++i) {
        for (int j=0; j<7; ++j) {
            vec1[i] = vec1[i] + cov[i][j]*(olds[j]-DATA[j][2]);
            vec2[i] = vec2[i] + cov[i][j]*(props[j]-DATA[j][2]);
        }
        chisq0 = chisq0 + vec1[i]*(olds[i]-DATA[i][2]);
        chisqp = chisqp + vec2[i]*(props[i]-DATA[i][2]);
    }
    double prior = exp(chisq0/2.0 - chisqp/2.0);
    
    // BEGIN PARALLEL
    //#pragma omp parallel sections private(BA,kf,J,mstar,K,L,h,mSigma,fp,p0,gs2,gw2,b,c,lambda,A,Z) 
    //{
    //#pragma omp section
    //{
    BA = olds[0]*DATA[0][3]; kf = olds[1]*DATA[1][3]; mstar = olds[2]*DATA[2][3]*mNuc;      // get bulk properties by rescaling 
    K = olds[3]*DATA[3][3]; J = olds[4]*DATA[4][3]; L = olds[5]*DATA[5][3]; h = olds[6]*DATA[6][3];
    mSigma = olds[7]*DATA[7][3]; fp = olds[8]*DATA[8][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0); // get density
    //cout << BA << "  " << kf << "  " << mstar << "  " << K << "  " << J << "  " << L << "  " << h << "  " << mSigma << "  " << fp << endl;
    get_parameters(BA,p0,J,mstar,K,L,h,mSigma,fp,oldparams,flag);     // solve for parameters given bulk properties
    
    gs2 = oldparams[0];
    gw2 = oldparams[1];
    gp2 = oldparams[2];
    b = oldparams[3];
    c = oldparams[4];
    h = oldparams[5];
    lambda = oldparams[6];
    mSigma = oldparams[7];
    fp = oldparams[8];
    
    A = 48; Z = 20; // Ca48
    hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma,mOmega,mRho,fp,A,Z,25,field_grid,meson_grid,2,Observables);
    RnRp48o = sqrt(Observables[1]) - sqrt(Observables[2]);

    A = 208; Z = 82; // Pb208
    hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma,mOmega,mRho,fp,A,Z,25,field_grid,meson_grid,2,Observables);
    RnRp208o = sqrt(Observables[1]) - sqrt(Observables[2]);

    //}
    //#pragma omp section
    //{
    // Proposed change
    BA = props[0]*DATA[0][3]; kf = props[1]*DATA[1][3]; mstar = props[2]*DATA[2][3]*mNuc; 
    K = props[3]*DATA[3][3]; J = props[4]*DATA[4][3]; L = props[5]*DATA[5][3]; h = props[6]*DATA[6][3];
    mSigma = props[7]*DATA[7][3]; fp = props[8]*DATA[8][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    get_parameters(BA,p0,J,mstar,K,L,h,mSigma,fp,newparams,flag);     // solve for parameters given bulk properties

    gs2 = newparams[0];
    gw2 = newparams[1];
    gp2 = newparams[2];
    b = newparams[3];
    c = newparams[4];
    h = newparams[5];
    lambda = newparams[6];
    mSigma = newparams[7];
    fp = newparams[8];
    cout << gs2 << "  " << gw2 << "  " << gp2 << "  " << b << "  " << c << "  " << h << "  " << lambda << "  " << mSigma << "  " << fp << endl;
    
    A = 48; Z = 20; // Ca48
    hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma,mOmega,mRho,fp,A,Z,25,field_grid,meson_grid,2,Observables);
    RnRp48p = sqrt(Observables[1]) - sqrt(Observables[2]);

    A = 208; Z = 82; // Pb208
    hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma,mOmega,mRho,fp,A,Z,25,field_grid,meson_grid,2,Observables);
    RnRp208p = sqrt(Observables[1]) - sqrt(Observables[2]);
    //}
    //}   // END PARALLEL

    lkl0 = lkl0*exp(-0.5*pow((0.121 - RnRp48o)/0.026,2.0))*exp(-0.5*pow((0.283 - RnRp208o)/0.071,2.0));
    lklp = lklp*exp(-0.5*pow((0.121 - RnRp48p)/0.026,2.0))*exp(-0.5*pow((0.283 - RnRp208p)/0.071,2.0));

    res = prior*lklp/lkl0;
    return res;
}

// DATA matrix specifies the current means of the model and appropriate scaling factors for the bulk properties along with good starting guesses and widths for MCMC
// output is of the form (BA, kf, mstar, K, J, L, h, mSigma, fp)
int Optimize(int runs, string covdata, double** DATA, int burnin, string outname) {
    srand(time(0));
    int counts[9]; double stds[9]; double arate[9]; double agoal[9]; double paramtest[9];
    double r, a, BA, kf, J, mstar, K, L, h, p0, mSigma, fp;
    ofstream out(outname);
    ofstream lout("lastpoint.txt");
    double *olds = new double[9];      // (bulks, radius1, radius2, ... , icp1, icp2, ... , TD1, TD2,  ... , Mmax, icpmax)
    double *cands = new double[9];     // (bulks, radius1, radius2, ... , icp1, icp2, ... , TD1, TD2,  ... , Mmax, icpmax)

    double** cov;
    dm1.importdata(covdata,cov);    // import covariance matrix

    // acceptance rate goals, starting point for bulk properties and their widths, acceptance counter, averages
    for (int i=0; i<9; ++i) {
        agoal[i] = 0.5;
        counts[i] = 0;
        olds[i] = DATA[i][0];
        stds[i] = DATA[i][1];
    }

    // rescale the bulk properties
    BA = olds[0]*DATA[0][3]; kf = olds[1]*DATA[1][3]; mstar = olds[2]*DATA[2][3]*mNuc; 
    K = olds[3]*DATA[3][3]; J = olds[4]*DATA[4][3]; L = olds[5]*DATA[5][3]; h = olds[6]*DATA[6][3];
    mSigma = olds[7]*DATA[7][3]; fp = olds[8]*DATA[8][3];
    
    // make sure theres no issue with initial guess for params
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    get_parameters(BA,p0,J,mstar,K,L,h,mSigma,fp,paramtest,true);
    if (paramtest[6] < 0 || paramtest[5] < 0) {
        cout << "problem" << endl;
        dm1.cleanup(cov,9);
        dm1.cleanup(DATA,9);
        delete olds;
        delete cands;
        return 0;
    } 

    // Burn in phase
    for (int i=0; i<burnin; ++i) {
        cout << i << "  " << endl;
        // go through each bulk parameter and suggest a random change
        for (int j=0; j<9; ++j) {
            // candidate points = old points
            for (int k=0; k<9; ++k) {
                cands[k] = olds[k];
            }

            cands[j] = rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = targetNuclei(cands,olds,cov,DATA,true);   // probability of accepance
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
            }

            // monitor the rate every 50 points
            if ((i+1)%10 == 0) {
                
                arate[j] = 1.0*counts[j]/10.0;
                counts[j] = 0;
                
                if (arate[j] < agoal[j]) {
                    stds[j] = 0.9*stds[j];                 // if acceptance rate is too low then decrease the range
                } else if (arate[j] > agoal[j]) {
                    stds[j] = 1.1*stds[j];                 // if acceptance rate is too high then increase the range
                }
                
                // save arate, stds, olds every 50 iterations (in case of code breaking or computer freezing can be deleted with better computer)
                for (int k=0; k<9; ++k) {
                    lout << "arate: " << arate[k] << "  ";
                }
                lout << endl;
                for (int k=0; k<9; ++k) {
                    lout << "stds: " << stds[k] << "  ";
                }
                lout << endl;
                for (int k=0; k<9; ++k) {
                    lout << "olds: " << olds[k] << "  ";
                }
                lout << endl;
            }
        }
        // print the rate, width, and bulk values every run
        cout << "arate: ";
        for (int k=0; k<9; ++k) {
            cout << arate[k] << "  ";
        }
        cout << endl;
        cout << "stds: ";
        for (int k=0; k<9; ++k) {
            cout << stds[k] << "  ";
        }
        cout << endl;
        cout << "olds: ";
        for (int k=0; k<9; ++k) {
            cout << olds[k] << "  ";
        }
        cout << endl;
        
        // save MCMC runs
        out << setprecision(10);
        for (int k=0; k<9; ++k) {                   // print params
            out << olds[k]*DATA[k][3] << "  ";
        }
        out << endl;
    }

    // print out the bulk properties values and widths at end of burn in
    cout << "--------------------------------------------------------------------" << endl;
    out << "---------------------------------------------------------------------" << endl;
    out << olds[0] << "  " << olds[1] << "  " << olds[2] << "  " << olds[3] << "  " << olds[4] << "  " << olds[5] << "  " << olds[6] << "  " << olds[7] << "  " << olds[8] << endl;
    out << stds[0] << "  " << stds[1] << "  " << stds[2] << "  " << stds[3] << "  " << stds[4] << "  " << stds[5] << "  " << stds[6] << "  " << stds[7] << "  " << stds[8] << endl;
    cout << "--------------------------------------------------------------------" << endl;
    out << "---------------------------------------------------------------------" << endl;
    
    // MCMC RUN
    for (int i=0; i<runs; ++i) {
        // go through each bulk parameter and suggest a random change
        for (int j=0; j<7; ++j) {
            // candidate points = old points
            for (int k=0; k<9; ++k) {
                cands[k] = olds[k];
            }

            cands[j] = rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = targetNuclei(cands,olds,cov,DATA,true);         // probability of acceptance
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
            }
        }
        
        cout << "olds: ";
        for (int k=0; k<9; ++k) {
            cout << olds[k]*DATA[k][3] << "  ";
        }
        cout << endl;

        // save MCMC runs
        out << setprecision(10);
        for (int k=0; k<9; ++k) {
            out << olds[k]*DATA[k][3] << "  ";
        }
    }
    
    // print out final rate
    cout << "arate: ";
    for (int k=0; k<9; ++k) {
        cout << 1.0*counts[k]/runs << "  ";
    }
    cout << endl;
    
    dm1.cleanup(cov,7);
    delete olds;
    delete cands;
    
    return 0;
}