#include "finitenuclei.hpp"
#include "NumMethods.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>

using namespace std;
data2 dm3;

int main() {

    auto start = chrono :: high_resolution_clock::now();     
    int A = 48;
    int Z = 20;

    double gs2, gw2, gp2, b, c, h, lambda, mSigma_mev, mOmega_mev, mRho_mev;
    double Observables[5];
    /*
    //FSUGarnet+R
    gs2 = 109.130;
	gw2 = 186.481;
	gp2 = 142.966;
	b = 3.25933/(2.0*939.0);
	c = -0.003285/6.0;
	h = 0.023812/6.0;
	lambda = 0.038274*2.0;
    mSigma_mev = 495.633;
    mOmega_mev = 782.5;
    mRho_mev = 763.0;
    */
    //FSUGarnet+R
    gs2 = 110.349;
	gw2 = 187.695;
	gp2 = 192.927;
	b = 3.26/(2.0*939.0);
	c = -0.003551/6.0;
	h = 0.0235/6.0;
	lambda = 0.043377*2.0;
    mSigma_mev = 496.939;
    mOmega_mev = 782.5;
    mRho_mev = 763.0;
    
    /*
    // FSUGold2+R
    gs2 = 103.760;
	gw2 = 169.410;
	gp2 = 128.301;
	b = 3.79239/(2.0*939.0);
	c = -0.010635/6.0;
	h = 0.011660/6.0;
	lambda = 0.0316212*2.0;
    mSigma_mev = 501.611;
    mOmega_mev = 782.5;
    mRho_mev = 763.0;
    */

    int field_grid = 4000;
    int meson_grid = 2000; // field_grid/meson_grid must be an integer

    hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma_mev,mOmega_mev,mRho_mev,A,Z,25,field_grid,meson_grid,2,Observables);
    /*
    // MCMC sample for charge radii
    int npoints = 1000;
    int k;
    double** couplings;
    ofstream out("MCMC_FSUGARNET_XEFT_Finite.txt");
    dm3.importdata("FSUGARNET_XEFT_params.txt",couplings);
    for (int i=1; i<=npoints; ++ i) {
        k = i*10-1;
        mSigma_mev = couplings[k][7];
        gs2 = couplings[k][0]*pow(mSigma_mev,2.0);
        gw2 = couplings[k][1]*pow(mOmega_mev,2.0);
        gp2 = couplings[k][2]*pow(mRho_mev,2.0);
        b = couplings[k][3];
        c = couplings[k][4];
        h = couplings[k][5];
        lambda = couplings[k][6];
        hartee_method(gs2,gw2,gp2,b,c,h,lambda,mSigma_mev,mOmega_mev,mRho_mev,A,Z,25,field_grid,meson_grid,2,Observables);

        // output (BA, Rn, Rp, Rch, Rwk, Rn-Rp, Rch-Rwk)
        cout << Observables[0] << "  " << Observables[1] << "  " << Observables[2] << "  " << Observables[3] << "  " << Observables[4] << endl;
        out << Observables[0] << "  " << Observables[1] << "  " << Observables[2] << "  " << Observables[3] << "  " << Observables[4]
            << "  " << Observables[1]-Observables[2] << "  " << Observables[3]-Observables[4] << endl;
    }
    dm3.cleanup(couplings,10000);
    */
    auto stop = chrono :: high_resolution_clock::now();
    auto duration = chrono :: duration_cast<chrono :: milliseconds>(stop - start);
    cout << setprecision(5) << duration.count()/1000.0 << "s" << endl;
    return 0;
    
}