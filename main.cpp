#include "finitenuclei.hpp"
#include "NumMethods.hpp"
#include "MCMC.hpp"
#include "infinitematter.hpp"
#include "Conversions.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <assert.h>
#include <omp.h>
#include "minpack.hpp"

using namespace std;
data2 dm3;
nummeth nmm;
bulks bulk1;
tools tool1;
equationofstate eosm;
Convert convm;

const double pi = 4.0*atan(1.0);

int main() {
    auto start = chrono :: high_resolution_clock::now();     
    
    double Observables[7]; 
    double inf_couplings[10]; 
    double fin_couplings[16];
    int gridsize = 401;
    srand(time(0));
    
    /*
    double** param_sets;
    dm3.importdata("param_sets_FSU.txt",param_sets);
    for (int i=0; i<16; ++i) {
        fin_couplings[i] = param_sets[49][i];
    }
    dm3.cleanup(param_sets,50);
    */

    //FSU230529
    //double params[8] = {502.2303545, 100.25574844165376, 159.94894689710125, 83.08017667927881, 4.475103985, -0.01870805824, 0.000399945274706, 0.004248329329141};
    
    //FSU BigApple
    //double params[8] = {492.7300000, 93.507400000000000, 151.68390000000000, 200.5562000000000, 5.203260000, -0.02173900000, 0.000700000000000, 0.047471000000000};
    
    //FSUGarnet
    double params[8] = {496.939, 110.349, 187.695, 192.927, 3.26, -0.003551, 0.0235, 0.043377};
    
    // FSUGold2
    //double params[8] = {497.479, 108.0943, 183.7893, 80.4656, 3.0029, -0.000533, 0.0256, 0.000823};

    // FSUGold2+R
    //double params[8] = {501.611, 103.760, 169.410, 128.301, 3.79239, -0.010635, 0.01166, 0.0316212};

    // FSUGarnet+R
    //double params[8] = {495.633, 109.130, 186.481, 142.966, 3.25933, -0.003285, 0.02381, 0.038274};

    /*
    // Experimental
    fin_couplings[0] = params[1]; // gs2
    fin_couplings[1] = params[2]; // gw2
    fin_couplings[2] = params[3]; // gp2
    fin_couplings[3] = 0.0; // gd2
    fin_couplings[4] = params[4]; // kappa
    fin_couplings[5] = params[5]; // lambda
    fin_couplings[6] = params[6]; // zeta
    fin_couplings[7] = 0.0; // xi
    fin_couplings[8] = params[7]; // Lambda_v
    fin_couplings[9] = 0.0; // Lambda_s
    fin_couplings[10] = 0.0; // fw
    fin_couplings[11] = 0.0;    // fp
    fin_couplings[12] = params[0];
    fin_couplings[13] = 782.5;
    fin_couplings[14] = 763.0;
    fin_couplings[15] = 980.0;
    */
    // Experimental
    fin_couplings[0] = params[1]; // gs2
    fin_couplings[1] = params[2]; // gw2
    fin_couplings[2] = params[3]*pow(2000.0/763.0,2.0); // gp2
    fin_couplings[3] = 0.0; // gd2
    fin_couplings[4] = params[4]; // kappa
    fin_couplings[5] = params[5]; // lambda
    fin_couplings[6] = params[6]; // zeta
    fin_couplings[7] = 0.0; // xi
    fin_couplings[8] = params[7]; // Lambda_v
    fin_couplings[9] = 0.0; // Lambda_s
    fin_couplings[10] = 0.0; // fw
    fin_couplings[11] = -0.5;    // fp
    fin_couplings[12] = params[0];
    fin_couplings[13] = 782.5;
    fin_couplings[14] = 2000.0;
    fin_couplings[15] = 980.0;
    
    for (int i=0; i<16; ++i) {
        cout << fin_couplings[i] << "  ";
    }
    cout << endl;

    tool1.convert_to_inf_couplings(fin_couplings, inf_couplings);
    bulk1.get_bulkproperties(inf_couplings);

    for (int i=0; i<10; ++i) {
        cout << inf_couplings[i] << "  ";
    }
    cout << endl;

    /*
    double BA, kf, p0, mstar, K, J_tilde, L, Ksym, zeta, xi, lambda_s, fw, fp;
    double masses[4];
    //           
    // Experimental Model
    masses[0] = 502.2303545;          
    masses[1] = 782.5; 
    masses[2] = 763.0; 
    masses[3] = 980.0;
    BA = -16.3253261;   // Binding energy (MeV) decreases binding proportionally (barely effects charge and skin)
    kf = 1.31138;                                    // larger values lower binding and (lower skin and charge radii equally)
    p0 = 0.1518791179; //2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);        
    J_tilde = 36.63400779;        // Symmetry energy at sat (MeV)
    mstar = 0.5981984388;      // Effective Mass (MeV)
    K = 246.6814529;             // Compressibility (MeV) 
    L = 105.67548212;               // Derivative of Symmetry Energy at sat (MeV)
    Ksym = -58.91130412652875;
    zeta = 0.000399945274706;
    xi = 0.0;                                            // Self interaction strength for w meson
    fw = 0.0;
    fp = 0.0;
    lambda_s = 0.0;
    bulk1.get_parameters(BA,p0,J_tilde,mstar*939.0,K,L,Ksym,zeta,xi,lambda_s,fw,fp,masses,fin_couplings,true,1,false);     // calculate coupling constants FSU Model  
    
    for (int i=0; i<16; ++i) {
        cout << fin_couplings[i] << "  ";
    }
    cout << endl;
    //cout << BA/-16.3 << "  " << kf/1.30 << "  " << mstar/(939*0.61) << "  " << K/230.0 << "  " << J/32.59 << "  " << L/60.50 << endl;
    tool1.convert_to_inf_couplings(fin_couplings, inf_couplings);
    //chisq(fin_couplings);
    */
    /*
    //double** Symm_EOS;
    double** PNM_EOS;
    int npoints = 200;
    eosm.get_PNMEOS(inf_couplings,PNM_EOS,npoints);
    //eosm.get_SymmetryEnergy(inf_couplings,Symm_EOS,npoints);
    //dm3.print(Symm_EOS,npoints,3,true,"DINOa_Symm.txt");
    dm3.print(PNM_EOS,npoints,2,true,"FSU-Garnet.txt");
    //dm3.cleanup(Symm_EOS,npoints);
    dm3.cleanup(PNM_EOS,npoints);
    */
    
    //hartree_method(fin_couplings,16,8,20,gridsize,3,Observables,1.3,false,false,0.0);
    //hartree_method(fin_couplings,40,20,20,gridsize,3,Observables,1.2,false,false,0.0);
    hartree_method(fin_couplings,48,20,20,gridsize,3,Observables,1.2,true,false,0.0);
    //hartree_method(fin_couplings,68,28,20,gridsize,3,Observables,1.2,false,false,0.0);
    //hartree_method(fin_couplings,90,40,20,gridsize,3,Observables,1.2,false,false,0.0);
    //hartree_method(fin_couplings,100,50,20,gridsize,3,Observables,1.4,false,false,0.0);
    //hartree_method(fin_couplings,116,50,20,gridsize,3,Observables,1.2,false,false,0.0);
    //hartree_method(fin_couplings,132,50,20,gridsize,3,Observables,1.2,false,false,0.0);
    //hartree_method(fin_couplings,144,62,20,gridsize,3,Observables,1.3,true,true,0.0);
    //hartree_method(fin_couplings,208,82,20,gridsize,3,Observables,1.2,true,false,0.0);
    
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
    /*
    double** init;
    dm3.importdata("startfile_FSUMid.txt",init);
    Optimize(200,"invcovmatrix_FSUGARNET.txt",init,"Optimal_Ksym50,fw10,fp.txt","exp_data_small.txt",1);
    dm3.cleanup(init,15);
    */
    /*
    double** init;
    dm3.importdata("excel_calib_cont.txt",init);
    excel_calibrate(init,"excel_analysis.txt",1,true);
    dm3.cleanup(init,15);
    */
    /*
    double ** temp_array; double p0, BA, K, J, L, mstar, Ksym, h, fw, fp, fc;
    double masses[5];
    string file = "Optimal_Ksymn100,fw10.txt";
    dm3.importdata(file,temp_array);
    int nrows = dm3.rowcount(file);
    int ncols = dm3.colcount(file);
    ofstream out("Compare_Ksymn100,fw10,gd2.txt");
    for (int i=0; i<nrows; ++i) {
        p0 = 2/(3.0*pow(pi,2.0))*pow(temp_array[i][1],3.0);
        BA = temp_array[i][0];
        J = temp_array[i][4];
        mstar = temp_array[i][2]*939.0;
        K = temp_array[i][3];
        L = temp_array[i][5];
        Ksym = temp_array[i][6];
        h = temp_array[i][7];
        fw = temp_array[i][8];
        fp = temp_array[i][9];
        fc = temp_array[i][10];
        masses[0] = temp_array[i][11];
        masses[1] = temp_array[i][12];
        masses[2] = temp_array[i][13];
        masses[3] = temp_array[i][14];
        bulk1.get_parameters(BA, p0,J,mstar,K,L,Ksym,h,fw,fp,fc,masses,fin_couplings,false,-1);
        for (int j=0; j<ncols; ++j) {
            out << temp_array[i][j] << "  ";
        }
        out << fin_couplings[3] << endl;
    }   
    */

    // ######################### Neutron Star Calcultions ##########################################
    
    // GET EOS FOR INFINITE MATTER
    /*
    double** CORE_EOS;
    int npoints = 150;
    eosm.get_EOS_NSM(inf_couplings,CORE_EOS,npoints,true,false);
    double Urca_dens = nmm.Urca_threshold(CORE_EOS,9,npoints,7,8,0);
    cout << "Urca_onset_dens = " << Urca_dens << endl;
    //dm3.cleanup(CORE_EOS,npoints);
    
    
    // ADD CRUST
    double** crust; double** EOS;
    dm3.importdata("dat_files/CRUSTEOS.txt",crust);
    int nrows = dm3.rowcount("dat_files/CRUSTEOS.txt");
    int n = tool1.ThermalCrust(crust,CORE_EOS,EOS,npoints,nrows,true,0,2,6);
    dm3.cleanup(crust,nrows);
    dm3.cleanup(CORE_EOS,npoints);
    dm3.cleanup(EOS,n);          // comment out if using the EOS to calc MR
    */
    
    /*
    // convert for RNS/Lorene
    double** array;
    dm3.importdata("BigAppleNStarEOS-167lines.txt",array);
    ofstream out("RNS_FSU_EOS.d");
    int N = dm3.rowcount("BigAppleNStarEOS-167lines.txt");
    for (int i=0; i<N; ++i) {
        out << convm.energyCONV(0,2)*array[i][1] << " " << convm.energyCONV(0,5)*array[i][2] << endl;
    }
    dm3.cleanup(array,N);
    */
    /*
    // CALCULATE MASS/RADIUS PROFILE
    // ---------------- Calculate Mass Radius Profile
    string mrfilename = "FSU";               // output MR filename title (adds _RM)
    int encol = 1;                          // input en col
    int prcol = 2;                          // input pr col
    int dpdecol = 3;
    double cv = convm.energyCONV(0,1);   // MeV/fm3 to unitless
    // --------------------------------------------------
    
    //int n = N;
    nmm.pretovconv(EOS, encol, prcol, cv, cv, n);                                 //convert to unitless en, pr
    double pr0 = 2.0*cv; // start at 2 MeV/fm3
    cout << pr0 << endl;

    double MR[2];
    //nmm.tovsolve(50.0*cv, 1e-4, EOS, n, 4, encol, prcol, dpdecol, MR, true);
    nmm.multitov(1e-4, EOS, n, 4, encol, prcol, dpdecol, 1000, mrfilename,pr0);    // calculate mass radius profile, output is (n,en,pr,r,m) (mev/fm3)
    dm3.cleanup(EOS,n);        // comment out if calc ILQ
    */
    /*
    // To specify what mass of star for the Tidal Deformability
    double** MRarr;
    string RMfile = mrfilename + "_RM.txt";
    double starmass = 1.40;
    dm3.importdata(RMfile,MRarr);
    int nrowsRM = dm3.rowcount(RMfile);
    double cpr = dm3.interpolate(nrowsRM,5,MRarr,starmass,4,2,true);
    double Urca_mass = dm3.interpolate(nrowsRM,5,MRarr,Urca_dens,0,4,true);
    cout << "Urca mass = " << Urca_mass << endl;
    dm3.cleanup(MRarr,nrowsRM);
    cout << cpr << " MeV/fm3" << endl;
    /*
    // ---------------------- Calculate Relativistic and Newtonian ILQ
    encol = 1;                          // input en col
    prcol = 2;                          // input pr col
    dpdecol = 3;
    int np = 100;                      // npoints for ILQ
    cv = convm.energyCONV(0,1);   // MeV/fm3 to unitless
    double ILQR[4];                   // (I,L,Q, R)
    //nm.pretovconv(EOS, encol, prcol, cv, cv, n);                     // convert to unitless en, pr if didnt convert already in MR
    // ----------------------------------------------------
    
    // Single ILQ
    double icp = cpr*cv;                                    // central density to calculate 
    cout << "Relativistic TD: " << endl;
    nmm.RloveMCMC(EOS,dpdecol,encol,prcol,1e-5,ILQR,icp,n,5);                      // Calculate I, L, Q
    dm3.cleanup(EOS,n);
    */

    // RBM Methods
    //RBM_generate_fields(208,82,"training_DINO.txt");
    //sample_param_space(200,"dat_files/startfile_dino.txt");
    int A[10] = {16,40,48,68,90,100,116,132,144,208};
    int Z[10] = {8 ,20,20,28,40,50 ,50 ,50 ,62 ,82 };
    //get_Observables("validation_DINO.txt",208,82);
    
    for (int i=0; i<6; ++i) {
        //RBM_generate_fields(A[i],Z[i],"training_DINO.txt");
        //get_Observables("validation_DINO.txt",A[i],Z[i]);
    }
    
    //RBM_error_check("RBM_samples.txt",8);


    // MCMC methods
    //MCMC_NS(1000,10000,"dat_files/invcovmatrix_RBM.txt","dat_files/CRUSTEOS.txt");
    //MCMC_FN(0,0,"dat_files/exp_data.txt");
    //MCMC_Observables("MCMC.txt","dat_files/CRUSTEOS.txt");
    
    // Get parameters from set of bulk properties
    /*
    double** MCMC_bulks;
    dm3.importdata("MCMC_burnin1.txt",MCMC_bulks);
    double masses[4] = {500.0,782.5,763.0,980.0};
    double kf, p0;
    ofstream out("DINO_calibration_set.txt");

    for (int i=0; i<36; ++i) {
        masses[0] = MCMC_bulks[i][8];          
        kf = MCMC_bulks[i][1];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);  
        bulk1.get_parameters(MCMC_bulks[i][0],p0,MCMC_bulks[i][4],MCMC_bulks[i][2]*939.0,MCMC_bulks[i][3],MCMC_bulks[i][5],MCMC_bulks[i][6],MCMC_bulks[i][7],0.0,0.0,0.0,0.0,masses,fin_couplings,true,1,true);
        for (int j=0; j<16; ++j) {
            out << fin_couplings[j] << "  ";
        }
        out << endl;
    }
    */

    // remove leftover files
    remove("Ap.txt"); remove("Bp.txt"); remove("Fn.txt"); remove("Gn.txt");
    remove("meson_fields.txt"); remove("neutron_spectrum.txt"); remove("proton_spectrum.txt");

    auto stop = chrono :: high_resolution_clock::now();
    auto duration = chrono :: duration_cast<chrono :: milliseconds>(stop - start);
    cout << setprecision(5) << duration.count()/1000.0 << "s" << endl;
    return 0;
}