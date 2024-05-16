#include "NumMethods.hpp"
#include "finitenuclei.hpp"
#include "infinitematter.hpp"
#include "MCMC.hpp"
#include "Conversions.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>

using namespace std;
data2 dm1;
nummeth nm1;
bulks bulkmc;
tools toolmc;
equationofstate eosmc;
Convert convmc;

const double pi = 4.0*atan(1.0);
const double mP = 939; //938.27231;    // Mass of proton (MeV)
const double mN = 939; //939.56542052;
const double mNuc_mev = (mP+mN)/2.0;
const double mE = 0.511;
const double mMU = 105.7;

const double hbar_mevfm = 197.32698; // MeV fm
const double r0_fm = 1.25; //fm (arbitrary)
const double enscale_mev = pow(hbar_mevfm,2)/(2.0*mNuc_mev*pow(r0_fm,2)); // arbitrary

const int gridsize = 401;

double rand_uniform(double cntr, double w) {
    double x,y;
    x = rand()*1.0/RAND_MAX;
    y = cntr + 2.0*w*x - w;
    return y;
}

//Box muller method
// returns a value from the normal distribution centered at some mean and stddev
double rand_normal(double mean, double stddev) {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double get_lkl_prior(double* proposed_params, double** prior_covariance, double** prior_DATA) {
    double vec[7];
    double chisq = 0;

    // initialize the temp arrays and get the prior distribution X^2 = (x-mu)^T Cov (x-mu)
    vec[0] = 0; vec[1] = 0; vec[2] = 0; vec[3] = 0; vec[4] = 0; vec[5] = 0; vec[6] = 0;
    for (int i=0; i<7; ++i) {
        for (int j=0; j<7; ++j) {
            vec[i] = vec[i] + prior_covariance[i][j]*(proposed_params[j]-prior_DATA[j][2]);
        }
        chisq = chisq + vec[i]*(proposed_params[i]-prior_DATA[i][2]);
    }
    
    double prior = exp(-chisq/2.0);
    return prior;
}

//Make sure proposed changes are physical (add Ksym)
int rho_mass_check(double* proposed_params, double** prior_DATA, double masses[4], int gd_sol_type,bool delta_coupling) {
    double BA, p0, kf, J, mstar, K, L, Ksym, zeta,xi,lambda_s, fw, fp;   // bulk parameters
    double fin_couplings[15];

    BA = proposed_params[0]*prior_DATA[0][3]; 
    kf = proposed_params[1]*prior_DATA[1][3]; 
    mstar = proposed_params[2]*prior_DATA[2][3]*mNuc_mev; 
    K = proposed_params[3]*prior_DATA[3][3]; 
    J = proposed_params[4]*prior_DATA[4][3]; 
    L = proposed_params[5]*prior_DATA[5][3]; 
    Ksym = proposed_params[6]*prior_DATA[6][3];
    zeta = proposed_params[7]*prior_DATA[7][3];
    xi = proposed_params[8]*prior_DATA[8][3];
    lambda_s = proposed_params[9]*prior_DATA[9][3];
    fw = proposed_params[10]*prior_DATA[10][3];
    fp = proposed_params[11]*prior_DATA[11][3];

    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    bulkmc.get_parameters(BA,p0,J,mstar,K,L,Ksym,zeta,xi,lambda_s,fw,fp,masses,fin_couplings,true,gd_sol_type,delta_coupling);     // solve for parameters given bulk properties
    // Check if lambda is negative (eff rho mass is imaginary)
    if (fin_couplings[7] < 0 || fin_couplings[6] < 0) {
        return -1;
    } else {
        return 0;
    }
}

// count the number of parameters that are going to be varied (just counts 1s that appear)
// return an array with the parameter indicies
int counter(double** priorDATA, int n_params, double* &index_array, int toggle_col) {
    int count = 0;
    for (int i=0; i<n_params; ++i) {
        count = count + priorDATA[i][toggle_col];
    }

    index_array = new double[count];
    int subcounter = 0;
    for (int i=0; i<n_params; ++i) {
        if (priorDATA[i][toggle_col] == 1) {
            index_array[subcounter] = i;
            subcounter = subcounter + 1;
        }
    }

    return count;
}

// #############################################################################################
// ################################################################################################
// RBM Calibration Stuff

// output set of wave functions and fields for a given parameter set and nucleus
int RBM_generate_fields(int A, int Z, string params_file) {

    // import parameters to be used in RBM sample
    double** param_set; string** narray; string** parray;
    dm1.importdata(params_file, param_set);
    int num_param_sets = dm1.rowcount(params_file);

    // variable declarations
    int nstates_n, nstates_p;
    double Observables[7]; 
    double fin_couplings[16];
    int ncols_meson = 8;
    int npoints_meson = gridsize;

    // initialize the arrays for meson fields, energies, and wave functions
    double** Fn_unitless; double** Gn_unitless;
    double** Ap_unitless; double** Bp_unitless;
    double** meson_fields_unitless;

    for (int k=0; k<16; ++k) {
        fin_couplings[k] = param_set[0][k];
    }

    int exit_code = hartree_method(fin_couplings,A,Z,20,gridsize,3,Observables,1.3,false,true);
    if(exit_code != 0) {
        exit(0);
    }
    

    nstates_n = dm1.rowcount("neutron_spectrum.txt");
    nstates_p = dm1.rowcount("proton_spectrum.txt");
    dm1.importdata_string("neutron_spectrum.txt",narray);
    dm1.importdata_string("proton_spectrum.txt",parray);

    for (int k=0; k<nstates_n; ++k) {
        ofstream fout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/f_wave/state" + narray[k][6] + ".txt",ios::out);
        ofstream gout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/g_wave/state" + narray[k][6] + ".txt",ios::out);
        fout.close();
        gout.close();
    }

    for (int k=0; k<nstates_p; ++k) {
        ofstream aout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/c_wave/state" + parray[k][6] + ".txt",ios::out);
        ofstream bout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/d_wave/state" + parray[k][6]+ ".txt",ios::out);
        aout.close();
        bout.close();
    }

    // initialize the output files
    ofstream sout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "sigma.txt");
    ofstream oout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "omega.txt");
    ofstream rout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "rho.txt");
    ofstream dout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "delta.txt");
    ofstream aout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "coulomb.txt");
    
    ofstream enout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/neutron/energies.txt");
    ofstream epout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/proton/energies.txt");
    ofstream rvecout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/rvec.txt");
    
    dm1.importdata("Fn.txt",Fn_unitless);
    rvecout << 0.0 << endl;
    for (int i=1; i<gridsize; ++i) {
        rvecout << scientific << setprecision(5) << Fn_unitless[i][0] << endl;
    }

    dm1.cleanup_string(narray,nstates_n);
    dm1.cleanup_string(parray,nstates_p);
    // loop through each nuclei and solve the high fidelity problem
    for (int i=0; i<num_param_sets; ++i) {
        for (int k=0; k<16; ++k) {
            fin_couplings[k] = param_set[i][k];
        }

        exit_code = hartree_method(fin_couplings,A,Z,20,gridsize,3,Observables,1.3,false,true);
        if(exit_code != 0) {
            exit(0);
        }

        nstates_n = dm1.rowcount("neutron_spectrum.txt");
        nstates_p = dm1.rowcount("proton_spectrum.txt");
        dm1.importdata_string("neutron_spectrum.txt",narray);
        dm1.importdata_string("proton_spectrum.txt",parray);
        dm1.importdata("Fn.txt",Fn_unitless);
        dm1.importdata("Gn.txt",Gn_unitless);
        dm1.importdata("Ap.txt",Ap_unitless);
        dm1.importdata("Bp.txt",Bp_unitless);

        for (int k=0; k<nstates_n; ++k) {
            ofstream fout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/f_wave/state" + narray[k][6] + ".txt", ios::app);
            ofstream gout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/g_wave/state" + narray[k][6] + ".txt", ios::app);
            gout << 0.0 << "  ";
            fout << 0.0 << "  ";
            for (int j=1; j<gridsize; ++j) {
                fout << scientific << setprecision(10) << Fn_unitless[j][k+1] << "  ";
                gout << scientific << setprecision(10) << Gn_unitless[j][k+1] << "  ";
            }
            fout << endl;
            gout << endl;
        }

        for (int k=0; k<nstates_p; ++k) {
            ofstream aout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/c_wave/state" + parray[k][6] + ".txt", ios::app);
            ofstream bout("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/d_wave/state" + parray[k][6] + ".txt", ios::app);
            aout << 0.0 << "  ";
            bout << 0.0 << "  ";
            for (int j=1; j<gridsize; ++j) {
                aout << scientific << setprecision(10) << Ap_unitless[j][k+1] << "  ";
                bout << scientific << setprecision(10) << Bp_unitless[j][k+1] << "  ";
            }
            aout << endl;
            bout << endl;
        }

        dm1.importdata("meson_fields.txt",meson_fields_unitless);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,0,1.0/r0_fm);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,1,1.0/enscale_mev);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,2,1.0/enscale_mev);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,3,1.0/enscale_mev);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,4,1.0/enscale_mev);
        dm1.convert_array(meson_fields_unitless,npoints_meson,ncols_meson,5,1.0/enscale_mev);
        remove("meson_fields.txt");
        cout << "pass: " << i+1 << endl;

        for (int j=0; j<npoints_meson; ++j) {
            sout << scientific << setprecision(10) << meson_fields_unitless[j][1] << "  ";
            oout << scientific << setprecision(10) << meson_fields_unitless[j][2] << "  ";
            rout << scientific << setprecision(10) << meson_fields_unitless[j][3] << "  ";
            dout << scientific << setprecision(10) << meson_fields_unitless[j][4] << "  ";
            aout << scientific << setprecision(10) << meson_fields_unitless[j][5] << "  ";
        }
        sout << endl;
        oout << endl;
        rout << endl;
        dout << endl;
        aout << endl;

        for (int k=0; k<nstates_n; ++k) {
            enout << scientific << setprecision(10) << narray[k][0] << "  ";
        }
        enout << endl;

        for (int k=0; k<nstates_p; ++k) {
            epout << scientific << setprecision(10) << parray[k][0] << "  ";
        }
        epout << endl;

        // cleanup
        dm1.cleanup(meson_fields_unitless,npoints_meson);
    
        // cleanup
        dm1.cleanup(Fn_unitless,gridsize);
        dm1.cleanup(Gn_unitless,gridsize);
        dm1.cleanup(Ap_unitless,gridsize);
        dm1.cleanup(Bp_unitless,gridsize); 
        dm1.cleanup_string(narray,nstates_n);
        dm1.cleanup_string(parray,nstates_p);  
    }

    dm1.importdata_string("neutron_spectrum.txt",narray);
    dm1.importdata_string("proton_spectrum.txt",parray);

    for (int k=0; k<nstates_n; ++k) {
        dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/f_wave/state" + narray[k][6] + ".txt");
        dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/neutron/g_wave/state" + narray[k][6] + ".txt");

    }

    for (int k=0; k<nstates_p; ++k) {
        dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/c_wave/state" + parray[k][6] + ".txt");
        dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data" + "/proton/d_wave/state" + parray[k][6] + ".txt");
    
    }      

    // transpose the meson fields into column vectors
    dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "sigma.txt");
    dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "omega.txt");
    dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "rho.txt");
    dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "delta.txt");
    dm1.transpose_file("/home/msals97/Desktop/FSU_FINITE/ReducedBasisMethods/" + to_string(A) + "," + to_string(Z) + "/" + to_string(A) + "," + to_string(Z) + ",Data/meson_fields/" + "coulomb.txt");

    dm1.cleanup(param_set,num_param_sets);
    dm1.cleanup_string(narray,nstates_n);
    dm1.cleanup_string(parray,nstates_p); 
    return 0;
}

// generate a random sample
int generate_sample(double params[16], double** start_data) {
    double bulks[16];

    // randomnly sample parameter space
    for (int i=0; i<16; ++i) {
        bulks[i] = start_data[i][0];
        if (start_data[i][1] == 1) {
            bulks[i] = rand_normal(start_data[i][0], start_data[i][2]);
        }
    }

    // meson masses
    double masses[4];
    for (int i=0; i<4; ++i) {
        masses[i] = bulks[12+i];
    }
    bulkmc.get_parameters(bulks[0],bulks[1],bulks[4],bulks[2]*mNuc_mev,bulks[3],bulks[5],bulks[6],bulks[7],bulks[8],bulks[9],bulks[10],bulks[11],masses,params,false,1,false);

    for (int i=0; i<16; ++i) {
        cout << params[i] << "  ";
    }
    cout << endl;
    
    for (int i=0; i<16; ++i) {
        cout << bulks[i] << "  ";
    }
    cout << endl;
    return 0;
}

// generate samples of parameters based on some distribution specified in file
double sample_param_space(int num_sets, string startfile) {
    double fin_couplings[16]; double Observables[7];
    int exit_code;
    srand(time(0)); // random seed

    ofstream out("param_sets.txt");

    double** start_data;
    dm1.importdata(startfile,start_data);

    for (int i=0; i<num_sets; ++i) {
        // Generate a random set of couplings and make sure gp2 and lambda_v > 0
        generate_sample(fin_couplings,start_data);
        while (fin_couplings[8]<0 || fin_couplings[2]<0) {
            generate_sample(fin_couplings,start_data);
        }

        exit_code = hartree_method(fin_couplings,16,8,20,gridsize,3,Observables,1.2,false,false);
        if (exit_code != 0) {
            dm1.cleanup(start_data,16);
            exit(0);
        }

        exit_code = hartree_method(fin_couplings,48,20,20,gridsize,3,Observables,1.2,false,false);
        if (exit_code != 0) {
            dm1.cleanup(start_data,16);
            exit(0);
        }

        exit_code = hartree_method(fin_couplings,208,82,20,gridsize,3,Observables,1.2,false,false);
        if (exit_code != 0) {
            dm1.cleanup(start_data,16);
            exit(0);
        }

        for (int j=0; j<16; ++j) {
            out << fin_couplings[j] << "  ";
        }
        out << endl;
    }
    dm1.cleanup(start_data,16);
    return 0;
}

// Used to compute ground state observables to test RBM error
void get_Observables(string param_set, int A, int Z) {
    ofstream observ_out(to_string(A)+ "," + to_string(Z) + "Observables.txt");
    ofstream enout(to_string(A) + "," + to_string(Z) + "energies_n.txt");
    ofstream epout(to_string(A) + "," + to_string(Z) + "energies_p.txt");

    double** param_set_arr;
    string** narray; string** parray; string** nref; string** pref;
    dm1.importdata(param_set,param_set_arr);
    int n_sets = dm1.rowcount(param_set);
    double fin_couplings[16];
    double Observables[7];

    int exit_code;
    int nstates_n, nstates_p;

    for (int i=0; i<n_sets; ++i) {
        for (int j=0; j<16; ++j) {
            fin_couplings[j] = param_set_arr[i][j];
        }
        
        exit_code = hartree_method(fin_couplings,A,Z,20,gridsize,3,Observables,1.2,false,false);
        nstates_n = dm1.rowcount("neutron_spectrum.txt");
        nstates_p = dm1.rowcount("proton_spectrum.txt");
        dm1.importdata_string("neutron_spectrum.txt",narray);
        dm1.importdata_string("proton_spectrum.txt",parray);
        dm1.importdata_string("spectrum_ref_files/" + to_string(A) + "," + to_string(Z) + "," + "neutron_spectrum.txt",nref);
        dm1.importdata_string("spectrum_ref_files/" + to_string(A) + "," + to_string(Z) + "," + "proton_spectrum.txt",pref);

        for (int k=0; k<nstates_n; ++k) {
            for (int l=0; l<nstates_n; ++l) {
                if (narray[l][6] == nref[k][0]) {
                    enout << scientific << setprecision(10) << narray[l][0] << "  ";
                    break;
                }
            }
        }
        enout << endl;

        for (int k=0; k<nstates_p; ++k) {
            for (int l=0; l<nstates_p; ++l) {
                if (parray[l][6] == pref[k][0]) {
                    epout << scientific << setprecision(10) << parray[l][0] << "  ";
                    break;
                }
            }
        }
        epout << endl;
        if(exit_code != 0) {
            dm1.cleanup(param_set_arr,n_sets);
            dm1.cleanup_string(narray,nstates_n);
            dm1.cleanup_string(parray,nstates_p);  
            exit(0);
        }
        observ_out << Observables[0] << "  " << Observables[3] << "  " << Observables[5] << endl;
        dm1.cleanup_string(narray,nstates_n);
        dm1.cleanup_string(parray,nstates_p);  
    }

    dm1.cleanup(param_set_arr,n_sets);
}

void RBM_error_check(string RBM_file, int n_params) {
    double** RBM_data;
    dm1.importdata(RBM_file,RBM_data);
    int nrows = dm1.rowcount(RBM_file);
    int A[10] = {16,40,48,68,90,100,116,132,144,208};
    int Z[10] = {8,20,20,28,40,50,50,50,62,82};
    double mw = 782.5; double mp = 763.0; double md = 980.0;
    double gd2 = 0.0; double xi = 0.0; double lambda_s = 0.0; double fw = 0.0; double fp = 0.0;
    double fin_couplings[16] = {0, 0, 0, gd2, 0, 0, 0, xi, 0, lambda_s, fw, fp, 0.0, mw, mp, md};
    double Observables[7];
    vector<double> hf_results;
    double error[22];

    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<10; ++j) {
            // ms, gs2, gw2, gp2, kappa, lambda, zeta, lambda_v
            fin_couplings[12] = RBM_data[i][0]; // ms
            fin_couplings[0] = RBM_data[i][1];  // gs2
            fin_couplings[1] = RBM_data[i][2];  // gw2
            fin_couplings[2] = RBM_data[i][3];  // gp2
            fin_couplings[4] = RBM_data[i][4];  // kappa
            fin_couplings[5] = RBM_data[i][5];  // lambda
            fin_couplings[6] = RBM_data[i][6];  // zeta
            fin_couplings[8] = RBM_data[i][7];  // lambda_v
            hartree_method(fin_couplings,A[j],Z[j],20,gridsize,3,Observables,1.2,false,false);
            if (A[j] == 48 || A[j] == 208) {
                hf_results.push_back(Observables[0]);
                hf_results.push_back(Observables[3]);
                hf_results.push_back(Observables[5]);
            } else {
                hf_results.push_back(Observables[0]);
                hf_results.push_back(Observables[3]);
            }
        }
        cout << i << "  ";
        for (int k=0; k<22; ++k) {
            error[k] = (i*error[k] + abs(RBM_data[i][n_params+k] - hf_results[k]))/(i+1);
            cout << error[k] << "  ";
        }
        cout << endl;
    }
}

// MCMC for Neutron Stars
// ################################################################################################
// ################################################################################################
// bulks has to be of form [BA, kf, mstar/m, K, J, L, zeta]
int param_change(int n_params, vector<double>& bulks_0, vector<double>& bulks_p, vector<double>& stds, double mw, double mp, int index, double inf_couplings[10]) {
    double ms = 500.0;
    double md = 980.0;
    double masses[4] = {ms,mw,mp,md};
    double fin_couplings[16];

    // copy old values
    for (int i=0; i<n_params; ++i) {
        bulks_p[i] = bulks_0[i];    
    }
    bulks_p[index] = rand_normal(bulks_0[index], stds[index]);

    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(bulks_p[1],3.0);
    int flag = bulkmc.get_parameters(bulks_p[0],p0,bulks_p[4],bulks_p[2]*mNuc_mev,bulks_p[3],bulks_p[5],0,bulks_p[6],0.0,0.0,0.0,0.0,masses,fin_couplings,true,1,false);
    if (flag == -1) {
        toolmc.convert_to_inf_couplings(fin_couplings,inf_couplings);
        return -1;
    } else {
        toolmc.convert_to_inf_couplings(fin_couplings,inf_couplings);
        return 0;
    }
}

double compute_prior(double** invcov, double means[7], vector<double>& bulks) {
    double vec1[7];         // temporary arrays for matrix mulitplication
    double chisq = 0.0;
    double scaling[7] = {-16.3,1.30,0.61,230.0,32.59,60.50,0.06};

    // initialize the temp arrays and get the prior distribution X^2 = (x-mu)^T Cov^-1 (x-mu)
    vec1[0] = 0; vec1[1] = 0; vec1[2] = 0; vec1[3] = 0; vec1[4] = 0; vec1[5] = 0; vec1[6] = 0;
    for (int i=0; i<7; ++i) {
        for (int j=0; j<7; ++j) {
            vec1[i] = vec1[i] + invcov[i][j]*(bulks[j]/scaling[j]-means[j]);
        }
        chisq = chisq + vec1[i]*(bulks[i]/scaling[i]-means[i]);
    }
    return exp(-chisq/2.0);
}

double compute_lkl(double inf_couplings[10], double** CRUST, int nrowscrust, int flag) {
    double lkl = 1.0;
    double chisq = 0.0;
    double Mmax_exp = 3.0;
    double** COREEOS; double** NS_EOS;
    int npoints = 250;
    double cv = convmc.energyCONV(0,1);   // MeV/fm3 to unitless
    vector<vector<double>> TOV_out;

    if (flag == -1) {
        cout << "flagged" << endl;
        return 0.0;
    }

    eosmc.get_EOS_NSM(inf_couplings,COREEOS,npoints,false,false);
    int n = toolmc.ThermalCrust(CRUST,COREEOS,NS_EOS,npoints,nrowscrust,false,0,2,6);
    nm1.pretovconv(NS_EOS,1,2,cv,cv,n);
    TOV_out = nm1.multitov(1e-4,NS_EOS,n,4,1,2,3,200,"NA",2.0*cv);
    double Mmax_th = dm1.findmax_vec(TOV_out,4,TOV_out.size(),5);
    cout << Mmax_th << endl;
    chisq = pow(Mmax_th - Mmax_exp, 2.0)/pow(0.05,2.0);
    lkl = lkl*exp(-chisq/2.0);

    dm1.cleanup(COREEOS,npoints);
    dm1.cleanup(NS_EOS,n);
    return lkl;
}

double metropolis(double post0, double postp, vector<double>& bulks_0, vector<double>& bulks_p, vector<int>& acc_counts, int index, int n_params) {
    double r = 1.0*rand()/RAND_MAX;
    double a = postp/post0;
    if (a>1 || post0==0) {
        a = 1.0;
    }
    if (r <= a) {
        post0 = postp;
        for (int i=0; i<n_params; ++i) {
            bulks_0[i] = bulks_p[i];
        }
        acc_counts[index] = acc_counts[index] + 1;
    }
    return post0;
}

void adaptive_width(int iter, int n_check, vector<double>& arate, vector<int>& acc_counts, vector<double>& stds, double agoal, int index) {
    if ((iter+1)%n_check == 0) {
        arate[index] = acc_counts[index]*1.0/n_check;
        acc_counts[index] = 0;
        if (arate[index] < agoal) {
            stds[index] = 0.9*stds[index];
        } else if (arate[index] > agoal) {
            stds[index] = 1.1*stds[index];
        }
    }
}

void MCMC_NS(int nburnin, int nruns, string covdata, string crust) {
    srand(time(0));
    int n_params = 7;
    double inf_couplings[10]; double fin_couplings[16];
    double ms = 500.0; double mw = 782.5; double mp = 763.0; double md = 980.0;
    double masses[4] = {ms,mw,mp,md};
    double lklp, postp;
    int n_check = 50;
    double agoal = 0.3;
    ofstream out("MCMC_burnin.txt");
    ofstream aout("MCMC.txt");
    int flag;

    // initialization
    vector<double> bulks_0 = {-16.280,1.306,0.593,238.0,37.62,112.80,0.0256};
    vector<double> bulks_p = {-16.280,1.306,0.593,238.0,37.62,112.80,0.0256};
    vector<double> stds = {0.02,0.005,0.004,2.8,1.11,16.1,0.002};
    vector<int> acc_counts = {0,0,0,0,0,0,0};
    vector<double> arate = {0,0,0,0,0,0,0};
    double prior_means[7] = {0.99878396660424529, 1.00491981968270430, 0.97168943949904518, 1.03461649348865610, 1.15443471952649210, 1.86516512236962420, 0.42744836577894457};
    double** invcov; double** CRUST;
    dm1.importdata(covdata,invcov);
    dm1.importdata(crust,CRUST);
    int nrowscrust = dm1.rowcount(crust);
    
    // MCMC start
    double prior = compute_prior(invcov,prior_means,bulks_0);
    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(bulks_0[1],3.0);
    flag = bulkmc.get_parameters(bulks_0[0],p0,bulks_0[4],bulks_0[2]*mNuc_mev,bulks_0[3],bulks_0[5],0,bulks_0[6],0.0,0.0,0.0,0.0,masses,fin_couplings,true,1,false);
    toolmc.convert_to_inf_couplings(fin_couplings,inf_couplings);
    double lkl0 = compute_lkl(inf_couplings,CRUST,nrowscrust,flag);
    double post0 = prior*lkl0;


    for (int i=0; i<nburnin; ++i) {
        for (int j=0; j<n_params; ++j) {
            // get new proposed params
            flag = param_change(n_params,bulks_0,bulks_p,stds,mw,mp,j,inf_couplings);
            cout << bulks_p[0] << "  " << bulks_p[1] << "  " << bulks_p[2] << "  " << bulks_p[3] << "  " << bulks_p[4] << "  " << bulks_p[5] << "  " << bulks_p[6] << endl;

            // get posterior
            prior = compute_prior(invcov,prior_means,bulks_p);
            lklp = compute_lkl(inf_couplings,CRUST,nrowscrust,flag);
            postp = prior*lklp;
            cout << post0 << "  " << postp << endl;
            
            // metropolis hastings step
            post0 = metropolis(post0,postp,bulks_0,bulks_p,acc_counts,j,n_params);

            // rate monitoring to adjust the width of the sampling
            adaptive_width(i,n_check,arate,acc_counts,stds,agoal,j);
        }
        for (int k=0; k<n_params; ++k) {
            out << bulks_0[k] << "  ";
        }
        out << endl;
        cout << i+1 << " completed" << endl;
    }

    for (int i=0; i<nruns; ++i) {
        for (int j=0; j<n_params; ++j) {
            // get new proposed params
            flag = param_change(n_params,bulks_0,bulks_p,stds,mw,mp,j,inf_couplings);

            // get posterior
            prior = compute_prior(invcov,prior_means,bulks_p);
            lklp = compute_lkl(inf_couplings,CRUST,nrowscrust,flag);
            postp = prior*lklp;
            
            // metropolis hastings step
            post0 = metropolis(post0,postp,bulks_0,bulks_p,acc_counts,j,n_params);
        }
        for (int k=0; k<n_params; ++k) {
            aout << bulks_0[k] << "  ";
        }
        aout << endl;
        cout << i+1 << " completed" << endl;
    }

    dm1.cleanup(invcov,n_params);
    dm1.cleanup(CRUST,nrowscrust);
}

// MCMC for Finite Nuclei
// ################################################################################################
// ################################################################################################
// bulks has to be of form [BA, kf, mstar/m, K, J, L, Ksym, zeta, ms]

double compute_prior_FN(vector<double>& bulks, int n_params) {
    double chisq = 0.0;
    double prior_data[9] = {-16.0, 0.15, 0.6, 230.0, 34.0, 80.0, 0.0,   0.03, 500.0};
    double prior_unct[9] = {1.0  , 0.04, 0.1, 10.0 , 4.0 , 40.0, 500.0, 0.03, 50.0 };

    // initialize the temp arrays and get the prior distribution X^2 = (x-mu)^T Cov^-1 (x-mu)
    for (int i=0; i<n_params; ++i) {
        chisq = chisq + pow((bulks[i]-prior_data[i])/prior_unct[i],2.0);
    }
    return exp(-chisq/2.0);
}

// define the lkl function for a single nucleus
double compute_lkl_single(double exp_data[6],double BA_mev_th, double Rch_th, double FchFwk_th) {
    double lkl = exp(-0.5*pow(exp_data[0]-BA_mev_th,2.0)/pow(exp_data[1],2.0));
    if (exp_data[2] != -1) {
        lkl = lkl*exp(-0.5*pow(exp_data[2]-Rch_th,2.0)/pow(exp_data[3],2.0));
    }
    if (exp_data[4] != -1) {
        lkl = lkl*exp(-0.5*pow(exp_data[4]-FchFwk_th,2.0)/pow(exp_data[5],2.0));
    }
    return lkl;
}

// define the total liklihood
double compute_nuclei_v2(int num_nuclei, double params[16], int flag, double** exp_data) {
    int A[10] ={16,40,48,68,90,100,116,132,144,208};
    int Z[10] = {8,20,20,28,40,50,50,50,62,82};
    double Observables[7];
    double exp_data_single[6];
    double conv_help = 1.2;
    int exit_code = -1;
    int count = 0;

    double lkl = 1.0;
    if (flag == -1) {
        cout << "flagged" << endl;
        return 0.0;
    }

    for (int i=0; i<num_nuclei; ++i) {
        for (int j=0; j<6; ++j) {
            exp_data_single[j] = exp_data[i][j];
        }
        while (exit_code == -1) {
            exit_code = hartree_method(params,A[i],Z[i],20,gridsize,3,Observables,pow(1.1,count*1.0),false,false);
            count = count + 1;
            if (count > 7) {
                break;
            }
        }
        if (exit_code == -1) {
            cout << "flagged" << endl;
            return 0;
        }
        exit_code = -1;
        count = 0;

        lkl = lkl*compute_lkl_single(exp_data_single,Observables[0],Observables[3],Observables[5]);
    }
    return lkl;
}

int param_change_FN(int n_params, vector<double>& bulks_0, vector<double>& bulks_p, vector<double>& stds, int index, double fin_couplings[16]) {
    double masses[4] = {500.0,782.5,763.0,980.0};

    // copy old values
    for (int i=0; i<n_params; ++i) {
        bulks_p[i] = bulks_0[i];    
    }
    bulks_p[index] = rand_normal(bulks_0[index], stds[index]);

    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(bulks_p[1],3.0);
    masses[0] = bulks_p[8];
    int flag = bulkmc.get_parameters(bulks_p[0],p0,bulks_p[4],bulks_p[2]*mNuc_mev,bulks_p[3],bulks_p[5],bulks_p[6],bulks_p[7],0.0,0.0,0.0,0.0,masses,fin_couplings,true,1,true);
    if (flag == -1) {
        return -1;
    } else {
        return 0;
    }
}

void MCMC_FN(int nburnin, int nruns, string exp_file) {
    srand(time(0));
    int n_params = 9;
    double fin_couplings[16];
    double ms = 500.0; double mw = 782.5; double mp = 763.0; double md = 980.0;
    double masses[4] = {ms,mw,mp,md};
    double lklp, postp;
    int n_check = 50;
    double agoal = 0.3;
    ofstream out("MCMC_burnin.txt");
    ofstream aout("MCMC.txt");
    int flag;
    double** exp_data;
    dm1.importdata(exp_file,exp_data);

    // initialization
    vector<double> bulks_0 = {-16.2535,1.30974,0.601303,207.712,34.6836,82.3931,344.916,0.0360663,485.044};
    vector<double> bulks_p = {-16.2535,1.30974,0.601303,207.712,34.6836,82.3931,344.916,0.0360663,485.044};
    vector<double> stds = {0.01,0.003,0.004,1.8,1.0,10.1,30.0,0.0015,0.8};
    vector<int> acc_counts = {0,0,0,0,0,0,0,0,0};
    vector<double> arate = {0,0,0,0,0,0,0,0,0};
    
    // MCMC start
    double prior = compute_prior_FN(bulks_0,n_params);
    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(bulks_0[1],3.0);
    masses[0] = bulks_0[8];
    flag = bulkmc.get_parameters(bulks_0[0],p0,bulks_0[4],bulks_0[2]*mNuc_mev,bulks_0[3],bulks_0[5],bulks_0[6],bulks_0[7],0.0,0.0,0.0,0.0,masses,fin_couplings,true,1,true);
    double lkl0 = compute_nuclei_v2(10,fin_couplings,flag,exp_data);
    double post0 = prior*lkl0;


    for (int i=0; i<nburnin; ++i) {
        for (int j=0; j<n_params; ++j) {
            // get new proposed params
            flag = param_change_FN(n_params,bulks_0,bulks_p,stds,j,fin_couplings);
            cout << bulks_p[0] << "  " << bulks_p[1] << "  " << bulks_p[2] << "  " << bulks_p[3] << "  " << bulks_p[4] << "  " << bulks_p[5] << "  " << bulks_p[6] << "  " << bulks_p[7] << "  " << bulks_p[8] << endl;

            // get posterior
            prior = compute_prior_FN(bulks_p,n_params);
            lklp = compute_nuclei_v2(10,fin_couplings,flag,exp_data);
            postp = prior*lklp;
            cout << post0 << "  " << postp << endl;
            
            // metropolis hastings step
            post0 = metropolis(post0,postp,bulks_0,bulks_p,acc_counts,j,n_params);

            // rate monitoring to adjust the width of the sampling
            adaptive_width(i,n_check,arate,acc_counts,stds,agoal,j);
        }
        for (int k=0; k<n_params; ++k) {
            out << bulks_0[k] << "  ";
        }
        out << endl;
        cout << i+1 << " completed" << endl;
    }

    for (int i=0; i<nruns; ++i) {
        for (int j=0; j<n_params; ++j) {
            // get new proposed params
            flag = param_change_FN(n_params,bulks_0,bulks_p,stds,j,fin_couplings);

            // get posterior
            prior = compute_prior_FN(bulks_p,n_params);
            lklp = compute_nuclei_v2(10,fin_couplings,flag,exp_data);
            postp = prior*lklp;
            cout << post0 << "  " << postp << endl;
            
            // metropolis hastings step
            post0 = metropolis(post0,postp,bulks_0,bulks_p,acc_counts,j,n_params);
        }
        for (int k=0; k<n_params; ++k) {
            aout << bulks_0[k] << "  ";
        }
        aout << endl;
        cout << i+1 << " completed" << endl;
    }
    dm1.cleanup(exp_data,10);
}
