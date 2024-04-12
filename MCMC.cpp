#include "NumMethods.hpp"
#include "finitenuclei.hpp"
#include "infinitematter.hpp"
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
bulks bulkmc;
tools toolmc;

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

int excel_calibrate(double** init_DATA, string outname, int gd_sol_type, bool delta_coupling) {
    // Initialization
    int n_params = 16;  //(BA,kf,mstar,K,J,L,Ksym,zeta,xi,fw,fp,lambda_s,ms,mw,mp,md)
    double* index_array;
    int n_varied_params = counter(init_DATA, n_params, index_array,2);   // get the number of varied params and fill index array
    int code;
    double stds[n_varied_params];
    ofstream out(outname);
    double fin_couplings[n_params];
    double Observables[7];

    double base_params[n_params];
    double low_params[n_params];
    double high_params[n_params];
    for (int i=0; i<n_params; ++i) {
        base_params[i] = init_DATA[i][0];
    }

    int index;
    for (int i=0; i<n_varied_params; ++i) {
        index = index_array[i];
        stds[i] = init_DATA[index][1];
    }

    double masses[4];
    masses[0] = base_params[0];
    masses[1] = base_params[1];
    masses[2] = base_params[2];
    masses[3] = base_params[3];

    double BA = base_params[4]; 
    double kf = base_params[5]; 
    double mstar = base_params[9]*mNuc_mev;
    double K = base_params[10]; 
    double J_tilde = base_params[6]; 
    double L = base_params[7]; 
    double Ksym = base_params[8];
    double zeta = base_params[14];
    double xi = base_params[15];
    double fw = base_params[11];
    double fp = base_params[12];
    double lambda_s = base_params[13];
    cout << masses[0] << "  " << masses[1] << "  " << masses[2] << "  " << masses[3] << endl;
    cout << BA << "  " << kf << "  " << mstar << "  " << K << "  " << J_tilde << "  " << L << "  " << Ksym << "  " << zeta << "  " << fw << "  " << fp << "  " << lambda_s << "  " << xi << endl;

    double p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0); // get density
    bulkmc.get_parameters(BA,p0,J_tilde,mstar,K,L,Ksym,zeta,xi,lambda_s,fw,fp,masses,fin_couplings,true,gd_sol_type,delta_coupling);
    cout << " gd2 = " << fin_couplings[3] << endl;

    // print out the excel format and calculate several configurations for base
    out << masses[0] << "," << BA << "," << kf << "," << Ksym << "," << mstar/mNuc_mev << "," << K << "," << fw << "," << fp << "," << lambda_s << endl;

    code = hartree_method(fin_couplings,40,20,20,gridsize,2,Observables,1.2,false,false);
    if (code == 0) {
        out << "0," << Observables[0] << ",";
    } else {
        exit(0);
    }
    code = hartree_method(fin_couplings,48,20,20,gridsize,2,Observables,1.2,false,false);
    if (code == 0) {
        out << Observables[0] << "," << Observables[3] << "," << Observables[5] << ",";
    } else {
        exit(0);
    }
    code = hartree_method(fin_couplings,208,82,20,gridsize,2,Observables,1.3,false,false);
    if (code == 0) {
        out << Observables[0] << "," << Observables[3] << "," << Observables[5] << "," << fin_couplings[3] << "," << Observables[6] << endl;
    } else {
        exit(0);
    }

    // do the same for the low and high configurations
    for (int i=0; i<n_varied_params; ++i) {
        for (int j=0; j<n_params; ++j) {
            low_params[j] = base_params[j];
            high_params[j] = base_params[j];
        }


        index = index_array[i]; // get index of the parameter to change
        low_params[index] = base_params[index] - stds[i]; // get low limit of param
        high_params[index] = base_params[index] + stds[i];

        masses[0] = low_params[0];
        masses[1] = low_params[1];
        masses[2] = low_params[2];
        masses[3] = low_params[3];
        BA = low_params[4]; 
        kf = low_params[5]; 
        mstar = low_params[9]*mNuc_mev;
        K = low_params[10]; 
        J_tilde = low_params[6]; 
        L = low_params[7]; 
        Ksym = low_params[8];
        zeta = low_params[14];
        xi = low_params[15];
        fw = low_params[11];
        fp = low_params[12];
        lambda_s = low_params[13];

        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0); // get density
        bulkmc.get_parameters(BA,p0,J_tilde,mstar,K,L,Ksym,zeta,xi,lambda_s,fw,fp,masses,fin_couplings,true,gd_sol_type,delta_coupling);

        // print out the excel format and calculate several configurations for base
        code = hartree_method(fin_couplings,40,20,20,gridsize,2,Observables,1.2,false,false);
        if (code == 0) {
            out << low_params[index] << "," << Observables[0] << ",";
        } else {
            exit(0);
        }
        code = hartree_method(fin_couplings,48,20,20,gridsize,2,Observables,1.2,false,false);
        if (code == 0) {
            out << Observables[0] << "," << Observables[3] << "," << Observables[5] << ",";
        } else {
            exit(0);
        }
        code = hartree_method(fin_couplings,208,82,20,gridsize,2,Observables,1.3,false,false);
        if (code == 0) {
            out << Observables[0] << "," << Observables[3] << "," << Observables[5] << "," << fin_couplings[3] << "," << Observables[6] << endl;
        } else {
            exit(0);
        }

        masses[0] = high_params[0];
        masses[1] = high_params[1];
        masses[2] = high_params[2];
        masses[3] = high_params[3];
        BA = high_params[4]; 
        kf = high_params[5]; 
        mstar = high_params[9]*mNuc_mev;
        K = high_params[10]; 
        J_tilde = high_params[6]; 
        L = high_params[7]; 
        Ksym = high_params[8];
        zeta = high_params[14];
        xi = high_params[15];
        fw = high_params[11];
        fp = high_params[12];
        lambda_s = high_params[13];

        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0); // get density
        bulkmc.get_parameters(BA,p0,J_tilde,mstar,K,L,Ksym,zeta,xi,lambda_s,fw,fp,masses,fin_couplings,true,gd_sol_type,delta_coupling);

        // print out the excel format and calculate several configurations for base
        // print out the excel format and calculate several configurations for base
        code = hartree_method(fin_couplings,40,20,20,gridsize,2,Observables,1.2,false,false);
        if (code == 0) {
            out << high_params[index] << "," << Observables[0] << ",";
        } else {
            exit(0);
        }
        code = hartree_method(fin_couplings,48,20,20,gridsize,2,Observables,1.2,false,false);
        if (code == 0) {
            out << Observables[0] << "," << Observables[3] << "," << Observables[5] << ",";
        } else {
            exit(0);
        }
        code = hartree_method(fin_couplings,208,82,20,gridsize,2,Observables,1.3,false,false);
        if (code == 0) {
            out << Observables[0] << "," << Observables[3] << "," << Observables[5] << "," << fin_couplings[3] << "," << Observables[6] << endl;
        } else {
            exit(0);
        }
    }
    return 0;
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

    int exit_code = hartree_method(fin_couplings,A,Z,20,gridsize,3,Observables,1.2,false,true);
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

        exit_code = hartree_method(fin_couplings,A,Z,20,gridsize,3,Observables,1.2,false,true);
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

int generate_sample(double params[16], double** start_data) {
    double BA,p0,J,mstar,K,L,Ksym,zeta,xi,lambda_s,fw,fp;
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

double sample_param_space(int num_sets, string startfile) {
    double fin_couplings[16]; double Observables[7];
    bool physical_check;
    int n_params = 8;
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
        dm1.importdata_string("nspectrum_ref.txt",nref);
        dm1.importdata_string("pspectrum_ref.txt",pref);

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