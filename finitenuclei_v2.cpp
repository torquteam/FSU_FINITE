#include <cmath>
#include "NumMethods.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

data2 dmf;

// set constants
const double pi = 4.0*atan(1.0);
const double mNuc_mev = 939.0;
const double e0_e2permevfm = 0.055263; // e^2/(MeV fm) 
const double e0_unitless = e0_e2permevfm*r0_fm*enscale_mev;
const double mNuc_unitless = mNuc_mev/enscale_mev;

// set conversions
const double hbar_mevfm = 197.32698; // MeV fm
const double r0_fm = 1.25; //fm (arbitrary)
const double enscale_mev = pow(hbar_mevfm,2)/(2.0*mNuc_mev*pow(r0_fm,2)); // arbitrary
const double fm_to_inversemev = 1.0/hbar_mevfm;
const double conv_r0_en = r0_fm*fm_to_inversemev*enscale_mev;

// wave function class used to compute 
class wavefunction {
    private:

        // normalize the dirac wave functions
        void normalize(double** &F_unitless, double** &G_unitless, int npoints) {
            double psi2 = 0;
            double integrand,integrand_mid,integrand_next,F_mid_unitless,G_mid_unitless;

            // range and step size to integrate
            double r_init_unitless = F_unitless[0][0];
            double r_final_unitless = F_unitless[npoints-1][0];
            double h = (r_final_unitless-r_init_unitless)/(npoints-1);

            // integrate f^2 + g^2 
            for (int i=0; i<(npoints-1); ++i) {
                integrand = pow(F_unitless[i][1],2.0) + pow(G_unitless[i][1],2.0);
                F_mid_unitless = 0.5*(F_unitless[i][1] + F_unitless[i+1][1]);
                G_mid_unitless = 0.5*(G_unitless[i][1] + G_unitless[i+1][1]);
                integrand_mid = pow(F_mid_unitless,2.0) + pow(G_mid_unitless,2.0);
                integrand_next = pow(F_unitless[i+1][1],2.0) + pow(G_unitless[i+1][1],2.0);
                psi2 = psi2 + h/6.0*(integrand + 4.0*integrand_mid + integrand_next);
            }
            double norm = sqrt(1.0/psi2); // normalization factor 

            // normalize the wave functions
            for (int i=0; i<npoints; ++i) {
                F_unitless[i][1] = norm*F_unitless[i][1];
                G_unitless[i][1] = norm*G_unitless[i][1];
            }
        }

        // differential equation for upper component of dirac equation
        double dFdr(double r_unitless, double alpha, double En_unitless, double scalar_r, double vector_r, double F_r_unitless, double G_r_unitless) {
            double res = conv_r0_en*(En_unitless + mNuc_unitless - scalar_r - vector_r)*G_r_unitless + alpha/r_unitless*F_r_unitless;
            return res;
        }

        // differential equation for lower component of dirac equation
        double dGdr(double r_unitless, double alpha, double En_unitless, double scalar_r, double vector_r, double F_r_unitless, double G_r_unitless) {
            double res = -conv_r0_en*(En_unitless - mNuc_unitless + scalar_r - vector_r)*F_r_unitless - alpha/r_unitless*G_r_unitless;
            return res;
        }

        // rk4 method from left to right
        void rk4_lr(double r_init_unitless, double r_final_unitless, int npoints, double** scalar, double** vector, double** &F_unitless, double** &G_unitless, double en_unitless, double alpha, double FG_unitless[2]) {
            double k1, k2, k3, k4, l1, l2, l3, l4;
            double h = (r_final_unitless-r_init_unitless)/(npoints-1);
            double r_unitless = r_init_unitless;
            double scalar_mid, vector_mid;
            cout << "check if stepsizes agree: " << h << endl;

            // set intial condiitons
            double mstar = mNuc_unitless - scalar[0][1];
            double estar = en_unitless - vector[0][1];
            if (alpha>0) {
                F_unitless[0][1] = 1e-10;
                G_unitless[0][1] = F_unitless[0][1]*r_unitless*(mstar-estar)/(2.0*alpha+1.0);
            } else {
                G_unitless[0][1] = 1e-10;
                F_unitless[0][1] = G_unitless[0][1]*r_unitless*(mstar+estar)/(1.0-2.0*alpha);
            }
                F_unitless[0][0] = r_unitless;
                G_unitless[0][0] = r_unitless;

            // solve the dirac equation using fourth order runge kutta
            for (int i=1; i<npoints; ++i) {
                k1 = h*dFdr(r_unitless, alpha, en_unitless, scalar[i-1][1], vector[i-1][1], F_unitless[i-1][1], G_unitless[i-1][1]);
                l1 = h*dGdr(r_unitless, alpha, en_unitless, scalar[i-1][1], vector[i-1][1], F_unitless[i-1][1], G_unitless[i-1][1]);
                scalar_mid = 0.5*(scalar[i-1][1] + scalar[i][1]);
                vector_mid = 0.5*(vector[i-1][1] + vector[i][1]);
                k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[i-1][1]+k1/2.0, G_unitless[i-1][1]+l1/2.0);
                l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[i-1][1]+k1/2.0, G_unitless[i-1][1]+l1/2.0);
                k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[i-1][1]+k2/2.0, G_unitless[i-1][1]+l2/2.0);
                l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[i-1][1]+k2/2.0, G_unitless[i-1][1]+l2/2.0);
                k4 = h*dFdr(r_unitless+h, alpha, en_unitless, scalar[i][1], vector[i][1], F_unitless[i-1][1]+k3, G_unitless[i-1][1]+l3);
                l4 = h*dGdr(r_unitless+h, alpha, en_unitless, scalar[i][1], vector[i][1], F_unitless[i-1][1]+k3, G_unitless[i-1][1]+l3);
                F_unitless[i][1] = F_unitless[i-1][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
                G_unitless[i][1] = G_unitless[i-1][1] + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
                r_unitless = r_unitless+h;
                F_unitless[i][0] = r_unitless;
                G_unitless[i][0] = r_unitless;

                // save the value of the wave functions at the matching point
                FG_unitless[0] = F_unitless[i][1];
                FG_unitless[1] = G_unitless[i][1];
            }
        }

        // rk4 method from right to left
        void rk4_rl(double r_init_unitless, double r_final_unitless, int npoints, double** scalar, double** vector, double** &F_unitless, double** &G_unitless, double en_unitless, double alpha, double FG_unitless[2]) {
            double k1, k2, k3, k4, l1, l2, l3, l4;
            double h = (r_final_unitless-r_init_unitless)/(npoints-1);
            double r_unitless = r_init_unitless;
            double scalar_mid, vector_mid;
            cout << "check if stepsizes agree: " << h << endl;

            // set intial condiitons
            double mstar = mNuc_unitless - scalar[npoints-1][1];
            double estar = en_unitless - vector[npoints-1][1];
            double w = sqrt(mstar*mstar-estar*estar);
            F_unitless[0][1] = 1e-20;
            G_unitless[0][1] = -F_unitless[0][1]*(conv_r0_en*w + alpha/r_unitless)/(estar+mstar);
            F_unitless[0][0] = r_unitless;
            G_unitless[0][0] = r_unitless;

            // solve the dirac equation using fourth order runge kutta
            h = -h;
            for (int i=1; i<npoints; ++i) {
                k1 = h*dFdr(r_unitless, alpha, en_unitless, scalar[npoints-i][1], vector[npoints-i][1], F_unitless[npoints-i][1], G_unitless[npoints-i][1]);
                l1 = h*dGdr(r_unitless, alpha, en_unitless, scalar[npoints-i][1], vector[npoints-i][1], F_unitless[npoints-i][1], G_unitless[npoints-i][1]);
                scalar_mid = 0.5*(scalar[npoints-i][1] + scalar[npoints-i-1][1]);
                vector_mid = 0.5*(vector[npoints-i][1] + vector[npoints-i-1][1]);
                k2 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[npoints-i][1]+k1/2.0, G_unitless[npoints-i][1]+l1/2.0);
                l2 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[npoints-i][1]+k1/2.0, G_unitless[npoints-i][1]+l1/2.0);
                k3 = h*dFdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[npoints-i][1]+k2/2.0, G_unitless[npoints-i][1]+l2/2.0);
                l3 = h*dGdr(r_unitless+h/2.0, alpha, en_unitless, scalar_mid, vector_mid, F_unitless[npoints-i][1]+k2/2.0, G_unitless[npoints-i][1]+l2/2.0);
                k4 = h*dFdr(r_unitless+h, alpha, en_unitless, scalar[npoints-i-1][1], vector[npoints-i-1][1], F_unitless[npoints-i][1]+k3, G_unitless[npoints-i][1]+l3);
                l4 = h*dGdr(r_unitless+h, alpha, en_unitless, scalar[npoints-i-1][1], vector[npoints-i-1][1], F_unitless[npoints-i][1]+k3, G_unitless[npoints-i][1]+l3);
                F_unitless[npoints-i-1][1] = F_unitless[npoints-i][1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
                G_unitless[npoints-i-1][1] = G_unitless[npoints-i][1] + 1.0/6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
                r_unitless = r_unitless+h;
                F_unitless[npoints-i-1][0] = r_unitless;
                G_unitless[npoints-i-1][0] = r_unitless;

                // save the value of the wave functions at the matching point
                FG_unitless[0] = F_unitless[npoints-i-1][1];
                FG_unitless[1] = G_unitless[npoints-i-1][1];
            }
        }
        
        // solve the dirac equation for a given energy and return the determinant of the left and right solutions for the shooting method
        double DIRAC(double r_init_unitless, double r_final_unitless, int npoints, double alpha, double en_unitless, double** scalar, double** vector, double** &F_unitless, double** &G_unitless) {
            double k1, k2, k3, k4, l1, l2, l3, l4;
            double h = (r_final_unitless-r_init_unitless)/(npoints-1);
            double r_unitless = r_init_unitless;
            double FG_unitless[2];
            cout << "check if stepsizes agree: " << h << endl;
            
            // specify the matching point and grids for the left and right sides of the wavefunction
            int npoints_l = int(floor(npoints/4.0));
            int npoints_r = npoints - npoints_l + 1;
            double r_match_unitless = scalar[npoints_l-1][0];
            double** F_int_unitless; double** F_ext_unitless; double** G_int_unitless; double** G_ext_unitless;
            dmf.create(F_int_unitless,npoints_l,2);
            dmf.create(F_ext_unitless,npoints_r,2);
            dmf.create(G_int_unitless,npoints_l,2);
            dmf.create(G_ext_unitless,npoints_r,2);
            double F_int_unitless_rmatch; double G_int_unitless_rmatch; double F_ext_unitless_rmatch; double G_ext_unitless_rmatch;

            // shooting method by integrating both sides of the wave function to the matching point
            rk4_lr(r_init_unitless,r_match_unitless,npoints_l,scalar,vector,F_int_unitless,G_int_unitless,en_unitless,alpha,FG_unitless);
            F_int_unitless_rmatch = FG_unitless[0];
            G_int_unitless_rmatch = FG_unitless[1];
            rk4_rl(r_final_unitless,r_match_unitless,npoints_r,scalar,vector,F_ext_unitless,G_ext_unitless,en_unitless,alpha,FG_unitless);
            F_ext_unitless_rmatch = FG_unitless[0];
            G_ext_unitless_rmatch = FG_unitless[1];

            // stitch the solutions
            for (int i=0; i<npoints_l; ++i) {
                F_unitless[i][0] = F_int_unitless[i][0]; 
                F_unitless[i][1] = F_int_unitless[i][1];
                G_unitless[i][0] = G_int_unitless[i][0]; 
                G_unitless[i][1] = G_int_unitless[i][1];
            }
            for (int i=npoints_l; i<npoints; ++i) {
                F_unitless[i][0] = F_ext_unitless[npoints-i][0]; 
                F_unitless[i][1] = F_ext_unitless[npoints-i][1];
                G_unitless[i][0] = G_ext_unitless[npoints-i][0]; 
                G_unitless[i][1] = G_ext_unitless[npoints-i][1];
            }

            dmf.cleanup(F_int_unitless,npoints_l);
            dmf.cleanup(F_ext_unitless,npoints_r);
            dmf.cleanup(G_int_unitless,npoints_l);
            dmf.cleanup(G_ext_unitless,npoints_r);
            return -G_ext_unitless_rmatch*F_int_unitless_rmatch + F_ext_unitless_rmatch*G_int_unitless_rmatch;
        }

        // perform the shooting method to solve for the energy eigenvalue
        double SHOOTING(double alpha, double** scalar, double** vector, int npoints, double en_min_mev, double en_max_mev, double** &F_unitless, double** &G_unitless) {
            double r_init_unitless = scalar[0][0];
            double r_final_unitless = scalar[npoints-1][0];
            
            int count = 0;
            double midx_mev = 0.0;
            double sol_mev = 0.0;
            double en_det;
            double error = 1e-8; // check if this is low enough
            double midy = error*2;

            // bisect the given range until convergence or too many bisections
            while (fabs(midy)>error) {
                count = count + 1;
                midx_mev = (en_min_mev+en_max_mev)/2.0;
                en_det = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,midx_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
                midy = en_det;
                en_det = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,en_min_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);

                if ((midy*en_det)<0) {
                    en_max_mev = midx_mev;
                } else {
                    en_min_mev = midx_mev;
                }
                sol_mev = midx_mev;

                // Check for divergence
                if (count>100) {
                    cout << "SHOOTING: too many iterations bisection" << endl;;
                    exit(0);
                }
            }
            DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,sol_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
            normalize(F_unitless,G_unitless,npoints);
            return sol_mev;
        }

        double vector_dens(double j, double frac, double F_r_unitless, double G_r_unitless, double r_unitless) {
            double vdens= frac*(2.0*j+1.0)/(4.0*pi*pow(r_unitless,2.0))*(pow(F_r_unitless,2.0) + pow(G_r_unitless,2.0));
            return vdens;
        }

        double scalar_dens(double j, double frac, double F_r_unitless, double G_r_unitless, double r_unitless) {
            double sdens= frac*(2.0*j+1.0)/(4.0*pi*pow(r_unitless,2.0))*(pow(F_r_unitless,2.0) - pow(G_r_unitless,2.0));
            return sdens;
        }

        double tensor_dens(double j, double frac, double F_r_unitless, double G_r_unitless, double r_unitless) {
            double tdens= 2.0*frac*(2.0*j+1.0)/(4.0*pi*pow(r_unitless,2.0))*(F_r_unitless*G_r_unitless);
            return tdens;
        }

    public:
        // get all the bound states for the dirac equation
        void getBOUND_STATES(double** scalar, double** vector, int npoints) {
            double en_det, en_det_prev, bound_state_mev, en_mev;
            double alpha, j;
            double h = 4.0; // energy step
            int node = 0;   // set nodes to zero

            double r_init_unitless = scalar[0][0];
            double r_final_unitless = scalar[npoints-1][0];

            // create the nucleon arrays
            double** F_unitless; double** G_unitless;
            dmf.create(F_unitless,npoints,2);
            dmf.create(G_unitless,npoints,2);

            // print out the neutron energy spectrum (energy, node, l, j, alpha)
            ofstream nout("neutron_spectrum.txt");

            // initiate the while loop
            int nstates_l = 1;
            int l=0;    // start at l=0 state

            // find neutron states for different l until no more bound states are found
            while(nstates_l>0) {
                nstates_l = 0;

                // check the j = l+1/2 states for bound states
                j = l*1.0+0.5;
                alpha = l*1.0+0.5+0.5; // j=l+1/2 state
                en_mev = mNuc_mev-1e-8-h*20.0; // start at a low energy

                // find where the energy determinant is zero and start at the zeroth node
                node = 0;
                en_det_prev = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,en_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
                while(en_mev<mNuc_mev) {
                    en_det = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,en_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
                    if (en_det_prev*en_det < 0) {
                        bound_state_mev = SHOOTING(alpha,scalar,vector,npoints,en_mev-h,en_mev,F_unitless,G_unitless);
                        nout << fixed << setprecision(15) << bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                        node = node+1;
                    }
                    en_det_prev = en_det;
                    en_mev = en_mev + h;
                }

                nstates_l = node;   
                
                // check the j = l-1/2 states for bound states for l>0
                if (l!=0) {
                    j = l*1.0-0.5;
                    alpha = -l*1.0-0.5+0.5; // j=l-1/2 state
                    en_mev = mNuc_mev-1e-8-h*20.0;

                    node = 0;
                    en_det_prev = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,en_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
                    while(en_mev<mNuc_mev) {
                        en_det = DIRAC(r_init_unitless,r_final_unitless,npoints,alpha,en_mev/enscale_mev,scalar,vector,F_unitless,G_unitless);
                        if (en_det_prev*en_det < 0) {
                            bound_state_mev = SHOOTING(alpha,scalar,vector,npoints,en_mev-h,en_mev,F_unitless,G_unitless);
                            nout << fixed << setprecision(15)<< bound_state_mev << "  " << node << "  " << l << "  " << j << "  " << alpha << endl;
                            node = node+1;
                        }
                        en_det_prev = en_det;
                        en_mev = en_mev + h;
                    }
                    nstates_l = nstates_l + node;
                }
                l = l + 1;
            }

            dmf.cleanup(F_unitless,npoints);
            dmf.cleanup(G_unitless,npoints);
        }   

        // returns exit code -1 if no states are found, 0 if success
        int FILL_STATES(int N_nucleons, int en_col, int j_col, string filename) {
            double quantum_j, fill_frac;
            int nstates;
            int N = N_nucleons;

            // get the number of bound states
            int nrows = dmf.rowcount(filename);
            int ncols = dmf.colcount(filename);
            double** state_array;
            if (nrows == 0) {
                cout << "No states exist." << endl;
                return -1;
            }

            // import the bound states
            dmf.importdata(filename,state_array);
            dmf.sortasc(state_array,0,nrows,ncols);

            // fill the orbitals from lowest to highest energy
            ofstream nout(filename);
            int i=0;
            while (N_nucleons>0) {
                // if all the possible bound states are exhausted and there are more nucleons left then throw error
                if (i>(nrows-1)) {
                    cout << "Not enough bound states to fill all orbitals" << endl;
                    cout << "States filled: " << N_nucleons-N << "/" << N_nucleons << endl;
                    dmf.cleanup(state_array,nrows);
                    return -1;
                }
                quantum_j = state_array[i][j_col];
                nstates = 2*quantum_j+1;
                
                // if all nucleons fit in the orbital then add them else the orbital is fractionally filled
                if (nstates <= N) {
                    fill_frac = 1.0;
                    N = N - nstates;
                } else {
                    fill_frac = 1.0*N/nstates;
                    N = 0;
                }
                
                // copy the state information to new file now with only occupied orbitals
                for (int k=0; k<ncols; ++k) {
                    nout << fixed << setprecision(15) << state_array[i][k] << "  ";
                }
                nout << fixed << setprecision(15) << fill_frac << endl;

                i = i+1;
            }
            
            dmf.cleanup(state_array,nrows);
            return 0;
        }

        void WAVEFCNS(int npoints, double** &Fn_unitless, double** &Gn_unitless, double** &Fp_unitless, double** &Gp_unitless, double** scalar_n, double** scalar_p, double** vector_n, double** vector_p) {

            // initialize the wave function arrays
            int nstates_n = dmf.rowcount("neutron_spectrum.txt");
            int nstates_p = dmf.rowcount("proton_spectrum.txt");
            dmf.create(Fn_unitless,npoints,nstates_n+1);
            dmf.create(Gn_unitless,npoints,nstates_n+1);
            dmf.create(Fp_unitless,npoints,nstates_p+1);
            dmf.create(Gp_unitless,npoints,nstates_p+1);
            
            // fill arrays with the radial coordinates
            for (int i=0; i<npoints; ++i) {
                Fn_unitless[i][0] = scalar_n[i][0];
                Gn_unitless[i][0] = scalar_n[i][0];
                Fp_unitless[i][0] = scalar_n[i][0];
                Gp_unitless[i][0] = scalar_n[i][0];
            }

            double** Ftemp_array; double** Gtemp_array;
            double r_init = scalar_n[0][0];
            double r_final = scalar_n[npoints-1][0];
            double** nspectrum; double** pspectrum;
            dmf.importdata("neutron_spectrum.txt",nspectrum);
            dmf.importdata("proron_spectrum.txt",pspectrum);
            double alpha, en_mev;

            // solve the DIRAC equation for the eigen energies to get the wave functions (NEUTRON)
            for (int i=0; i<nstates_n; ++i) {
                alpha = nspectrum[i][4];
                en_mev = nspectrum[i][0];
                DIRAC(r_init,r_final,npoints,alpha,en_mev/enscale_mev,scalar_n,vector_n,Ftemp_array,Gtemp_array);
                normalize(Ftemp_array,Gtemp_array,npoints);

                // copy the solutions into the wf array
                for (int j=0; j<npoints; ++j) {
                    Fn_unitless[j][i+1] = Ftemp_array[j][1];
                    Gn_unitless[j][i+1] = Gtemp_array[j][1];
                }
            }

            // solve the DIRAC equation for the eigen energies to get the wave functions (PROTON)
            for (int i=0; i<nstates_p; ++i) {
                alpha = pspectrum[i][4];
                en_mev = pspectrum[i][0];
                DIRAC(r_init,r_final,npoints,alpha,en_mev/enscale_mev,scalar_p,vector_p,Ftemp_array,Gtemp_array);
                normalize(Ftemp_array,Gtemp_array,npoints);

                // copy the solutions into the wf array
                for (int j=0; j<npoints; ++j) {
                    Fp_unitless[j][i+1] = Ftemp_array[j][1];
                    Gp_unitless[j][i+1] = Gtemp_array[j][1];
                }
            }
            dmf.cleanup(nspectrum,nstates_n);
            dmf.cleanup(pspectrum,nstates_p);
            dmf.cleanup(Ftemp_array,npoints);
            dmf.cleanup(Gtemp_array,npoints);
        }

        void DENSITIES(double** densities_unitless, int npoints, double** Fn_unitless, double** Gn_unitless, double** Fp_unitless, double** Gp_unitless) {

            // initialie the density array and import the spectrum
            dmf.create(densities_unitless,npoints,9); //(r,sn,sp,vn,vp,tn,tp,c,w)
            dmf.zero(densities_unitless,npoints,9);
            double** nspectrum; double** pspectrum;
            dmf.importdata("neutron_spectrum.txt",nspectrum);
            dmf.importdata("proron_spectrum.txt",pspectrum);
            int nstates_n = dmf.rowcount("neutron_spectrum.txt");
            int nstates_p = dmf.rowcount("proton_spectrum.txt");

            // compute the densities
            double quantum_j, frac, r_unitless;
            for (int i=0; i<npoints; ++i) {
                r_unitless = Fn_unitless[i][0];
                densities_unitless[i][0] = r_unitless;
                for (int k=0; k<nstates_n; ++k) {
                    quantum_j = nspectrum[k][3];
                    frac = nspectrum[k][5];
                    densities_unitless[i][1] = densities_unitless[i][1] + scalar_dens(quantum_j,frac,Fn_unitless[i][k+1], Gn_unitless[i][k+1],r_unitless);
                    densities_unitless[i][3] = densities_unitless[i][3] + vector_dens(quantum_j,frac,Fn_unitless[i][k+1], Gn_unitless[i][k+1],r_unitless);
                    densities_unitless[i][5] = densities_unitless[i][5] + tensor_dens(quantum_j,frac,Fn_unitless[i][k+1], Gn_unitless[i][k+1],r_unitless);
                }
                for (int k=0; k<nstates_p; ++k) {
                    quantum_j = pspectrum[k][3];
                    frac = pspectrum[k][5];
                    densities_unitless[i][2] = densities_unitless[i][2] + scalar_dens(quantum_j,frac,Fp_unitless[i][k+1], Gp_unitless[i][k+1],r_unitless);
                    densities_unitless[i][4] = densities_unitless[i][4] + vector_dens(quantum_j,frac,Fp_unitless[i][k+1], Gp_unitless[i][k+1],r_unitless);
                    densities_unitless[i][6] = densities_unitless[i][6] + tensor_dens(quantum_j,frac,Fp_unitless[i][k+1], Gp_unitless[i][k+1],r_unitless);
                }
            }
        }
};

class meson_field {
    private:
        double greens_meson(double r_unitless, double rp_unitless, double meson_mass_unitless) {
            double res = 0;
            if (r_unitless>rp_unitless) {
                res = 1.0/meson_mass_unitless*rp_unitless/r_unitless*exp(-meson_mass_unitless*r_unitless*conv_r0_en)*sinh(meson_mass_unitless*rp_unitless*conv_r0_en);
            } else {
                res = 1.0/meson_mass_unitless*rp_unitless/r_unitless*exp(-meson_mass_unitless*rp_unitless*conv_r0_en)*sinh(meson_mass_unitless*r_unitless*conv_r0_en);
            }
            return res;
        }
    public:
        // V0 controls the depth, R controls the width, and a controls the smoothness of the tail
        void initialize(int npoints, double r_init_fm, double r_final_fm, int A, double couplings[4], double** &meson_array) {
            
            // set initial potentials (s,v,b,d)
            double V0_mev[4] = {42.0, 25.0, 0.0, 0.0};
            double R_fm = 1.2*pow(A,1.0/3.0);
            double a_fm = 0.65;

            // convert to unitless
            double R_unitless = R_fm/r0_fm;
            double a_unitless = a_fm/r0_fm;
            double r_init_unitless = r_init_fm/r0_fm;
            double r_final_unitless = r_final_fm/r0_fm;

            double h = (r_final_unitless-r_init_unitless)/(npoints-1);

            // initialize the grid
            double r_unitless = r_init_unitless;
            for (int i=0; i<npoints; ++i) {
                meson_array[i][0] = r_unitless;
                r_unitless = r_unitless+h;
            }

            // initialize each meson field
            for (int i=0; i<4; ++i) {
                double V_unitless = V0_mev[i]/enscale_mev;
                V_unitless = V_unitless*sqrt(couplings[i]);

                r_unitless = r_init_unitless;
                for (int j=0; j<npoints; ++j) {
                    meson_array[j][i+1] = V_unitless/(1.0+exp((r_unitless-R_unitless)/a_unitless));
                    r_unitless = r_unitless+h;
                }
            }
        }

        void split_meson(int npoints,double** meson_fields_unitless, double** &scalar_n, double** &scalar_p, double** &vector_n, double** &vector_p) {
            dmf.create(scalar_n,npoints,2);
            dmf.create(scalar_p,npoints,2);
            dmf.create(vector_n,npoints,2);
            dmf.create(vector_p,npoints,2);

            for (int i=0; i<npoints; ++i) {
                scalar_n[i][0] = meson_fields_unitless[i][0];
                scalar_p[i][0] = meson_fields_unitless[i][0];
                vector_n[i][0] = meson_fields_unitless[i][0];
                vector_p[i][0] = meson_fields_unitless[i][0];

                scalar_n[i][1] = meson_fields_unitless[i][1] - 0.5*meson_fields_unitless[i][4];
                scalar_p[i][1] = meson_fields_unitless[i][1] + 0.5*meson_fields_unitless[i][4];
                vector_n[i][1] = meson_fields_unitless[i][2] - 0.5*meson_fields_unitless[i][3];
                vector_p[i][1] = meson_fields_unitless[i][2] + 0.5*meson_fields_unitless[i][3] + meson_fields_unitless[i][5];
            }
        }

        void scalar() {

        }

        void vector() {
            
        }

};

class coulomb {
    public:
        void initialize(int npoints, double r_init_fm, double r_final_fm, int A, int Z, double** &meson_array) {
            double R_fm = 1.2*pow(A,1.0/3.0);
            double R_unitless = R_fm/r0_fm;   // convert unitless
            double r_init_unitless = r_init_fm/r0_fm;     // convert unitless
            double r_final_unitless = r_final_fm/r0_fm;   // convert unitless

            double h = (r_final_unitless-r_init_unitless)/(npoints-1);
            double r_unitless = r_init_unitless;

            // fill grid for coulomb field
            int i=0;
            while (r_unitless<R_unitless) {
                meson_array[i][5] = Z/(4.0*pi*e0_unitless)/(2.0*pow(R_unitless,3))*(3.0*pow(R_unitless,2) - pow(r_unitless,2));
                i = i+1;
                r_unitless = r_unitless+h;
            }
            while (i<npoints) {
                meson_array[i][5] = Z/(4.0*pi*e0_unitless)*1.0/(r_unitless);
                i = i+1;
                r_unitless = r_unitless+h;
            }
        }
};

int hartree(double fin_couplings[16], int A, int Z, int iterations, int npoints, double Observables[7], bool print_densities, bool print_meson_fields, bool print_wfs) {
    double r_init_fm = 1e-10;
    double r_final_fm = 1e-10;
    double yukawa_couplings[4] = {fin_couplings[0],fin_couplings[1],fin_couplings[2],fin_couplings[3]};

    // initialize the hartree with some meson fields (wood saxon like)
    double** meson_fields_unitless; double** scalar_n; double** scalar_p; double** vector_n; double** vector_p;
    dmf.create(meson_fields_unitless,npoints,6); // create a meson array of 6 columns (r,s,v,b,d,a)
    meson_field MESON;
    MESON.initialize(npoints,r_init_fm,r_final_fm,A,yukawa_couplings,meson_fields_unitless);
    MESON.split_meson(npoints,meson_fields_unitless,scalar_n,scalar_p,vector_n,vector_p);

    // get the energies
    wavefunction WF;
    WF.getBOUND_STATES(scalar_n,vector_n,npoints);
    WF.getBOUND_STATES(scalar_p,vector_p,npoints);
    WF.FILL_STATES(A-Z,0,3,"neutron_spectrum.txt");
    WF.FILL_STATES(Z,0,3,"proton_spectrum.txt");

    // get wavefunctions
    double** Fn_unitless; double** Gn_unitless;
    double** Fp_unitless; double** Gp_unitless;
    WF.WAVEFCNS(npoints,Fn_unitless,Gn_unitless,Fp_unitless,Gp_unitless,scalar_n,scalar_p,vector_n,vector_p);

    // start self consistent hartree and set the precision
    double precision = 1e-8;
    double err = precision*2.0;
    double** densities_unitless; // initialize density array of 9 columns (r,sn,sp,vn,vp,tn,tp,ch,wk)
    while (err<precision) {
        WF.DENSITIES(densities_unitless,npoints,Fn_unitless,Gn_unitless,Fp_unitless,Gp_unitless);
    }

}