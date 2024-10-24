
using namespace std;

double rand_uniform(double cntr, double w);
double rand_normal(double mean, double stddev);
double sample_param_space(int num_sets, string startfile);
int generate_sample(double params[19], double** start_data);
int RBM_generate_fields(int A, int Z, string params_file);
void get_Observables(string param_set, int A, int Z);

int param_change(int n_params, vector<double>& bulks_0, vector<double>& bulks_p, vector<double>& stds, double mw, double mp, int index, double inf_couplings[10]);
double compute_prior(double** invcov, double means[7], vector<double>& bulks);
double compute_lkl(double inf_couplings[10], double** CRUST, int nrowscrust, int flag);
double metropolis(double lkl0, double lklp, vector<double>& bulks_0, vector<double>& bulks_p, vector<int>& acc_counts, int index, int n_params);
void adaptive_width(int iter, int n_check, vector<double>& arate, vector<int>& acc_counts, vector<double>& stds, double agoal, int index);
void MCMC_NS(int nburnin, int nruns, string covdata, string crust);
void RBM_error_check(string RBM_file, int n_params);
int MCMC_Observables(string MCMC_data, string crust);


double compute_prior_FN(vector<double>& bulks, int n_params);
double compute_lkl_single(double exp_data[6],double BA_mev_th, double Rch_th, double FchFwk_th);
double compute_nuclei_v2(int num_nuclei, double params[19], int flag, double** exp_data);
int param_change_FN(int n_params, vector<double>& bulks_0, vector<double>& bulks_p, vector<double>& stds, int index, double fin_couplings[19]);
void MCMC_FN(int nburnin, int nruns, string exp_file);