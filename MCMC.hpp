
using namespace std;

double rand_uniform(double cntr, double w);
double rand_normal(double mean, double stddev);
int excel_calibrate(double** init_DATA, string outname, int gd_sol_type, bool delta_coupling);
double sample_param_space(int num_sets, string startfile);
int generate_sample(double params[16], double** start_data);
int RBM_generate_fields(int A, int Z, string params_file);
void get_Observables(string param_set, int A, int Z);