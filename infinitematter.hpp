class tools {
private:
    double gss_FE(double kf, double gss, double gdd, double couplings[10], double t);
    double gdd_FE(double kf, double gss, double gdd, double couplings[10], double t);
    double gssFE_derivative(double kf, double gss, double gdd, double couplings[10], double t);
    double gddFE_derivative(double kf, double gss, double gdd, double couplings[10], double t);
public:
    double endensity_integral(double k, double mstar);
    double pr_integral(double k, double mstar);
    double scalardens(double k, double mstar);
    double common_integral(double k, double mstar);
    double scalarfield_2D_NR(double kf, double couplings[10], double t, double eps, double fields[2]);
    double FSUgssnewton(double dens, double gdd, double couplings[10], double t, double eps);
    double FSUgddnewton(double dens, double couplings[10], double t, double eps);
    double FSUgwwnewton(double couplings[10], double dens, double t);
    double get_gpomp2(double kf, double J, double L, double gss, double gww, double gsoms2, double gwomw2, double gdomd2, double kappa, double lambda, double zeta, double lambda_s);
    double get_gdomd2(double kf, double J, double L, double Ksym, double gss, double gww, double gsoms2, double gwomw2, double kappa, double lambda, double zeta, double lambda_s, int sol);
    double get_lambda_v(double kf, double J, double gss, double gww, double gdomd2, double gpomp2, double lambda_s);
    void convert_to_inf_couplings(double fin_couplings[16], double inf_couplings[10]);
    double get_t_betaeq(double kf, double couplings[10], double max, double fields_v[2]);
    double effK2(double kf, double couplings[10], double gss, double gww, double gpp, double gdd, double t);
    int ThermalCrust(double** crust, double** core, double** &neweos, int nrowscore, int nrowscrust, bool print, int nbcol, int prcol, int Keffcol);
};

class bulks {
public:
    double get_K(double kf, double gss, double gww, double couplings[10]);
    double get_J(double kf, double gss, double gww, double couplings[10]);
    double get_L(double kf, double gss, double gww, double couplings[10]);
    double get_Ksym(double kf, double gss, double gww, double couplings[10]);
    void get_bulkproperties(double couplings[10]);
    int get_parameters(double BA, double p0, double Jtilde, double mstar, double K, double L, double Ksym, double zeta, double xi, double lambda_s, double fw, double fp, double masses[4], double fin_couplings[16], bool flag, int gd_sol_type, bool delta_coupling);
};

class equationofstate {
public:
    double get_en(double kf, double t, double gss, double gww, double gpp, double gdd, double couplings[10]);
    double get_pr(double kf, double t, double gss, double gww, double gpp, double gdd, double couplings[10]);
    double qken(double cp, double mass, double deg);
    double qkpr(double cp, double mass, double deg);
    double qknb(double cp, double mass, double deg);
    int get_SNMEOS(double couplings[10], double** &eos, int npoints);
    int get_PNMEOS(double couplings[10], double** &eos, int npoints);
    int get_EOS_NSM(double couplings[10], double** &eos, int npoints, bool print, bool unstable);
    int get_SymmetryEnergy(double couplings[10], double** &eos, int npoints);
};