void chkder ( int m, int n, double x[], double fvec[], double fjac[], 
  int ldfjac, double xp[], double fvecp[], int mode, double err[] );
void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
  double delta, double x[], double wa1[], double wa2[] );
double enorm ( int n, double x[] );
void fdjac1 ( void fcn ( int n, double x[], double f[], int &iflag, double p[] ), 
  int n, double x[], double fvec[], double fjac[], int ldfjac, int &iflag,
  int ml, int mu, double epsfcn, double wa1[], double wa2[], double p[] );
void fdjac2 ( void fcn ( int m, int n, double x[], double fvec[], int &iflag ),
  int m, int n, double x[], double fvec[], double fjac[], int ldfjac,
  int &iflag, double epsfcn, double wa[] );
int hybrd ( void fcn ( int n, double x[], double fvec[], int &iflag, double p[] ), 
  int n, double x[], double fvec[], double xtol, int maxfev, int ml, 
  int mu, double epsfcn, double diag[], int mode, double factor, int nprint, 
  int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[], 
  double wa1[], double wa2[], double wa3[], double wa4[], double p[] );
int hybrd1 ( void fcn ( int n, double x[], double fvec[], int &iflag, double p[] ), int n, 
  double x[], double fvec[], double tol, double wa[], int lwa, double p[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void qform ( int m, int n, double q[], int ldq );
void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
  int lipvt, double rdiag[], double acnorm[] );
void qrsolv ( int n, double r[], int ldr, int ipvt[], double diag[], 
  double qtb[], double x[], double sdiag[] );
void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] );
bool r1updt ( int m, int n, double s[], int ls, double u[], double v[], double w[] );
double r8_tiny ( );
double r8_uniform_01 ( int &seed );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
