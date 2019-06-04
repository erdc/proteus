void Solve(double **a, double *b, int m, int n, double *solution, int MP, int NP)
{
int i;
double *w, **v, wmax, wmin;
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void dsvbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

w=dvector(1,NP);
v=dmatrix(1,NP,1,NP);

// Perform decomposition

dsvdcmp(a,m,n,w,v);

// Set up: see p65 of Press et al.

wmax = 0.;
for ( i=1 ; i<=n ; i++) if ( w[i] >= wmax ) wmax = w[i];
wmin = wmax * 1.e-12;
for ( i=1 ; i<=n ; i++) if ( w[i] <= wmin ) w[i] = 0.;

// Back substitute

dsvbksb(a,w,v,m,n,b,solution);

free_dmatrix(v,1,NP,1,NP);
free_dvector(w,1,NP);
}

