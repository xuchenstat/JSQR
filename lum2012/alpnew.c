#include "R.h"
#include "Rmath.h"           //Repackages functions from the 'nmath' package
#include "R_ext/Applic.h"
#include "unistd.h" // notice this! you need it!

//----------------------------Function prototypes-------------------------//
// (name, parameter, return type) let program know to expect these functions later

double log2(double x);
double *vect(int n);
int *ivect(int n);
double **mymatrix(int nr, int nc);
void Rprintvec(char *a, char *format, double *x, int n);
void Rprintmat(char *a, char *format, double **x, int m, int n, int flip);
void Rprintveci(char *a, char *format, int *x, int n);
double sumsquares(double *x, int n);
double inprod(double *x, double *y, int n);
double vmax(double *x, int n);
double vmin(double *x, int n);
double logsum(double *lx, int n);
double logmean(double *lx, int n);
double sum(double *x, int n);
int rdraw(int n, double *lprob, int inlog);
void locator_string(int *ix, int n, char *a);
void locator_string_inverse(char *a, int *ix);
void mmprod(double **a, double **b, double **c, int m, int k, int n, int atrans, int btrans, int ctrans);
void mvprod(double **a, double *b, double *c, int m, int k, int atrans);
void backsolve(int n, double **R, double **Rinv);
void mInprod(int n, int p, double **X, double **S, int transpose);
void set_lower_tri_zero(double **A, int n, int m);
void spchol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
            double gpcov(int, int, double*, int*, double **, int), double *covpar, int *include,
            double **K0, int addK, int max_scan, double *d, int dopivoting,
            int dostopping, int padzero, double *lpen, double *d2);
void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
          double *d, double **A, int dopivoting, int padzero, double eps);
void trisolve(double **R, int m, double *b, double *x, int transpose);
void triprod(double **R, int m, int n, double *x, double *b, int transpose);
double wtsumsquares(double *x, double *wt, int n);
double logpost(int phi, double alpha, double *z, int n);
// new utility functions
double ***array3d(int d1, int d2, int d3);

//global integers
int n, p, nphi, niter, thin, nphi, updateAlpha;
double tau, **x, *y, a0, b0, sigsq_inv, sdlalpha, tau1, tau2;
double *Ltilde, *Ztilde, **Lphigrid, ***Gphigrid;

int phi, *pivot, *rank;
double alpha, eta, *xi, a, b, *d, *zsamp;

//----------------------------Helper Functions-----------------------------//
// called inside of other functions


void alp(double *tauVar, double *xVar, double *yVar, double *SVar, int *dim, int *phiVar, double *init, double *dhyper, int *updateAlphaVar,
            double *betasamp, double *etasamp, double *xisamp, int *phisamp, double *alphasamp, double *acptsamp){

	int i, j, k, l, reach = 0; //reach reused to move pointers along via reach++
	n = dim[reach++]; p = dim[reach++]; nphi = dim[reach++]; niter = dim[reach++]; thin = dim[reach++];

  tau = tauVar[0]; updateAlpha = updateAlphaVar[0]; phi = phiVar[0];

  reach = 0;
	x = mymatrix(n, p);
	for(j = 0; j < p; j++) for(i = 0; i < n; i++) x[i][j] = xVar[reach++];
	y = yVar;

  reach = 0;
  alpha = init[reach++], eta = init[reach++];
  xi = vect(n);
  for(i = 0; i < n; i++) xi[i] = init[reach++];

  reach = 0;
  sigsq_inv = 1.0/dhyper[reach++]; a0 = dhyper[reach++]; b0 = dhyper[reach++]; sdlalpha = dhyper[reach++];

  reach = 0;
	Lphigrid = mymatrix(nphi, n);
	Gphigrid = array3d(nphi, n, n);

	for(i = 0; i < nphi; i++) {
		for(l = 0; l < n; l++) Lphigrid[i][l] = SVar[reach++];
		for(k = 0; k < n; k++) for(l = 0; l < n; l++) Gphigrid[i][k][l] = SVar[reach++];
	}

  tau1 = (1.0 - 2.0*tau)/tau/(1.0 - tau);
  tau2 = 2.0/tau/(1.0 - tau);


  a = a0 + n*3.0/2.0;
  double **D_hinv = mymatrix(2, n);
  int Dpar = 0, Dparnew = 1;
  for(i = 0; i < n; i++) {
    D_hinv[Dpar][i] = sqrt(1.0/tau2/xi[i]);
    D_hinv[Dparnew][i] = D_hinv[Dpar][i];
  }

  for(i = 0; i < (n + 1); i++) acptsamp[i] = 0.0;

  pivot = ivect(p); rank = ivect(1); d = vect(p); zsamp = vect(p); Ltilde = vect(n); Ztilde = vect(n);
  double alpha_new, cod, cod_new, lpr_cur, lpr_new, lp_cur, lp_new, lp_diff, ll_cur, ll_new, *beta = vect(p), *z = vect(n), *loglik = vect(nphi), *mu = vect(p), *xi_new = vect(n), *xb = vect(n), **u = mymatrix(2, n);
  double *y_res1 = vect(n), *y_res2 = vect(n), *Sigma_hz = vect(p), *Lambda_h = vect(n), **L_hinv_Gt = mymatrix(n, n), ***CD_hinv = array3d(2, n, n), **CD_hinvx = mymatrix(n, p);
  double **Sigma_inv = mymatrix(p, p), **Sigma_hinv = mymatrix(p, p), **Sigma_h = mymatrix(p, p), **Sigma = mymatrix(p, p);
  double *CD_hinvu = vect(n), *CD_hinvu_new = vect(n), *CD_hinvy = vect(n), *xCD_invy = vect(p);
  int iter, store_beta = 0, store_phi = 0, store_eta = 0, store_xi = 0, store_alpha = 0;
  int CDpar = 0, CDparnew = 1, upar = 0, uparnew = 1;
  lpr_cur = log(alpha) + 6.0*log1p(-alpha);

  // Rprintf("n = %d, p = %d, iter = %d, niter = %d, thin = %d", n, p, iter, niter, thin);
  GetRNGstate();
  for(iter = 0; iter < niter; iter++){
    for(i = 0; i < n; i++) Lambda_h[i] = sqrt(alpha * Lphigrid[phi][i] + 1.0 - alpha);
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++){
        L_hinv_Gt[i][j] = Gphigrid[phi][i][j]/Lambda_h[i];
        CD_hinv[CDpar][i][j] = L_hinv_Gt[i][j]*D_hinv[Dpar][j];
        CD_hinv[CDparnew][i][j] = CD_hinv[CDpar][i][j];
      }
    }
    mmprod(CD_hinv[CDpar], x, CD_hinvx, n, n, p, 0, 0, 0);

    for(i = 0; i < p; i++){
      for(Sigma_inv[i][i] = 0.0, k = 0; k < n; k++){
        Sigma_inv[i][i] += CD_hinvx[k][i] * CD_hinvx[k][i];
      }
      Sigma_inv[i][i] = eta*Sigma_inv[i][i] + sigsq_inv;
      for(j = i + 1; j < p; j++){
        for(Sigma_inv[i][j] = 0.0, k = 0; k < n; k++){
          Sigma_inv[i][j] += CD_hinvx[k][i] * CD_hinvx[k][j];
        }
        Sigma_inv[i][j] = eta*Sigma_inv[i][j];
        Sigma_inv[j][i] = Sigma_inv[i][j];
      }
    }

    chol(Sigma_hinv, p, 0.0, pivot, rank, p, d, Sigma_inv, 0, 0, 1e-10);
    backsolve(p, Sigma_hinv, Sigma_h);
    mInprod(p, p, Sigma_h, Sigma, 1);

    for(i = 0; i < n; i++) y_res1[i] = y[i] - tau1*xi[i];
    mvprod(CD_hinv[CDpar], y_res1, CD_hinvy, n, n, 0);
    mvprod(CD_hinvx, CD_hinvy, xCD_invy, p, n, 1);
    mvprod(Sigma, xCD_invy, mu, p, p, 0);
    for(i = 0; i < p; i++) mu[i] = eta*mu[i];
    for(i = 0; i < p; i++) zsamp[i] = rnorm(0.0, 1.0);
    triprod(Sigma_h, p, p, zsamp, Sigma_hz, 0);
    for(i = 0; i < p; i++) beta[i] = mu[i] + Sigma_hz[i];

    mvprod(x, beta, xb, n, p, 0);
    for(i = 0; i < n; i++) {
      u[upar][i] = y_res1[i] - xb[i];
      u[uparnew][i] = u[upar][i];
    }
    mvprod(CD_hinv[CDpar], u[upar], CD_hinvu, n, n, 0);
    cod = sumsquares(CD_hinvu, n)/2.0;
    b = b0 + cod + sum(xi, n);
    eta = rgamma(a, 1.0/b);

    for(i = 0; i < n; i++) {
      y_res2[i] = y[i] - xb[i];
      xi_new[i] = rexp(1.0/eta);
      u[uparnew][i] = y_res2[i] - tau1*xi_new[i];
      D_hinv[Dparnew][i] = sqrt(1.0/tau2/xi_new[i]);
      for(j = 0; j < n; j++) CD_hinv[CDparnew][j][i] = L_hinv_Gt[j][i]*D_hinv[Dparnew][i];
      mvprod(CD_hinv[CDparnew], u[uparnew], CD_hinvu_new, n, n, 0);
      cod_new = sumsquares(CD_hinvu_new, n)/2.0;

      lp_cur = -log(xi[i])/2.0 - eta*cod;
      lp_new = -log(xi_new[i])/2.0 - eta*cod_new;
      lp_diff = lp_new - lp_cur;
      if(runif(0.0, 1.0) < exp(lp_diff)){
        acptsamp[i] = acptsamp[i] + 1.0;
        u[upar][i] = u[uparnew][i];
        D_hinv[Dpar][i] = D_hinv[Dparnew][i];
        for(j = 0; j < n; j++) CD_hinv[CDpar][j][i] = CD_hinv[CDparnew][j][i];
        xi[i] = xi_new[i];
        cod = cod_new;
      } else {
        u[uparnew][i] = u[upar][i];
        D_hinv[Dparnew][i] = D_hinv[Dpar][i];
        for(j = 0; j < n; j++) CD_hinv[CDparnew][j][i] = CD_hinv[CDpar][j][i];
      }
    }

    for(i = 0; i < n; i++) z[i] = u[upar][i]/sqrt(xi[i]*tau2/eta);
    for(i = 0; i < nphi; i++) loglik[i] = logpost(i, alpha, z, n);
    phi = rdraw(nphi, loglik, 1);
    ll_cur = loglik[phi];

    if(updateAlpha){
      zsamp[0] = rnorm(qlogis(alpha, 0.0, 1.0, 1, 0), sdlalpha);
      alpha_new = plogis(zsamp[0], 0.0, 1.0, 1, 0);
      ll_new = logpost(phi, alpha_new, z, n);
      lpr_new = log(alpha_new) + 6.0*log1p(-alpha_new);

      if(runif(0.0, 1.0) < exp(ll_new - ll_cur + lpr_new - lpr_cur)){
        acptsamp[n] = acptsamp[n] + 1.0;
        alpha = alpha_new;
        lpr_cur = lpr_new;
      }
    }

    if((iter + 1) % thin == 0){
      for(i = 0; i < p; i++) betasamp[store_beta++] = beta[i];
      etasamp[store_eta++] = eta;
      for(i = 0; i < n; i++) xisamp[store_xi++] = xi[i];
      phisamp[store_phi++] = phi;
      alphasamp[store_alpha++] = alpha;
    }

    if((iter + 1) % 1000 == 0) Rprintf("tau = %g, iter = %d\n", tau, iter + 1);
    // if(iter > 1990) {Rprintf("iter = %d, beta = %g\n", iter, beta[1]); sleep(1);} // format is sleep(x); where x is # of seconds.
  }
  PutRNGstate();
  for(i = 0; i < (n + 1); i++) acptsamp[i] /= niter;
}


double logpost(int phi, double alpha, double *z, int n){
  double result = 0.0;
  int i;
  for(i = 0; i < n; i++) {
    Ltilde[i] = alpha * Lphigrid[phi][i] + 1.0 - alpha;
    result -= log(Ltilde[i]);
  }

  mvprod(Gphigrid[phi], z, Ztilde, n, n, 0);
  result -= wtsumsquares(Ztilde, Ltilde, n);
  return result/2.0;
}

//----------------------------------------------------------------------------//
//--------Functions for commonly used routines "mydefs" by Surya Tokdar-------//
//----------------------------------------------------------------------------//


double wtsumsquares(double *x, double *wt, int n){
	double wss = 0.0;
	int i;

	for(i = 0; i < n; i++) wss += x[i] * x[i] / wt[i];
	return wss;
}

// log base two
double log2(double x){
	return log(x)/log(2.0);
}

// Functions to create pointers and allocating memory to objects
// releases memory at end of call to .c

// Makes a pointer and allocates memory for vector of doubles
double * vect(int n){
	return (double *)R_alloc(n, sizeof(double));
}

// Makes a pointer and allocates memory for vector of integers
int * ivect(int n){
	return (int *)R_alloc(n, sizeof(int));
}

// Makes set of pointers and allocates memory for set of vectors of doubles (ie a matrix)
double ** mymatrix(int nr, int nc){
	int i;
	double **m;
	m = (double **) R_alloc(nr, sizeof(double *));
	for (i = 0; i < nr; i++)
		m[i] = (double *) R_alloc(nc, sizeof(double));
	return m;
}


// Prints a vector on screen with starting *a as character tag.
void Rprintvec(char *a, char *format, double *x, int n){
	int i;
	Rprintf("%s", a);
	for(i = 0; i < n; i++)
		Rprintf(format, x[i]);
	Rprintf("\n");
}

// Prints a matrix on screen with starting *a as character tag.
void Rprintmat(char *a, char *format, double **x, int m, int n, int flip){
	int i, j;
	Rprintf("%s\n", a);
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++)
			Rprintf(format, x[i][j]);
		Rprintf("\n");
	}
}


// Prints a vector of integers on screen with starting *a as character tag.
void Rprintveci(char *a, char *format, int *x, int n){
	int i;
	Rprintf("%s", a);
	for(i = 0; i < n; i++)
		Rprintf(format, x[i]);
	Rprintf("\n");
}


// A function to calculate the sum of squares of a vector x
double sumsquares(double *x, int n){
	double ss = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ss += x[i] * x[i];
	return ss;
}

// A function to calculate the inner product (sum of elementwise products) of x and y
double inprod(double *x, double *y, int n){
	double ip = 0.0;
	int i;
	for(i = 0; i < n; i++)
		ip += x[i] * y[i];
	return ip;
}


// Uses the inverse CDF method to generate a trunc-norm realization
// Note: not used in current code
// Note: pnorm gives the normal CDF; comes from Rmath package
// pnorm(x, mean, standard deviation, lower tail 1==TRUE  , log scale => 0 is false)
double rnormtrunc(double mu, double sigma, double lo, double hi){
	double u = runif(0.0, 1.0);
	double p = u * pnorm(hi, mu, sigma, 1, 0) + (1.0 - u) * pnorm(lo, mu, sigma, 1, 0);
	if(p <= 0.0) p = 1.0e-10;
	if(p >= 1.0) p = 1.0 - 1.0e-10;
	return qnorm(p, mu, sigma, 1, 0);
}


// Function to find the max value of vector
double vmax(double *x, int n){
	int i;
	double xmax = x[0];
	for(i = 1; i < n; i++) if(x[i] > xmax) xmax = x[i];
	return xmax;
}

// Numerically stable way to get sum of logs (log of product)
// lx being passed in is alredy log of something
double logsum(double *lx, int n){
	double lxmax = vmax(lx, n), a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += exp(lx[i] - lxmax);
	return lxmax + log(a);
}

// Mean version of above function
double logmean(double *lx, int n){
	return logsum(lx, n) - log((double)n);
}

// sum a vector
double sum(double *x, int n){
	double a = 0.0;
	int i;
	for(i = 0; i < n; i++) a += x[i];
	return a;
}

// Random draw from an atomic distribution (n=0....) with known/given probabilities prob
// Note: not used
int rdraw(int n, double *prob, int inlog){
	double psum, u = runif(0.0, 1.0), cprob;
	int j = 0;

	if(inlog) {
		psum = logsum(prob, n);
		cprob = exp(prob[0] - psum);
		while(u > cprob && j < n - 1) {
			j++;
			cprob += exp(prob[j] - psum);
		}
	} else {
		psum = sum(prob, n);
		cprob = prob[0] / psum;
		while(u > cprob && j < n - 1) {
			j++;
			if(prob[j] > 0.0) cprob += prob[j] / psum;
		}
	}
	return j;
}


// Unused function
void locator_string(int *ix, int n, char *a){
	const char *fmt[2]; fmt[0] = "%d"; fmt[1] = ".%d";
	int i, skip = 0;
	for(i = 0; i < n; i++) {
		if(ix[i]) {
			sprintf(a + skip, fmt[skip > 0], i + 1);
			skip = strlen(a);
		}
	}
}


// Unused function
void locator_string_inverse(char *a, int *ix){
	const char s[2] = ".";
	char *token;
	token = strtok(a, s);
	while( token != NULL ) {
		ix[atoi(token) - 1] = 1;
		token = strtok(a, s);
	}
}


// Matrix-matrix product.  a is the matrix, b is a matrix, c the matrix output,
// A is m x k in dimension which means b should be k x n; however there appear to be
// ways to take the transpose of A, B, or C with the logicals atrans, btrans, ctrans
void mmprod(double **a, double **b, double **c, int m, int k, int n, int atrans, int btrans, int ctrans){
	int i, j, l;
	if(!ctrans) {
		if(atrans && btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[l][i] * b[j][l];
		} else if (!atrans && btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[i][l] * b[j][l];
		} else if (atrans && !btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[l][i] * b[l][j];
		} else {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[i][j] += a[i][l] * b[l][j];
		}
	} else {
		if(atrans && btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[l][i] * b[j][l];
		} else if (!atrans && btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[i][l] * b[j][l];
		} else if (atrans && !btrans) {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[l][i] * b[l][j];
		} else {
			for(i = 0; i < m; i++)
				for(j = 0; j < n; j++)
					for(c[i][j] = 0.0, l = 0; l < k; l++) c[j][i] += a[i][l] * b[l][j];
		}
	}
}

// Matrix-vector product.  a is the matrix, b the product, c the vector output,
// A is m x k in dimension.  atrans tells whether to do A'b instead of Ab
void mvprod(double **a, double *b, double *c, int m, int k, int atrans){
	int i, l;
	if(atrans) {
		for(i = 0; i < m; i++)
			for(c[i] = 0.0, l = 0; l < k; l++) c[i] += a[l][i] * b[l];
	} else {
		for(i = 0; i < m; i++)
			for(c[i] = 0.0, l = 0; l < k; l++) c[i] += a[i][l] * b[l];
	}
}

// Function to find the min value of vector
double vmin(double *x, int n){
	int i;
	double xmin = x[0];
	for(i = 1; i < n; i++) if(x[i] < xmin) xmin = x[i];
	return xmin;
}

//--------Cholesky Decompositions------------//

// Function to set lower-left triangle of matrix (not necessarily square) = 0
void set_lower_tri_zero(double **A, int n, int m ){
	int i, j;
	for(i = 0; i < n; i++)
		for(j = i + 1; j < m; j++)
			A[j][i] = 0.0;
}


// Function to calculate the sparse cholesky factorization of matrix K0
void spchol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
            double gpcov(int, int, double*, int*, double**, int), double *covpar, int *include,
            double **K0, int addK, int max_scan, double *d, int dopivoting,
            int dostopping, int padzero, double *lpen, double *d2){


	// sparse cholesky factorization with pivoting and diagonal augmentation
	// accepts an empty matrix R and a function gpcov(i, j, ...) that is called
	// to compute the (i,j)-th element of the original covariance matrix


	set_lower_tri_zero(R, N, max_rank);

	int i, a, l;
	double u, b;

	if(dopivoting) {
		for(i = 0; i < N; i++)
			pivot[i] = i;
	}
	for(i = 0; i < N; i++) d[i] = gpcov(pivot[i], pivot[i], covpar, include, K0, addK);
	for(i = 0; i < max_scan; i++) d2[i] = lpen[pivot[i]];

	int k = 0, max_diag;
	for(max_diag = k, i = k + 1; i < max_scan; i++)
		if(d[i] > d[max_diag] * exp(d2[max_diag] - d2[i]))
			max_diag = i;
	tol *= d[max_diag];
	int flag = (k < max_rank);
	if(dostopping)
		flag = (d[max_diag] > tol);

	while(flag) {
		if(dopivoting) {
			if(max_diag > k) {
				a = pivot[k];
				pivot[k] = pivot[max_diag];
				pivot[max_diag] = a;

				b = d[k];
				d[k] = d[max_diag];
				d[max_diag] = b;

				b = d2[k];
				d2[k] = d2[max_diag];
				d2[max_diag] = b;

				for(i = 0; i < k; i++) {
					b = R[i][k];
					R[i][k] = R[i][max_diag];
					R[i][max_diag] = b;
				}
			}
		}

		R[k][k] = sqrt(d[k]);

		for(i = k + 1; i < N; i++) {
			u = gpcov(pivot[i], pivot[k], covpar, include, K0, addK);
			for(R[k][i] = u, l = 0; l < k; l++)
				R[k][i] -= R[l][i] * R[l][k];
			R[k][i] /= R[k][k];
			d[i] -= R[k][i] * R[k][i];
		}

		k++;
		flag = (k < max_rank);
		if(flag && dostopping) {
			for(max_diag = k, i = k + 1; i < max_scan; i++)
				if(d[i] > d[max_diag] * exp(d2[max_diag] - d2[i]))
					max_diag = i;
			flag = (d[max_diag] > tol);
		}
	}

	rank[0] = k;
	if(padzero) {
		for(l = k; l < N; l++)
			d[l] = 0.0;
	}
}

// Function to calculate Cholesky factorization of matrix A
void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
          double *d, double **A, int dopivoting, int padzero, double eps){

	set_lower_tri_zero(R, N, max_rank);

	int i, a, l;
	double u, b;

	for(i = 0; i < N; i++) {
		pivot[i] = i;
		d[i] = A[i][i] + eps * (1.0 + A[i][i]);
	}

	int k = 0, max_diag;
	for(max_diag = k, i = k + 1; i < N; i++)
		if(d[i] > d[max_diag])
			max_diag = i;
	int flag = (d[max_diag] > tol);


	while(flag) {
		if(dopivoting) {
			if(max_diag > k) {
				a = pivot[k];
				pivot[k] = pivot[max_diag];
				pivot[max_diag] = a;

				b = d[k];
				d[k] = d[max_diag];
				d[max_diag] = b;

				for(i = 0; i < k; i++) {
					b = R[i][k];
					R[i][k] = R[i][max_diag];
					R[i][max_diag] = b;
				}
			}
		}

		R[k][k] = sqrt(d[k]);

		for(i = k + 1; i < N; i++) {
			u = A[pivot[i]][pivot[k]];
			for(R[k][i] = u, l = 0; l < k; l++)
				R[k][i] -= R[l][i] * R[l][k];
			R[k][i] /= R[k][k];
			d[i] -= R[k][i] * R[k][i];
		}

		k++;
		flag = (k < max_rank);
		if(flag) {
			for(max_diag = k, i = k + 1; i < N; i++)
				if(d[i] > d[max_diag])
					max_diag = i;
			flag = (d[max_diag] > tol);
		}
	}

	rank[0] = k;
	if(padzero) {
		for(l = k; l < N; l++)
			d[l] = 0.0;
	}
}


// Function to solve an uppertriangular system of equations Rx=b
// R: upper triangular matrix
// m: size of the matrix
// b: right hand side
// transpose: option to solve R'x=b

void trisolve(double **R, int m, double *b, double *x, int transpose){

	int i, j;
	if(transpose) {
		for(j = 0; j < m; j++) {
			for(x[j] = b[j], i = 0; i < j; i++)
				x[j] -= x[i] * R[i][j];
			x[j] /= R[j][j];
		}
	} else {
		for(j = m - 1; j >= 0; j--) {
			for(x[j] = b[j], i = j + 1; i < m; i++)
				x[j] -= R[j][i] * x[i];
			x[j] /= R[j][j];
		}
	}
}

void triprod(double **R, int m, int n, double *x, double *b, int transpose){

	int i, j;
	if(transpose) {
		for(i = 0; i < m; i++)
			for(b[i] = 0.0, j = 0; j <= i; j++)
				b[i] += R[j][i] * x[j];
		for(; i < n; i++)
			for(b[i] = 0.0, j = 0; j < m; j++)
				b[i] += R[j][i] * x[j];
	} else{
		for(i = 0; i < m; i++)
			for(b[i] = 0.0, j = i; j < n; j++)
				b[i] += R[i][j] * x[j];
	}
}

// get the inverse of an upper triangular (n by n) matrix R using back substitution
void backsolve(int n, double **R, double **Rinv){
  int i, j, k;
  for (i = 0; i < n; i++) Rinv[i][i] = 1.0/R[i][i];
  for (j = 1; j < n; j++) {
    for (i = j - 1; i >= 0; i--){
      Rinv[j][i] = 0.0;
      for (Rinv[i][j] = 0.0, k = i + 1; k <= j; k++){
        Rinv[i][j] -= R[i][k]*Rinv[k][j];
      }
      Rinv[i][j] /= R[i][i];
    }
  }
}

// matrix inner product: S = X^\top X where X is a n by p matrix
// transpose = 1: S = XX^\top
void mInprod(int n, int p, double **X, double **S, int transpose){
  int i, j, k;
  if (!transpose){
    for (i = 0; i < p; i++){
      for (S[i][i] = 0.0, j = 0; j < n; j++) S[i][i] += X[j][i] * X[j][i];
      for (j = i + 1; j < p; j++){
        for (S[i][j] = 0.0, k = 0; k < n; k++){
          S[i][j] += X[k][i] * X[k][j];
        }
        S[j][i] = S[i][j];
      }
    }
  } else {
    for (i = 0; i < n; i++){
      for (S[i][i] = 0.0, j = 0; j < p; j++) S[i][i] += X[i][j] * X[i][j];
      for (j = i + 1; j < n; j++){
        for (S[i][j] = 0.0, k = 0; k < p; k++){
          S[i][j] += X[i][k] * X[j][k];
        }
        S[j][i] = S[i][j];
      }
    }
  }
}

double ***array3d(int d1, int d2, int d3){
	int i, j;
	double ***m = (double ***) R_alloc(d1, sizeof(double **));
	for (i = 0; i < d1; i++) {
		m[i] = (double **) R_alloc(d2, sizeof(double *));
		for (j = 0; j < d2; j++) {
			m[i][j] = (double *) R_alloc(d3, sizeof(double));
		}
	}
	return m;
}
