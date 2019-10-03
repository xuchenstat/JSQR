#include "R.h"
#include "Rmath.h"           //Repackages functions from the 'nmath' package
#include "R_ext/Applic.h"
#include "time.h"
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
double rnormtrunc(double mu, double sigma, double lo, double hi);
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
void set_lower_tri_zero(double **A, int n, int m);
void spchol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
            double gpcov(int, int, double*, int*, double **, int), double *covpar, int *include,
            double **K0, int addK, int max_scan, double *d, int dopivoting,
            int dostopping, int padzero, double *lpen, double *d2);
void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
          double *d, double **A, int dopivoting, int padzero, double eps);
void trisolve(double **R, int m, double *b, double *x, int transpose);
void triprod(double **R, int m, int n, double *x, double *b, int transpose);
void adMCMCgp(int niter, int thin, int n, int nparmg, double *parmg, int phi, double rho, int nblocks, int **blocks, int *blocks_size, double **mu, double ***S, double theta, double eta, double *acpt_target, double acpt_hptarget, double decay, int refresh, int refreshhp, double *lm, double lmhp, double temp,
              int *refresh_counter, int refreshhp_counter, int verbose, int ticker,
              double *parmgsamp, double *acptsamp, double *lpsamp, double *lprsamp, double *lkappasamp, int *phisamp, double *acpthpsamp, double *lmglhsamp, double *usamp);
void transform_grid(double *w, double *v, int *ticks, double *dists);

double find_tau_lo(double target, double baseline, double a, double b, double taua, double taub);
double find_tau_up(double target, double baseline, double a, double b, double taua, double taub);
double part_trape_rp(double loc, double a, double b, double taua, double taub);

// new utility functions
double ***array3d(int d1, int d2, int d3);
double wtsumsquares(double *x, int phi, double rho, int n);
double lgauCopDen(int p, double *x, int phi, double rho);
double lprior(int p, int phi, double lkappa, double ltau);
// double lpost1(double sigma, double *y, int *cens, double lphi);
// double lpost2(double lphi);
// void Max_Search_Golden_Section1( double (*f)(double, double*, int*, double), double* a, double *fa, double* b, double* fb, double *y, int *cens, double lphi, double tolerance);
// void Max_Search_Golden_Section2( double (*f)(double), double* a, double *fa, double* b, double* fb, double tolerance);
// void INIT_COP(double *par, double *xVar, double *yVar, double *distVar, double *kappaVar, int *status, double *weights, int *toShrink, double *hyper, int *dim, double *gridpars, double *tauG, double *siglim,
//               int *covfun, int *distribution, double *lphiVar, double *decaypar);
//global integers
int n, p, L, mid, m, nkap, ngrid, shrink, dist, nphi;

//global constants
double *taugrid, *akap, *bkap, *lpkap, asig, bsig, arho, brho, shrinkFactor, ***Agrid, ***Rgrid, *ldRgrid, *lpgrid, **x, *y, *wt, **Lphigrid, ***Gphigrid;
int *cens, *zeta0_tick;

//global variables (i.e memory overwritten multiple times)
double *lb;
double **wgrid, *llgrid, *zknot, *lw;
double **wMat, **vMat, **vTilde, *w0, *zeta0dot, *zeta0, *vNormSq, **a, *aX, *gam, *xLin, *resLin, *b0dot, **bdot;
double *Q0Pos, **bPos, *Q0Neg, **bNeg;
double *zeta0_dist;


double *llvec, *pgvec, *rpvec, *lp, *zcop, *zprod;
int *pivot, *rank;

//----------------------------Helper Functions-----------------------------//
// called inside of other functions

// Function that returns e^(z/2); Allows sigma to be sampled on full real and returned to positive
double sigFn(double z){
	return exp(z/2.0);
}

// Function that returns 2*log(s); Return sigma from positive real back to full real
double sigFn_inv(double s) {
	return 2.0*log(s);
}

// Allows z, sampled on full real line, to be returned to values on [0.5,inf) when evaluation in dt needed
// 0.5 added to avoid over and underflow
double nuFn(double z)  {
	return 0.5 + 5.5*exp(z/2.0);
}

// Inverse nuFn; unused in C code.
double nuFn_inv(double nu) {
	return 2.0*log((nu - 0.5)/5.5);
}

// Two functions that only get used if shrink==TRUE
double dbase_joint_scl(double b, double *gam){
	double a = 0.0;
	int j;
	for(j = 0; j < p; j++) a += dt(gam[j] / b, 1.0, 1) - log(b);
	return a;
}
double dbasejoint(double *gam){
	int i;
	double pp = 0.525;
	for(i = 0; i < 10; i++) {
		lb[i] = dbase_joint_scl(qt(pp, 1.0, 1, 0), gam);
		pp += 0.05;
	}
	return logmean(lb, 10);
}

// A function to ensure that everything is between 0.000000000000001 and .999999999999999
double unitFn(double u){
	if(u < 1.0e-15)
		u = 1.0e-15;
	else if(u > 1.0 - 1.0e-15)
		u = 1.0 - 1.0e-15;
	return u;
}

//-----------------------------Base Functions-----------------------------//
// Functions related to the user-specified conditional distribution of Y|X=0

// Derivative of base quantile function Q0, equivalent to 1/f(Q(t)))
// Supports two cases: default is a t-distribution, scaled
// Note: qt(quantile, df, 1=lower tail, 0=not logged)
// gets pulled in through the Rmath.h package, which in turn calls the
// nmath.h package -- more info about functions such as dt and qt at
// https://svn.r-project.org/R/trunk/src/nmath/qt.c
double q0(double u, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = 1.0 / dlogis(qlogis(unitFn(u), 0.0, 1.0, 1, 0), 0.0, 1.0, 0);
		//val = 1.0 / dnorm(qnorm(unitFn(u), 0.0, 1.0, 1, 0), 0.0, 1.0, 0);
		break;
	case 3:
		val = 1.0 / dunif(qunif(u, -1.0, 1.0, 1, 0), -1.0, 1.0, 0);
		break;
	default:
		val = 1.0 / (dt(qt(unitFn(u), nu, 1, 0), nu, 0) * qt(0.9, nu, 1, 0));
		// Note: Scaling by qt helps disentangle tail from variance of distribution (depends on current nu)
		break;
	}
	return val;
	//return 1.0 / dnorm(qnorm(unitFn(u), 0.0, 1.0, 1, 0), 0.0, 1.0, 0);
	//return 1.0 / (unitFn(u) * unitFn(1.0 - u));
}

//Base quantile function
double Q0(double u, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = qlogis(unitFn(u), 0.0, 1.0, 1, 0);
		//val = qnorm(unitFn(u), 0.0, 1.0, 1, 0);
		break;
	case 3:
		val = qunif(u, -1.0, 1.0, 1, 0);
		break;
	default:
		val = qt(unitFn(u), nu, 1, 0) / qt(0.9, nu, 1, 0);
		break;
	}
	return val;
	//return qnorm(unitFn(u), 0.0, 1.0, 1, 0);
	//return qlogis(unitFn(u), 0.0, 1.0, 1, 0);
}

//Base CDF
double F0(double x, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = plogis(x, 0.0, 1.0, 1, 0);
		//val = pnorm(x, 0.0, 1.0, 1, 0);
		break;
	case 3:
		val = punif(x, -1.0, 1.0, 1, 0);
		break;
	default:
		val = pt(x * qt(0.9, nu, 1, 0), nu, 1, 0);
		break;
	}
	return val;
	//return pnorm(x, 0.0, 1.0, 1, 0);
	//return plogis(x, 0.0, 1.0, 1, 0);
}

// Adjustments to tails of base function's q0
// In tails of q0, replace with very heavy tailed t-distribution
double q0tail(double u, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = 1.0/dlogis(qlogis(u, 0.0, 1.0, 1, 0), 0.0, 1.0, 0);
		// val = 1.0/dnorm(qnorm(u, 0.0, 1.0, 1, 0), 0.0, 1.0, 0);
		break;
	case 3:
		val = 1.0/dt(qt(u, 0.5, 1, 0), 0.5, 0);
		break;
	default:
		val = 1.0/dt(qt(u, nu, 1, 0), nu, 0);
		break;
	}
	return val;
}


// Adjustments to F0 in the tail
double F0tail(double x, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = unitFn(plogis(x, 0.0, 1.0, 1, 0));
		// val = unitFn(pnorm(x, 0.0, 1.0, 1, 0));
		break;
	case 3:
		//val = punif(u, -1.0, 1.0, 1, 0);
		val = unitFn(pt(x, nu, 1, 0));
		break;
	default:
		val = unitFn(pt(x, nu, 1, 0));
		break;
	}
	return val;
}

// Adjustments to Q0 in the tail
double Q0tail(double u, double nu) {
	double val;
	switch (dist) {
	case 2:
		val = qlogis(u, 0.0, 1.0, 1, 0);
		// val = qnorm(u, 0.0, 1.0, 1, 0);
		break;
	case 3:
		val = qt(u, nu, 1, 0);
		break;
	default:
		//val = qt(u, 0.1, 1, 0);
		val = qt(u, nu, 1, 0);
		break;
	}
	return val;
}

// log likelihood of base density's tail-portion of distribution
double lf0tail(double x, double nu){
	double val;
	switch (dist) {
	case 2:
		val = dlogis(x, 0.0, 1.0, 1);
		//val = dnorm(x, 0.0, 1.0, 1);
		//val = dunif(x, -1.0, 1.0, 1);
		//val = dt(x, 0.1, 1);
		break;
	case 3:
		val = dt(x, nu, 1); // dt(value, df, 1=log)
		break;
	default:
		//val = dt(x, 0.1, 1);
		val = dt(x, nu, 1); // dt(value, df, 1=log)
		break;
	}
	return val;
}

//----------------------Additional Helper Functions-----------------------------//

// Function for approximating a definite integral using the trapezoid rule
// 1 x: input function (heights of curve to be integrated)
// 2 h: input grid locations
// 3 n: number of bins (between grid locations)
// 4 c: modified (or "returned") vector of cumulative summations of integral
//      after each additional trapezoid
void trape(double *x, double *h, int n, double *c, int reverse){
	int i, j = 0;
	c[0] = 0.0;
	if(reverse) {
		for(i = 1; i < n; i++) {
			c[i] = c[i-1] + 0.5 * (h[j] - h[j-1]) * (x[j] + x[j-1]);
			j--;
		}
	} else {
		for(i = 1; i < n; i++) {
			c[i] = c[i-1] + 0.5 * (h[j+1] - h[j]) * (x[j] + x[j+1]);
			j++;
		}
	}
}

// shrinkage function, not used in default settings
double shrinkFn(double x){
	return 1.0; ///(1.0 + log(x));
}

//-----------Functions for evaluating log likelihood or log posterior--------//

// Function for obtaining part of the log likelihood; this part specific to
// multivariate t-dist for covariate=intercept (a discretized t-process
// after integrating out kappa_0 [IG(0.1,0.1) scale mixture] from GP)
// Note: moves from low-rank GP approx to full set of L (all lambdas)
//       using global constants Agrid and Rgrid
// Arguments
// 1 wknot    : input, not modified
// 2 w        : created/modified/output (equivalent to w-tilde in paper)
// 3 postgrid : partial log posterior evaluation, vector across lambdas
double ppFn0(double *wknot, double *w, double *postgrid){
	int i, l;
	double akapm, zss;
	for(i = 0; i < ngrid; i++) { // i indexes lambdas
		mvprod(Agrid[i], wknot, wgrid[i], L, m, 0); // Agrid is Lxm constant matrix. Multiply Agrid by wknot to produce wgrid.
		trisolve(Rgrid[i], m, wknot, zknot, 1); // Rgrid constant. Solve for zknot (in place, creates vector solution) Rgrid'%*%zknot = wknot
		zss = sumsquares(zknot, m); // obtain sums of squares, save in zss
		//akapm = 1.5 + 0.5 * (double)m;
		//llgrid[i] = -akapm * log1p(0.5 * zss / 1.5);
		akapm = 0.1 + 0.5 * (double)m; // m is the number of knots. For intercept, akap is fixed at 0.1
		llgrid[i] = -akapm * log1p(0.5 * zss / 0.1); // log1p = ln(1+x)
		postgrid[i] = llgrid[i] - ldRgrid[i] + lpgrid[i]; // References objects made in previous steps (Agrid, Rgrid, ldRgrid AND lpgrid)
		                                                  // lpgrid holds log priors for lambdas
	}

	double lps = logsum(postgrid, ngrid); // Sum across G/ngrid items (over lambda grid)
	for(i = 0; i < ngrid; i++) postgrid[i] = exp(postgrid[i] - lps);
	for(l = 0; l < L; l++) for(w[l] = 0.0, i = 0; i < ngrid; i++) w[l] += wgrid[i][l] * postgrid[i];
	return lps;
}


// Another function for obtaining part of the log likelihood; this one similar
// to previous function but determines contributions of each covariate over grid
// of lambdas
// Note: kappa_j prior params not fixed at 0.1, 0.1; can vary for each covariate
double ppFn(double *wknot, double *w, double *postgrid){
	int i, j, l;
	double akapm, zss;
	for(i = 0; i < ngrid; i++) {
		mvprod(Agrid[i], wknot, wgrid[i], L, m, 0);
		trisolve(Rgrid[i], m, wknot, zknot, 1);
		zss = sumsquares(zknot, m);
		for(j = 0; j < nkap; j++) {
			akapm = akap[j] + 0.5 * (double)m;
			lw[j] = -akapm * log1p(0.5 * zss/ bkap[j]) + lgamma(akapm) - lgamma(akap[j]) - .5 * (double)m * log(bkap[j]) + lpkap[j];
		}
		llgrid[i] = logsum(lw, nkap);
		postgrid[i] = llgrid[i] - ldRgrid[i] + lpgrid[i];
	}

	double lps = logsum(postgrid, ngrid);
	for(i = 0; i < ngrid; i++) postgrid[i] = exp(postgrid[i] - lps);
	for(l = 0; l < L; l++) for(w[l] = 0.0, i = 0; i < ngrid; i++) w[l] += wgrid[i][l] * postgrid[i];
	return lps;
}


void logpostFn_noX(double *par, double temp, int llonly, double *ll, double *pg, double *y, int *cens, double *lp, double *rp){

	int i, l, reach = 0, reach2 = 0;
	double w0max, zeta0tot, lps0, gam0, sigma, nu, QPos, QPosold, QNeg, QNegold, sigmat1, sigmat2, den0;

	lps0 = ppFn0(par, w0, pg);
	reach += m;
	reach2 += ngrid;

	w0max = vmax(w0, L);
	for(l = 0; l < L; l++) zeta0dot[l] = exp(w0[l] - w0max);
	trape(zeta0dot + 1, taugrid + 1, L-2, zeta0 + 1, 0);
	zeta0tot = zeta0[L-2];
	zeta0[0] = 0.0; zeta0[L-1] = 1.0;
	for(l = 1; l < L-1; l++) zeta0[l] = taugrid[1] + (taugrid[L-2] - taugrid[1]) * zeta0[l] / zeta0tot;
	zeta0dot[0] = 0.0; zeta0dot[L-1] = 0.0;
	for(l = 1; l < L-1; l++) zeta0dot[l] = (taugrid[L-2] - taugrid[1]) * zeta0dot[l] / zeta0tot;

	if(temp > 0.0) {
		for(i = 0; i < n; i++) {
			ll[i] = log(0.0);
			rp[i] = taugrid[mid];
		}

		gam0 = par[reach++];
		sigma = sigFn(par[reach++]);
		nu = nuFn(par[reach++]);
		//Rprintf("sigma = %g, nu = %g\n", sigma, nu);

		for(i = 0; i < n; i++) resLin[i] = y[i] - gam0;
		//Rprintvec("resLin = ", "%g ", resLin, n);

		for(l = 0; l < L; l++) b0dot[l] = sigma * q0(zeta0[l], nu) * zeta0dot[l];

		trape(b0dot + mid, taugrid + mid, L - mid, Q0Pos, 0); Q0Pos[L-mid] = qt(1.0, 1.0, 1, 0);
		trape(b0dot + mid, taugrid + mid, mid + 1, Q0Neg, 1); Q0Neg[mid+1] = qt(1.0, 1.0, 1, 0);
		//Rprintvec("Q0Pos = ", "%g ", Q0Pos, L - mid + 1);
		//Rprintvec("Q0Neg = ", "%g ", Q0Neg, mid + 2);

		sigmat1 = sigmat2 = sigma;

		//double h=0.5;
		for(i = 0; i < n; i++) {
			if(resLin[i] == 0.0) {
				den0 = b0dot[mid];
				ll[i] = -log(den0);
			} else if(resLin[i] > 0.0) {
				l = 0;
				QPosold = 0.0;
				QPos = Q0Pos[l];
				while(resLin[i] > QPos && l < L-mid-1) {
					QPosold = QPos;
					l++;
					QPos = Q0Pos[l];
				}
				if(l == L - mid - 1) {
					rp[i] = F0tail(Q0tail(taugrid[L-2], nu) + (resLin[i] - QPosold)/sigmat1, nu);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = lf0tail(Q0tail(taugrid[L-2], nu) + (resLin[i] - QPosold)/sigmat1, nu) - log(sigmat1); break;
					}
				} else {
					rp[i] = find_tau_lo(resLin[i], QPosold, b0dot[mid+l-1], b0dot[mid+l], taugrid[mid+l-1], taugrid[mid+l]);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = -log(part_trape_rp(rp[i], b0dot[mid+l-1], b0dot[mid+l], taugrid[mid+l-1], taugrid[mid+l])); break;
					}
				}
			} else {
				l = 0;
				QNegold = 0.0;
				QNeg = Q0Neg[l];
				while(resLin[i] < -QNeg && l < mid) {
					QNegold = QNeg;
					l++;
					QNeg = Q0Neg[l];
				}
				if(l == mid) {
					rp[i] = F0tail(Q0tail(taugrid[1], nu) + (resLin[i] + QNegold)/sigmat2, nu);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = lf0tail(Q0tail(taugrid[1], nu) + (resLin[i] + QNegold)/sigmat2, nu) - log(sigmat2); break;
					}
				} else {
					rp[i] = find_tau_up(resLin[i], -QNegold, b0dot[mid-l], b0dot[mid-l+1], taugrid[mid-l], taugrid[mid-l+1]);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = -log(part_trape_rp(rp[i], b0dot[mid-l], b0dot[mid-l+1], taugrid[mid-l], taugrid[mid-l+1])); break;
					}
				}
			}
			if(ll[i] == qt(1.0, 1.0, 1, 0)) Rprintf("i = %d, ll[i] = %g, resLin[i] = %g, l = %d\n", i, ll[i], resLin[i], l);
		}
	} else {
		for(i = 0; i < n; i++) ll[i] = 0.0;
	}
	//Rprintvec("ll = ", "%g ", ll, n);
	lp[0] = temp * inprod(ll, wt, n);

	if(!llonly) {
		//lp[0] += lps0 + dt(par[m*(p+1)], 1.0, 1) + dlogis(par[(m+1)*(p+1)], 0.0, 1.0, 1) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		//lp[0] += lps0 + dlogis(par[(m+1)*(p+1)], 0.0, 1.0, 1) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		lp[0] += lps0 + dlogis(par[m+2], 0.0, 1.0, 1);
		//else lp[0] -= 0.5*sumsquares(par + m*(p+1) + 1, p)/9.0;
	}
}


// Function for evaluating the log posterior (case for log-likelihood only)
// Arguments/Inputs/Outputs
// 1 *par      Pointer to vector of all parameters
// 2 *temp
// 3 *llonly   1=evaluate ll only, leaving out nu prior, 0=evaluate log posterior including nu prior
// 4 *ll       Vector of log-likelihood evaluations for each of n observations
// 5 *pg       Vector of contributions of gamma priors to log likelihood
// 6 *rp       Vector of taus associated with each of n observations
// Output:
// double lp   Log posterior
void logpostFn(double *par, double temp, int llonly, double *ll, double *pg, double *y, int *cens, double *lp,  double *rp){

	int i, j, l, reach = 0, reach2 = 0;
	double w0max, zeta0tot, lps0, gam0, sigma, nu, QPos, QPosold, QNeg, QNegold, sigmat1, sigmat2, den0;

	lps0 = ppFn0(par, w0, pg); // evaluate log posterior for the intercept contribution... here the w0 gets updated
	reach += m; // increment over knots
	reach2 += ngrid; // increment over gamma grid

	w0max = vmax(w0, L);
	for(l = 0; l < L; l++) zeta0dot[l] = exp(shrinkFactor * (w0[l] - w0max)); //Subtract w0max for numerical stability
	trape(zeta0dot + 1, taugrid + 1, L-2, zeta0 + 1, 0); // Get 'functional' zeta0 by integrating its derivative zeta0dot along taugrid.
	zeta0tot = zeta0[L-2];        // Integral under full curve
	zeta0[0] = 0.0; zeta0[L-1] = 1.0; // Fix zeta at 0 to be 0 and zeta at 1 to be 1.
	for(l = 1; l < L-1; l++) zeta0[l] = taugrid[1] + (taugrid[L-2] - taugrid[1]) * zeta0[l] / zeta0tot; // move everything sozeta(tau)=tau on first and last delta intervals and its appropriately scaled self in between
	zeta0dot[0] = 0.0; zeta0dot[L-1] = 0.0; // implies that zeta0 is constant on these first and last delta intervals
	for(l = 1; l < L-1; l++) zeta0dot[l] = (taugrid[L-2] - taugrid[1]) * zeta0dot[l] / zeta0tot;

	// Setup for subsequnt linear interpolation to get vMat = zeta0(wMat)
	// creates zeta0_tick (vector of indices) and zeta0_dist (proportion of way between subsequent tau gridpoint)
	zeta0_tick[0] = 0; zeta0_dist[0] = 0.0;
	i = 1;
	for(l = 1; l < L-1; l++) {
		while(zeta0[l] >= taugrid[i] && i < L) i++;
		if(i > L-1) i = L-1;
		zeta0_tick[l] = i-1;
		zeta0_dist[l] = (zeta0[l] - taugrid[i-1]) / (taugrid[i] - taugrid[i-1]);
	}
	zeta0_tick[L-1] = L-2; zeta0_dist[L-1] = 1.0;
	//Rprintvec("zeta0 = ", "%g ", zeta0, L);
	//Rprintveci("zeta0_tick = ", "%d ", zeta0_tick, L);
	//Rprintvec("zeta0_dist = ", "%g ", zeta0_dist, L);

	for(j = 0; j < p; j++) { // add to log posterior the contributions from other covariates.  first m in par relate to intercept, remaining to other covariates.  vmat updates here
		lps0 += ppFn(par + reach, wMat[j], pg + reach2); // wMat comes from the ppFn function
		// Linear interpolation to get vMat = zeta0(wMat)
		transform_grid(wMat[j], vMat[j], zeta0_tick, zeta0_dist);
		reach += m; reach2 += ngrid;
	}

	if(temp > 0.0) {
		for(i = 0; i < n; i++) {
			ll[i] = log(0.0);
			rp[i] = taugrid[mid];
		}
		mmprod(vMat, x, a, L, p, n, 1, 1, 0); // vMat%*%x  (L x p0)%*%(p x n) gives 'a' (L x n) as output
		for(l = 0; l < L; l++) {
			for(vNormSq[l] = 0.0, j = 0; j < p; j++) vNormSq[l] += vMat[j][l] * vMat[j][l];
			if(vNormSq[l] > 0.0) {
				aX[l] = -vmin(a[l], n) / sqrt(vNormSq[l]);
				for(j = 0; j < p; j++) vTilde[j][l] = vMat[j][l] / (aX[l] * sqrt(1.0 + vNormSq[l])); // calculate v-tidle, the attenuator
			} else {
				for(j = 0; j < p; j++) vTilde[j][l] = 0.0;
			}
		}

// Median values, sigma, and nu are stored as last parameters is par.
		gam0 = par[reach++]; // pointer to gamma_0 (median intercept)
		gam = par + reach; reach += p; // pointer to gamma's (median effect for each covariate)
		sigma = sigFn(par[reach++]); // sigma, stored in logged form, exponentiated here
		nu = nuFn(par[reach++]); // nu, transformed
//Rprintf("sigma = %g, nu = %g\n", sigma, nu);

// Obtain prediction and residuals for each observation if predictions
// were solely based reference
		for(i = 0; i < n; i++) xLin[i] = gam0 + inprod(x[i], gam, p);
		for(i = 0; i < n; i++) resLin[i] = y[i] - xLin[i];
//Rprintvec("resLin = ", "%g ", resLin, n);

// Obtain beta_0(tau)'s derivative at each tau location
		// Then obtain beta(tau)'s derivative at each tau location
		for(l = 0; l < L; l++) b0dot[l] = sigma * q0(zeta0[l], nu) * zeta0dot[l];
		for(j = 0; j < p; j++) for(l = 0; l < L; l++) bdot[j][l] = b0dot[l] * vTilde[j][l];

// Obtain beta_j(tau) for j=0 to p, by integrating each derivative.
// Each integration dones in two parts cumulatively moving away from the median va\lues.
		trape(b0dot + mid, taugrid + mid, L - mid, Q0Pos, 0); Q0Pos[L-mid] = qt(1.0, 1.0, 1, 0);
		trape(b0dot + mid, taugrid + mid, mid + 1, Q0Neg, 1); Q0Neg[mid+1] = qt(1.0, 1.0, 1, 0);
//Rprintvec("Q0Pos = ", "%g ", Q0Pos, L - mid + 1);
//Rprintvec("Q0Neg = ", "%g ", Q0Neg, mid + 2);

		for(j = 0; j < p; j++) {
			trape(bdot[j] + mid, taugrid + mid, L - mid, bPos[j], 0);
			trape(bdot[j] + mid, taugrid + mid, mid + 1, bNeg[j], 1);
		}

//sigmat1 = sigma * q0(taugrid[L-2],nu) / q0tail(taugrid[L-2]);
//sigmat2 = sigma * q0(taugrid[1],nu) / q0tail(taugrid[1]);
		sigmat1 = sigmat2 = sigma;

// Contribution to log-likelihood of point y_i
		double qdot_lo = 0.0, qdot_up = 0.0; //, h=0.5;
		for(i = 0; i < n; i++) { // add in contributions from each observation
			if(resLin[i] == 0.0) { // case Y_i=median data quantile
				for(den0 = b0dot[mid], j = 0; j < p; j++) den0 += x[i][j] * bdot[j][mid];
				ll[i] = -log(den0);
			} else if(resLin[i] > 0.0) { // case Y_i > median data quantile,
				l = 0;
				QPosold = 0.0;
				// Find conditional quantile of median (on x_i) by adding contributions of all covariates
				for(QPos = Q0Pos[l], j = 0; j < p; j++) QPos += x[i][j] * bPos[j][l];
				// Then check where Y_i is in relation to quantile.
				// If it is still above the quantile, calc Conditional quantile of next tau up on grid
				// Repeat until locating which tau corresponds to Q(tau) ~ Y_i
				while(resLin[i] > QPos && l < L-mid-1) {
					QPosold = QPos;
					l++;
					for(QPos = Q0Pos[l], j = 0; j < p; j++) QPos += x[i][j] * bPos[j][l];
				}
				// if point located above largest grid tau...
				if(l == L - mid - 1) {
					rp[i] = F0tail(Q0tail(taugrid[L-2], nu) + (resLin[i] - QPosold)/sigmat1, nu);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = lf0tail(Q0tail(taugrid[L-2], nu) + (resLin[i] - QPosold)/sigmat1, nu) - log(sigmat1); break;
					}
				}  else {
					// if located anywhere in estimated grid...
					for(qdot_lo = b0dot[mid+l-1], j = 0; j < p; j++) qdot_lo += bdot[j][mid+l-1] * x[i][j];
					for(qdot_up = b0dot[mid+l], j = 0; j < p; j++) qdot_up += bdot[j][mid+l] * x[i][j];
					rp[i] = find_tau_lo(resLin[i], QPosold, qdot_lo, qdot_up, taugrid[mid+l-1], taugrid[mid+l]);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = -log(part_trape_rp(rp[i], qdot_lo, qdot_up, taugrid[mid+l-1], taugrid[mid+l])); break;
					}
				}
			} else {
				l = 0;
				QNegold = 0.0;
				// Follow similar algorithm but working way from median to lower quantiles
				for(QNeg = Q0Neg[l], j = 0; j < p; j++) QNeg += x[i][j] * bNeg[j][l];
				while(resLin[i] < -QNeg && l < mid) {
					QNegold = QNeg;
					l++;
					for(QNeg = Q0Neg[l], j = 0; j < p; j++) QNeg += x[i][j] * bNeg[j][l];
				}
				if(l == mid) { // tail
					rp[i] = F0tail(Q0tail(taugrid[1], nu) + (resLin[i] + QNegold)/sigmat2, nu);
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = lf0tail(Q0tail(taugrid[1], nu) + (resLin[i] + QNegold)/sigmat2, nu) - log(sigmat2); break;
					}
				} else {
					for(qdot_lo = b0dot[mid-l], j = 0; j < p; j++) qdot_lo += bdot[j][mid-l] * x[i][j];
					for(qdot_up = b0dot[mid-l+1], j = 0; j < p; j++) qdot_up += bdot[j][mid-l+1] * x[i][j];
					rp[i] = find_tau_up(resLin[i], -QNegold, qdot_lo, qdot_up, taugrid[mid-l], taugrid[mid-l+1]); //2
					switch (cens[i]) {
					case 1: ll[i] = log(1.0 - rp[i]); break;
					case 2: ll[i] = log(rp[i]); break;
					default: ll[i] = -log(part_trape_rp(rp[i], qdot_lo, qdot_up, taugrid[mid-l], taugrid[mid-l+1])); break;
					}
				}
			}
			if(ll[i] == qt(1.0, 1.0, 1, 0)) Rprintf("i = %d, ll[i] = %g, resLin[i] = %g, l = %d\n", i, ll[i], resLin[i], l);
		}
	} else {
		for(i = 0; i < n; i++) ll[i] = 0.0;
	}
//Rprintvec("ll = ", "%g ", ll, n);
	lp[0] = temp * inprod(ll, wt, n); //Weighting each log-likelihood contribution

// Add contribution of log prior on nu/6 (a standard logistic distribution)
	if(!llonly) {
		//lp[0] += lps0 + dt(par[m*(p+1)], 1.0, 1) + dlogis(par[(m+1)*(p+1)], 0.0, 1.0, 1) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		//lp[0] += lps0 + dlogis(par[(m+1)*(p+1)], 0.0, 1.0, 1) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		lp[0] += lps0 + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		// lp[0] += lps0 - 6.0*par[(m+1)*(p+1)] - 4.0*exp(-par[(m+1)*(p+1)]) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		// lp[0] += lps0 - par[(m+1)*(p+1)] - exp(-par[(m+1)*(p+1)]) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		// lp[0] += lps0 + par[(m+1)*(p+1)] - exp(par[(m+1)*(p+1)]) + dlogis(par[(m+1)*(p+1)+1], 0.0, 1.0, 1);
		if(shrink) lp[0] += dbasejoint(par + m*(p+1) + 1);
		//else lp[0] -= 0.5*sumsquares(par + m*(p+1) + 1, p)/9.0;
	}
}
// End of evaluation of log posterior function


// Define global pointers



//---------Function for performing Bayesian Joint Quantile Regression--------//

// Main function for evaluation of Bayesian Joint Quantile Regression
// It sets up pointers and allocates memory to the places they point, it
// unpackages concatenated vectors and matrices from R, and calls the
// MCMC routine
// Arguments/Inputs/Outputs in C and their equivalences in R:
//  1 *par          R par
//  2 *xVar         R x
//  3 *yVar         R y
//  4 *status       R cens
//  5 *weights      R wt
//  6 *toShrink     R shrink
//  7 *hyper        R hyper      (asig, bsig, akap[1], bkap[1], lpkap[1], akap[2], bkap[2], lpkap[3],...)
//  8 *dim          R dimpars    (n, p, L, mid index of median, m, ngrid, nkap, niter, thin)
//  9 *gridpars     R gridmats   (one long vech of matrix (each column is vech A, vech R followed by ldRgrid and lpgrid corresponding to one lambda))
// 10 *tauG         R tau.g
// 11 *muVar        R muV=blocks.mu
// 12 *SVar         R SV=blocks.S
// 13 *blocksVar    R blocks=blocks.ix
// 13 *blocks_size  R blocks.size
// 14 *dmcmcpar     R dmcmc.par  ()
// 15 *imcmcpar     R imcmcpar=imcmc.par
// 16 *parsamp      R parsamp=nsamp*length(par)
// 17 *acptsamp     R acptsamp=hsamp*nblocks
// 18 *lpsamp       R lpsamp=length(nsamp)
// 19 *distribution R fbase.choice=as.integer(fbase.choice)

void BJQRgp(double *parmg, int *phiVar, double *rhoVar, double *phigrid, double *xVar, double *yVar, int *status, double *weights, int *toShrink, double *hyper, int *dim, double *gridpars,
            double *tauG, double *muVar, double *SVar, double *thetaVar, double *etaVar, int *blocksVar, int *blocks_size, double *dmcmcpar, int *imcmcpar,
            double *parmgsamp, double *acptsamp, double *lpsamp, double *lprsamp, double *lkappasamp, int *phisamp, double *acpthpsamp, double *lmglhsamp, double *usamp, int *distribution){

	// local constants setup
	// n: number of observations
	// p: number of predictors is design matrix, not inluding intercept
	// L: number of tau grid points
	// mid: index of tau in grid acting as median.  (in C indexed as one fewer than R)
	// m: number of knots to act in the interpolation approximation of each w_j
	// ngrid: number of points on lambda grid
	// nkap: number of columns in kappa

	int i, j, k, l;

	int reach = 0; //reach reused to move pointers along via reach++
	shrink = toShrink[0];
	n = dim[reach++]; p = dim[reach++]; L = dim[reach++]; mid = dim[reach++];
	m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++]; nphi = dim[reach++];/***/
	dist = distribution[0];
	int niter = dim[reach++], thin = dim[reach++], nparmg = (m+1)*(p+1) + 2;
	taugrid = tauG;

	asig = hyper[0]; bsig = hyper[1];
	akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
	for(reach = 2, i = 0; i < nkap; i++) {
		akap[i] = hyper[reach++];
		bkap[i] = hyper[reach++];
		lpkap[i] = hyper[reach++];
	}
	shrinkFactor = shrinkFn((double)p);

	arho = 1; brho = 1;
	// Create pointers and allocate memory to the places they point
	// Organize constants Agrid and Rgrid from vector to consummable array format
	reach = 0;
	Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
	Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
	ldRgrid = vect(ngrid);
	lpgrid = vect(ngrid);

	for(i = 0; i < ngrid; i++) {
		Agrid[i] = mymatrix(L, m);
		for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++]; // put Agrid from vector form into array format

		Rgrid[i] = mymatrix(m, m);
		for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++]; // put Rgrid from vector form into array format

		ldRgrid[i] = gridpars[reach++];
		lpgrid[i] = gridpars[reach++];
	}

	reach = 0;
	Lphigrid = mymatrix(nphi, n);
	Gphigrid = array3d(nphi, n, n);

	for(i = 0; i < nphi; i++) {
		for(l = 0; l < n; l++) Lphigrid[i][l] = phigrid[reach++];
		for(k = 0; k < n; k++) for(l = 0; l < n; l++) Gphigrid[i][l][k] = phigrid[reach++];
	}


	reach = 0;
	x = mymatrix(n, p);
	for(j = 0; j < p; j++) for(i = 0; i < n; i++) x[i][j] = xVar[reach++];
	y = yVar;
	cens = status;
	wt = weights;
	int phi = phiVar[0];
	double rho = rhoVar[0];

	// Create pointers and space for objects that will be evaluated in lpFn
	lb = vect(10);
	wgrid = mymatrix(ngrid, L);
	lw = vect(nkap);
	llgrid = vect(ngrid);
	zknot = vect(m);
	wMat = mymatrix(p, L);
	vMat = mymatrix(p, L);
	vTilde = mymatrix(p, L);
	w0 = vect(L);
	zeta0dot = vect(L);
	zeta0 = vect(L);
	vNormSq = vect(L);
	a = mymatrix(L, n);
	aX = vect(L);
	gam = vect(p);
	xLin = vect(n);
	resLin = vect(n);
	b0dot = vect(L);
	bdot = mymatrix(p, L);
	Q0Pos = vect(L);
	bPos = mymatrix(p, L);
	Q0Neg = vect(L);
	bNeg = mymatrix(p, L);
	llvec = vect(n);
	rpvec = vect(n);
	pgvec = vect(ngrid*(p+1));
	zeta0_tick = ivect(L);
	zeta0_dist = vect(L);
	lp = vect(1);
	zcop = vect(n);
	zprod = vect(n);

	// Unpack parameters related to running of MCMC
	int b;
	int nblocks = imcmcpar[0], refresh = imcmcpar[1], verbose = imcmcpar[2], ticker = imcmcpar[3], refreshhp = imcmcpar[4], refreshhp_counter = imcmcpar[5];
	int *refresh_counter = imcmcpar + 6;

	double temp = dmcmcpar[0], decay = dmcmcpar[1], acpt_hptarget = dmcmcpar[2], lmhp = dmcmcpar[3];
	double *acpt_target = dmcmcpar + 4, *lm = acpt_target + nblocks;

	//Initialize mu[b] and S[b] which get used in adaptive MCMC for each block's update
	double **mu = (double **)R_alloc(nblocks, sizeof(double *));
	double ***S = (double ***)R_alloc(nblocks, sizeof(double **));
	int **blocks = (int **)R_alloc(nblocks, sizeof(int *));

	int mu_point = 0, S_point = 0, blocks_point = 0;
	for(b = 0; b < nblocks; b++) {
		mu[b] = muVar + mu_point; mu_point += blocks_size[b]; // Initialize mu with values passed from R
		S[b] = (double **)R_alloc(blocks_size[b], sizeof(double *)); // Initialize  S with values passed from R
		for(i = 0; i < blocks_size[b]; i++) {
			S[b][i] = SVar + S_point;
			S_point += blocks_size[b];
		}
		blocks[b] = blocksVar + blocks_point; blocks_point += blocks_size[b]; // blocks[b] holds indices of params needing to be updated in each block
	}

	double theta = thetaVar[0], eta = etaVar[0];

	//Call to adaptive MCMC code
	adMCMCgp(niter, thin, n, nparmg, parmg, phi, rho, nblocks, blocks, blocks_size, mu, S, theta, eta, acpt_target, acpt_hptarget, decay, refresh, refreshhp, lm, lmhp, temp,
	         refresh_counter, refreshhp_counter, verbose, ticker, parmgsamp, acptsamp, lpsamp, lprsamp, lkappasamp, phisamp, acpthpsamp, lmglhsamp, usamp);
}


//-------------Search Golden Selection - Maximization Routine-----------------//

// Define squareroot 5 which gets used in Max_Search_Golden_Selection
#define sqrt5 2.236067977499789696

// Function to stop if within tolerence
// Returns 1 if absolute difference is less than tolerance and 0 if it is greater
static int Stopping_Rule(double x0, double x1, double tolerance){
	double xm = 0.5 * fabs( x1 + x0 ); //fabs is absolute value

	if ( xm <= 1.0 ) return ( fabs( x1 - x0 ) < tolerance ) ? 1 : 0;
	return ( fabs( x1 - x0 ) < tolerance * xm ) ? 1 : 0;
}

// Function to find max via the Search_Golden_Section algorithm
void Max_Search_Golden_Section( double (*f)(double, double*, int*), double* a, double *fa, double* b, double* fb, double *y, int *cens, double tolerance){
	static const double lambda = 0.5 * (sqrt5 - 1.0);
	static const double mu = 0.5 * (3.0 - sqrt5); // = 1 - lambda
	double x1;
	double x2;
	double fx1;
	double fx2;


	x1 = *b - lambda * (*b - *a);
	x2 = *a + lambda * (*b - *a);
	fx1 = f(x1, y, cens); //pointer to input function
	fx2 = f(x2, y, cens); //pointer to input function

	// If tolerance specificed as negative, set to default, minimal tolerance
	if (tolerance <= 0.0) tolerance = 1.0e-5 * (*b - *a);

	while ( !Stopping_Rule( *a, *b, tolerance) ) {
		if (fx1 < fx2) {
			*a = x1;
			*fa = fx1;
			if ( Stopping_Rule( *a, *b, tolerance) ) break; // if within absolute tolerence, stop
			x1 = x2; // otherwise move closer & continue
			fx1 = fx2;
			x2 = *b - mu * (*b - *a);
			fx2 = f(x2, y, cens);
		} else {
			*b = x2;
			*fb = fx2;
			if ( Stopping_Rule( *a, *b, tolerance) ) break;
			x2 = x1;
			fx2 = fx1;
			x1 = *a + mu * (*b - *a);
			fx1 = f(x1, y, cens);
		}
	}
	return;
}



//-------------------------------------------------------------------------//
//---------------Functions for initializing MCMC paramters-----------------//

int sig_pos;
double *par0;

// Shortened evaluation of logpostFn; defaults to temp==1, exclusion of prior on nu
// Function to be maximized over sigma when finding starting sigma values
double lpFn2_noX(double sigma, double *y, int *cens){
	par0[sig_pos] = sigma;
	logpostFn_noX(par0, 1.0, 0, llvec, pgvec, y, cens, lp, rpvec);
	return(lp[0]);
}

double lpFn2(double sigma, double *y, int *cens){
	par0[sig_pos] = sigma;
	logpostFn(par0, 1.0, 0, llvec, pgvec, y, cens, lp, rpvec);
	return(lp[0]);
}


// Initialization function
// Performs many of same unpacking functionalities that BJQR does but additionally
// identifies starting value for sigma by maximizing the log posterior value
// over range of possible sigma values
// Arguments/Inputs/Outputs (se BJQR)
// 11 *siglim      R siglim  Upper and lower limits of sigma to search for optimal starting value

void INIT(double *par, double *xVar, double *yVar, int *status, double *weights, int *toShrink, double *hyper, int *dim, double *gridpars, double *tauG, double *siglim, int *distribution){

	int i, j, k, l;

	int reach = 0;
	shrink = toShrink[0];
	n = dim[reach++]; p = dim[reach++]; L = dim[reach++]; mid = dim[reach++];
	m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++];
	dist = distribution[0];

	taugrid = tauG;
	asig = hyper[0]; bsig = hyper[1];
	akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
	for(reach = 2, i = 0; i < nkap; i++) {
		akap[i] = hyper[reach++];
		bkap[i] = hyper[reach++];
		lpkap[i] = hyper[reach++];
	}
	shrinkFactor = shrinkFn((double)p);

	reach = 0;
	Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
	Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
	ldRgrid = vect(ngrid);
	lpgrid = vect(ngrid);

	for(i = 0; i < ngrid; i++) {
		Agrid[i] = mymatrix(L, m);
		for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++];

		Rgrid[i] = mymatrix(m, m);
		for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++];

		ldRgrid[i] = gridpars[reach++];
		lpgrid[i] = gridpars[reach++];
	}

	reach = 0;
	x = mymatrix(n, p);
	for(j = 0; j < p; j++) for(i = 0; i < n; i++) x[i][j] = xVar[reach++];
	cens = status;
	wt = weights;

	lb = vect(10);
	wgrid = mymatrix(ngrid, L);
	lw = vect(nkap);
	llgrid = vect(ngrid);
	zknot = vect(m);
	wMat = mymatrix(p, L);
	vMat = mymatrix(p, L);
	vTilde = mymatrix(p, L);
	w0 = vect(L);
	zeta0dot = vect(L);
	zeta0 = vect(L);
	vNormSq = vect(L);
	a = mymatrix(L, n);
	aX = vect(L);
	gam = vect(p);
	xLin = vect(n);
	resLin = vect(n);
	b0dot = vect(L);
	bdot = mymatrix(p, L);
	Q0Pos = vect(L);
	bPos = mymatrix(p, L);
	Q0Neg = vect(L);
	bNeg = mymatrix(p, L);
	llvec = vect(n);
	rpvec = vect(n);
	pgvec = vect(ngrid*(p+1));
	zeta0_tick = ivect(L);
	zeta0_dist = vect(L);
	lp = vect(1);

	sig_pos = (m+1)*(p+1); //position in par of sigma
	par0 = par;

	double sig_a = siglim[0], sig_b = siglim[1];
	double fa = lpFn2(sig_a, yVar, cens);
	double fb = lpFn2(sig_b, yVar, cens);

	// Use Max Golden Section Search algorithm to maximize log likelihood over sigma
	Max_Search_Golden_Section(lpFn2, &sig_a, &fa, &sig_b, &fb, yVar, cens, 1.0e-5);
	par[sig_pos] = (sig_a + sig_b) / 2.0; // midpoint of two points between which max lpFn2 occurs
}


// Function to calculate post-hoc the Deviance at each iteration;
// called within the summary.qrjoint function in R
void DEV(double *par, double *xVar, double *yVar, int *status, double *weights, int *toShrink, double *hyper, int *dim, double *gridpars, double *tauG, double *devsamp, double *llsamp, double *pgsamp, double *rpsamp, int *distribution){

	int i, j, k, l;

	int reach = 0;
	shrink = toShrink[0];
	n = dim[reach++]; p = dim[reach++]; L = dim[reach++]; mid = dim[reach++];
	m = dim[reach++]; ngrid = dim[reach++]; nkap = dim[reach++]; nphi = dim[reach++];
	dist = distribution[0];
	int niter = dim[reach++], nparmg = (m+1)*(p+1) + 2;

	taugrid = tauG;
	asig = hyper[0]; bsig = hyper[1];
	akap = vect(nkap); bkap = vect(nkap); lpkap = vect(nkap);
	for(reach = 2, i = 0; i < nkap; i++) {
		akap[i] = hyper[reach++];
		bkap[i] = hyper[reach++];
		lpkap[i] = hyper[reach++];
	}
	shrinkFactor = shrinkFn((double)p);

	reach = 0;
	Agrid = (double ***)R_alloc(ngrid, sizeof(double **));
	Rgrid = (double ***)R_alloc(ngrid, sizeof(double **));
	ldRgrid = vect(ngrid);
	lpgrid = vect(ngrid);

	for(i = 0; i < ngrid; i++) {
		Agrid[i] = mymatrix(L, m);
		for(l = 0; l < L; l++) for(k = 0; k < m; k++) Agrid[i][l][k] = gridpars[reach++];

		Rgrid[i] = mymatrix(m, m);
		for(k = 0; k < m; k++) for(l = 0; l < m; l++) Rgrid[i][l][k] = gridpars[reach++];

		ldRgrid[i] = gridpars[reach++];
		lpgrid[i] = gridpars[reach++];
	}

	reach = 0;
	x = mymatrix(n, p);
	for(j = 0; j < p; j++) for(i = 0; i < n; i++) x[i][j] = xVar[reach++];
	cens = status;
	wt = weights;

	lb = vect(10);
	wgrid = mymatrix(ngrid, L);
	lw = vect(nkap);
	llgrid = vect(ngrid);
	zknot = vect(m);
	wMat = mymatrix(p, L);
	vMat = mymatrix(p, L);
	vTilde = mymatrix(p, L);
	w0 = vect(L);
	zeta0dot = vect(L);
	zeta0 = vect(L);
	vNormSq = vect(L);
	a = mymatrix(L, n);
	aX = vect(L);
	gam = vect(p);
	xLin = vect(n);
	resLin = vect(n);
	b0dot = vect(L);
	bdot = mymatrix(p, L);
	Q0Pos = vect(L);
	bPos = mymatrix(p, L);
	Q0Neg = vect(L);
	bNeg = mymatrix(p, L);
	zeta0_tick = ivect(L);
	zeta0_dist = vect(L);
	lp = vect(1);

	reach = 0;
	int iter, reach2 = 0, reach3 = 0;
	for(iter = 0; iter < niter; iter++) {
		logpostFn(par + reach, 1.0, 1, llsamp + reach2, pgsamp + reach3, yVar, cens, lp, rpsamp + reach2);
		devsamp[iter] = -2.0 * lp[0];
		reach += nparmg; reach2 += n; reach3 += ngrid * (p+1);
	}
}


//----------------------------------------------------------------------------//
//--------Functions for commonly used routines "mydefs" by Surya Tokdar-------//
//----------------------------------------------------------------------------//

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



// Note: not used in current code
// double matern(double x, double phi, double kappa){
//         /*Returns the Matern function for x, phi and kappa.*/
//
//         /* Variables */
//         double ans, cte;
//         double uphi=x/phi;
//
//         /* Matern */
//
//         if (uphi==0.0) return 1;
//         else{
//                 if (kappa==0.5)
//                         ans = exp(-uphi);
//                 else {
//                         cte = R_pow(2, (-(kappa-1)))/gammafn(kappa);
//                         ans = cte * R_pow(uphi, kappa) * bessel_k(uphi,kappa,1);
//                 }
//         }
//
//         return ans;
// }





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

//----------------------------------------------------------------------------//
//-----------------------Adaptive MCMC Using MVN proposals--------------------//
//-------------------------and overlapping block updates----------------------//
//----------------------------------------------------------------------------//


void adMCMCgp(int niter, int thin, int n, int nparmg, double *parmg, int phi, double rho, int nblocks, int **blocks, int *blocks_size, double **mu, double ***S, double theta, double eta, double *acpt_target, double acpt_hptarget, double decay, int refresh, int refreshhp, double *lm, double lmhp, double temp,
              int *refresh_counter, int refreshhp_counter, int verbose, int ticker,
              double *parmgsamp, double *acptsamp, double *lpsamp, double *lprsamp, double *lkappasamp, int *phisamp, double *acpthpsamp, double *lmglhsamp, double *usamp){

	// Get cholesky decomps of each block's S (proposal Gaussian covariance), store in R[b]
// Reserve memory for blocks_rank, plock_pivot, and blocks_d.  Space needed for decomps
	int b, i, j;
	double ***R = (double ***)R_alloc(nblocks, sizeof(double **));
	int **blocks_pivot = (int **)R_alloc(nblocks, sizeof(double *)), *blocks_rank = ivect(nblocks);
	double *blocks_d = vect(nparmg);
	for(b = 0; b < nblocks; b++) {
		R[b] = mymatrix(blocks_size[b], blocks_size[b]);
		blocks_pivot[b] = ivect(blocks_size[b]);
		chol(R[b], blocks_size[b], 1.0e-8, blocks_pivot[b], blocks_rank + b, blocks_size[b], blocks_d, S[b], 0, 0, 1.0e-10);
	}


// Reserve memory (via pointers) for
	int *chunk_size = ivect(nblocks); // holds period of adaptive update for each block
	double *acpt_chunk = vect(nblocks); // acceptance rate
	double **parbar_chunk = (double **)R_alloc(nblocks, sizeof(double *));
	double *alpha = vect(nblocks);
	double *frac = vect(nblocks);
// copula part
	int hp_size = 0;
	double acpt_hp = 0.0, hpbar = 0.0, alphahp, frachp = sqrt(1.0 / ((double)refreshhp_counter + 1.0));

	for(b = 0; b < nblocks; b++) {
		chunk_size[b] = 0;
		acpt_chunk[b] = 0.0;
		parbar_chunk[b] = vect(blocks_size[b]);
		for(i = 0; i < blocks_size[b]; i++) parbar_chunk[b][i] = 0.0;
		frac[b] = sqrt(1.0 / ((double)refresh_counter[b] + 1.0));
	}

	double lp_diff, lmglh, lmglhnew, lcop, lcopnew, lpr, lprnew, lkappanew, kappasq, rhonew;
	double **parmgstore = mymatrix(2, nparmg), *par_incr = vect(nparmg), *zsamp = vect(nparmg);
	double *lprgrid = vect(nphi), *lcopgrid = vect(nphi), *lphigrid = vect(nphi), *parmgtemp = vect(nparmg);
	double **U = mymatrix(2, n);
	int ipar = 0, iparnew = 1;
	double lkappa = log(rho) + parmg[nparmg - 2];

	for(i = 0; i < nparmg; i++) {
		parmgstore[ipar][i] = parmg[i];
		parmgtemp[i] = parmg[i];
	}
	parmgstore[ipar][nparmg - 2] = log(1.0 - rho) + parmg[nparmg - 2];

	logpostFn(parmgtemp, temp, 0, llvec, pgvec, y, cens, lp, rpvec);
	lmglh = lp[0];
	for(i = 0; i < n; i++) U[ipar][i] = rpvec[i];

	lcop = lgauCopDen(n, U[ipar], phi, rho);
	lpr = lprior(n, phi, lkappa, parmgstore[ipar][nparmg - 2]);


	if(verbose) Rprintf("Initial lp = %g\n", lmglh + lcop + lpr);
// Rprintf("Initial lmglh = %g\n", lmglh);
// Rprintf("Initial lcop = %g\n", lcop);
// Rprintf("Initial lpr = %g\n", lpr);
// Rprintf("Initial phi = %d\n", phi);
// Rprintf("Initial lrho = %g\n", lrho);

	int iter, store_lp = 0, store_parmg = 0, store_acpt = 0, store_lmglh = 0, store_lpr = 0, store_lkappa = 0, store_phi = 0, store_acpthp = 0, store_u = 0;
	double *alpha_run = vect(nblocks), alphahp_run = 0.0, chs, lambda;
	int *run_counter = ivect(nblocks), runhp_counter = 0;
	for(b = 0; b < nblocks; b++) {alpha_run[b] = 0.0; run_counter[b] = 0;}

// joint part parameters
	int joint_size = 0, runjoint_counter = 0, refreshjoint = 3, refreshjoint_counter = 0, *pivotjoint = ivect(2), *rankjoint = ivect(1);
	double **Sj = mymatrix(2,2), **Rj = mymatrix(2,2), *joint_incr = vect(2), lmjoint = 2.38/sqrt(2.0), alphajoint, alphajoint_run = 0.0, acpt_joint = 0.0, *jointbar = vect(2), fracjoint = sqrt(1.0 / ((double)refreshjoint_counter + 1.0)), *muj = vect(2), *joint_d= vect(2), acpt_jtarget = 0.15;

	Sj[0][0] = 1.0; Sj[1][1] = 1.0; Sj[0][1] = 0.0, Sj[1][0] = 0.0;
	Rj[0][0] = 1.0; Rj[1][1] = 1.0; Rj[0][1] = 0.0, Rj[1][0] = 0.0;
	jointbar[0] = 0.0; jointbar[1] = 0.0;
	muj[0] = 0.0, muj[1] = 0.0;

	clock_t begin = clock();
	GetRNGstate();
	for(iter = 0; iter < niter; iter++) {
		for(b = 0; b < nblocks; b++) {
			chunk_size[b]++;
			for(i = 0; i < blocks_size[b]; i++) zsamp[i] = rnorm(0.0, 1.0);
			triprod(R[b], blocks_size[b], blocks_size[b], zsamp, par_incr, 1);
			for(i = 0; i < nparmg; i++) parmgstore[iparnew][i] = parmgstore[ipar][i];
			lambda = lm[b] * sqrt(3.0 / rgamma(3.0, 1.0));
			for(i = 0; i < blocks_size[b]; i++) parmgstore[iparnew][blocks[b][i]] += lambda * par_incr[i];

			for(i = 0; i < nparmg; i++) parmgtemp[i] = parmgstore[iparnew][i];
			kappasq = exp(lkappa);
			rhonew = kappasq/(kappasq + exp(parmgstore[iparnew][nparmg - 2]));
			parmgtemp[nparmg - 2] = parmgstore[iparnew][nparmg - 2] - log(1.0 - rhonew);

			logpostFn(parmgtemp, temp, 0, llvec, pgvec, y, cens, lp, rpvec);
			lmglhnew = lp[0];

			for(i = 0; i < n; i++) U[iparnew][i] = rpvec[i];
			lcopnew = lgauCopDen(n, U[iparnew], phi, rhonew);
			lprnew = lprior(n, phi, lkappa, parmgstore[iparnew][nparmg - 2]);

			lp_diff = lmglhnew - lmglh + lcopnew - lcop + lprnew - lpr;
			alpha[b] = exp(lp_diff); if(alpha[b] > 1.0) alpha[b] = 1.0;
			if(runif(0.0, 1.0) < alpha[b]) {
				ipar = iparnew;
				iparnew = !ipar;
				lmglh = lmglhnew;
				lcop = lcopnew;
				lpr = lprnew;
				rho = rhonew;
			}
			alpha_run[b] = ((double)run_counter[b] * alpha_run[b] + alpha[b]) / ((double)(run_counter[b] + 1.0));
			run_counter[b]++;
		}


		hp_size++;
		zsamp[0] = rnorm(0.0, 1.0);
		lambda = lmhp * sqrt(3.0 / rgamma(3.0, 1.0));
		lkappanew = lkappa + lambda*zsamp[0]*eta;
		kappasq = exp(lkappanew);
		rhonew = kappasq/(kappasq + exp(parmgstore[ipar][nparmg - 2]));

		for(i = 0; i < nparmg; i++) parmgtemp[i] = parmgstore[ipar][i];
		parmgtemp[nparmg - 2] = parmgstore[ipar][nparmg - 2] - log(1.0 - rhonew);
		logpostFn(parmgtemp, temp, 0, llvec, pgvec, y, cens, lp, rpvec);
		lmglhnew = lp[0];

		lcopnew = lgauCopDen(n, rpvec, phi, rhonew);
		lprnew = lprior(n, phi, lkappanew, parmgstore[ipar][nparmg - 2]);

		// Rprintf("lgtr = %g", lgtrstore[iparrnew[g]][g]);
		// Rprintf("lcopnew = %g", lcopnew[g]);
		// Rprintf("lprnew = %g", lprnew[g]);
		// Rprintf("\n");
		lp_diff = lmglhnew - lmglh + lcopnew - lcop + lprnew - lpr;
		alphahp = exp(lp_diff); if(alphahp > 1.0) alphahp = 1.0;
		if(runif(0.0, 1.0) < alphahp) {
			lcop = lcopnew;
			lpr = lprnew;
			lmglh = lmglhnew;
			lkappa = lkappanew;
			rho = rhonew;
			for(i = 0; i < n; i++) U[ipar][i] = rpvec[i];
		}
		alphahp_run = ((double)runhp_counter * alphahp_run + alphahp) / ((double)(runhp_counter + 1.0));
		runhp_counter++;


		for(i = 0; i < nphi; i++) {
			lprgrid[i] = lprior(n, i, lkappa, parmgstore[ipar][nparmg - 2]);
			lcopgrid[i] = lgauCopDen(n, U[ipar], i, rho);
			lphigrid[i] = lprgrid[i] + lcopgrid[i];
		}
		// Rprintf("lphi1 = %g, lphi2 = %g ", lphigrid[0], lphigrid[1]);
		phi = rdraw(nphi, lphigrid, 1);
		lpr = lprgrid[phi];
		lcop = lcopgrid[phi];

		joint_size++;
		for(i = 0; i < 2; i++) zsamp[i] = rnorm(0.0, 1.0);
		triprod(Rj, 2, 2, zsamp, joint_incr, 1);
		for(i = 0; i < nparmg; i++) parmgstore[iparnew][i] = parmgstore[ipar][i];


		parmgstore[iparnew][nparmg - 2] += lmjoint * joint_incr[0];
		lkappanew = lkappa + lmjoint * joint_incr[1];

		for(i = 0; i < nparmg; i++) parmgtemp[i] = parmgstore[iparnew][i];
		kappasq = exp(lkappanew);
		rhonew = kappasq/(kappasq + exp(parmgstore[iparnew][nparmg - 2]));
		parmgtemp[nparmg - 2] = parmgstore[iparnew][nparmg - 2] - log(1.0 - rhonew);

		logpostFn(parmgtemp, temp, 0, llvec, pgvec, y, cens, lp, rpvec);
		lmglhnew = lp[0];
		for(i = 0; i < n; i++) U[iparnew][i] = rpvec[i];
		lcopnew = lgauCopDen(n, U[iparnew], phi, rhonew);
		lprnew = lprior(n, phi, lkappanew, parmgstore[iparnew][nparmg - 2]);
		lp_diff = lmglhnew - lmglh + lcopnew - lcop + lprnew - lpr;
		alphajoint = exp(lp_diff); if(alphajoint > 1.0) alphajoint = 1.0;
		if(runif(0.0, 1.0) < alphajoint) {
			ipar = iparnew;
			iparnew = !ipar;
			lcop = lcopnew;
			lpr = lprnew;
			lkappa = lkappanew;
			lmglh = lmglhnew;
			rho = rhonew;
		}
		alphajoint_run = ((double)runjoint_counter * alphajoint_run + alphajoint) / ((double)(runjoint_counter + 1.0));
		runjoint_counter++;


		if((iter + 1) % thin == 0) {
			lpsamp[store_lp++] = lmglh + lcop + lpr;
			lmglhsamp[store_lmglh++] = lmglh;
			lprsamp[store_lpr++] = lpr;
			for(i = 0; i < nparmg; i++) {
				if (i != (nparmg - 2)) {
					parmgsamp[store_parmg++] = parmgstore[ipar][i];
				} else {
					parmgsamp[store_parmg++] = parmgstore[ipar][i] - log(1.0 - rho);
				}
			}
			for(b = 0; b < nblocks; b++) acptsamp[store_acpt++] = alpha[b];
			phisamp[store_phi++] = phi;
			lkappasamp[store_lkappa++] = lkappa;
			acpthpsamp[store_acpthp++] = alphahp;
			for(i = 0; i < n; i++) usamp[store_u++] = U[ipar][i];
		}

		if(verbose) {
			if(niter < ticker || (iter + 1) % (niter / ticker) == 0) {
				Rprintf("iter = %d, lp = %g \n", iter + 1, lmglh + lcop + lpr);
				// Rprintf("iter = %d, lcop = %g ", iter + 1, lcop);
				// Rprintf("iter = %d, lpr = %g ", iter + 1, lpr);

				//Rprintvec("acpt = ", "%0.2f ", alpha_run, nblocks);
				for(b = 0; b < nblocks; b++) {
					alpha_run[b] = 0.0;
					run_counter[b] = 0;
				}
				alphajoint_run = 0.0;
				runjoint_counter = 0;
				alphahp_run = 0.0;
				runhp_counter = 0;
			}
		}

		for(b = 0; b < nblocks; b++) {
			chs = (double) chunk_size[b]; if(chs < 1.0) chs = 1.0;
			acpt_chunk[b] += (alpha[b] - acpt_chunk[b]) / chs;
			for(i = 0; i < blocks_size[b]; i++)
				parbar_chunk[b][i] += (parmgstore[ipar][blocks[b][i]] - parbar_chunk[b][i]) / chs;
			if(chunk_size[b] == refresh * blocks_size[b]) {
				refresh_counter[b]++;
				frac[b] = sqrt(1.0 / ((double)refresh_counter[b] + 1.0));
				for(i = 0; i < blocks_size[b]; i++) {
					for(j = 0; j < i; j++) {
						S[b][i][j] = (1.0 - frac[b]) * S[b][i][j] + frac[b] * (parbar_chunk[b][i] - mu[b][i]) * (parbar_chunk[b][j] - mu[b][j]);
						S[b][j][i] = S[b][i][j];
					}
					S[b][i][i] = (1.0 - frac[b]) * S[b][i][i] + frac[b] * (parbar_chunk[b][i] - mu[b][i]) * (parbar_chunk[b][i] - mu[b][i]);
				}
				chol(R[b], blocks_size[b], 1.0e-8, blocks_pivot[b], blocks_rank + b, blocks_size[b], blocks_d, S[b], 0, 0, 1.0e-10);
				for(i = 0; i < blocks_size[b]; i++) {
					mu[b][i] += frac[b] * (parbar_chunk[b][i] - mu[b][i]);
					parbar_chunk[b][i] = 0.0;
				}
				lm[b] *= exp(frac[b] * (acpt_chunk[b] - acpt_target[b]));
				acpt_chunk[b] = 0.0;
				chunk_size[b] = 0;
			}
		}


		chs = (double) joint_size; if(chs < 1.0) chs = 1.0;
		acpt_joint += (alphajoint - acpt_joint) / chs;
		jointbar[0] += (parmgstore[ipar][nparmg - 2] - jointbar[0]) / chs;
		jointbar[1] += (lkappa - jointbar[1]) / chs;
		if(joint_size == refreshjoint * 2) {
			refreshjoint_counter++;
			fracjoint = sqrt(1.0 / ((double)refreshjoint_counter + 1.0));
			Sj[0][1] = (1.0 - fracjoint) * Sj[0][1] + fracjoint * (jointbar[0] - muj[0]) * (jointbar[1] - muj[1]);
			Sj[1][0] = Sj[0][1];
			Sj[0][0] = (1.0 - fracjoint) * Sj[0][0] + fracjoint * (jointbar[0] - muj[0]) * (jointbar[0] - muj[0]);
			Sj[1][1] = (1.0 - fracjoint) * Sj[1][1] + fracjoint * (jointbar[1] - muj[1]) * (jointbar[1] - muj[1]);

			chol(Rj, 2, 1.0e-8, pivotjoint, rankjoint, 2, joint_d, Sj, 0, 0, 1.0e-10);
			for(i = 0; i < 2; i++) {
				muj[i] += fracjoint * (jointbar[i] - muj[i]);
				jointbar[i] = 0.0;
			}
			lmjoint *= exp(fracjoint * (acpt_joint - acpt_jtarget));
			acpt_joint = 0.0;
			joint_size = 0;
		}

		chs = (double) hp_size; if(chs < 1.0) chs = 1.0;
		acpt_hp += (alphahp - acpt_hp) / chs;
		hpbar += (lkappa - hpbar) / chs;
		if(hp_size == refreshhp) {
			refreshhp_counter++;
			frachp = sqrt(1.0 / ((double)refreshhp_counter + 1.0));
			eta = sqrt((1.0 - frachp) * eta * eta + frachp * (hpbar - theta) * (hpbar - theta));
			theta += frachp * (hpbar - theta);
			hpbar = 0.0;
			lmhp *= exp(frachp * (acpt_hp - acpt_hptarget));
			acpt_hp = 0.0;
			hp_size = 0;
		}
	}
	PutRNGstate();
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	Rprintf("time = %g ", time_spent);
}

void transform_grid(double *w, double *v, int *ticks, double *dists){
	int l;
	for(l = 0; l < L; l++) v[l] = (1.0 - dists[l]) * w[ticks[l]] + dists[l] * w[ticks[l]+1];
}

// Locates value "target" (in relation to baseline). Uses quadratic polynomial to approximate Q.
// Employs constraint such that 'a' and 'b' are derivative of Q at lower and upper delta bounds.
// And baseline is Q(taua)
double find_tau_lo(double target, double baseline, double a, double b, double taua, double taub){
	double loc, Delta = taub - taua;
	if (fabs(b-a)>1.0e-15) {
		loc = ((b*taua - a*taub) + sqrt(a*a * Delta*Delta + 2*Delta*(b-a)*(target - baseline)))/(b-a);
	} else { // when slopes equal, linear interpolant
		loc = taua + (target - baseline)/a;
	}
	return(loc);
}

// Version for approximating Q in the QNeg region. Also employs constraints such that
// 'a' and 'b' are derivative of Q at lower and upper grid points. This one additionally
// constrains so that baseline is Q(taub).  Will need target = res, baseline= -QNeg
double find_tau_up(double target, double baseline, double a, double b, double taua, double taub){
	double loc, Delta = taub - taua;
	if (fabs(b-a)>1.0e-15) {
		loc = ((b*taua - a*taub) + sqrt(b*b * Delta*Delta + 2*Delta*(b-a)*(target - baseline)))/(b-a);
	} else { // when slopes equal, linear interpolant
		loc = taub - (target - baseline)/b;
	}
	return(loc);
}

// Given proportion and lower & upper derivatives at delta grid, approximates derivative of
// quantile function at y for likelihood contribution; use in conjunction with prop_tau
double part_trape_rp(double loc, double a, double b, double taua, double taub){
	return ((b*(loc - taua)/(taub - taua) + a*(taub - loc)/(taub - taua)));
}


double lprior(int p, int phi, double lkappa, double ltau){
	double result = 0.0, kapsq = exp(lkappa), tausq = exp(ltau);
	int i;
	double sigsq = kapsq + tausq, rho = kapsq/sigsq;

	for(i = 0; i < p; i++) result += log1p(rho*(Lphigrid[phi][i]-1.0));
	result = -result/2.0 + (arho - 1.0) * log(rho) + (brho - 1.0) * log(1.0 - rho) - 2.0 * log(sigsq) + log(kapsq) + log(tausq);
	return result;
}


double wtsumsquares(double *x, int phi, double rho, int n){
	double wss = 0.0;
	int i;

	for(i = 0; i < n; i++)
		wss += x[i] * x[i] * (1.0 / (rho*(Lphigrid[phi][i] - 1.0) + 1.0) - 1.0);
	return wss;
}

double lgauCopDen(int p, double *x, int phi, double rho){
	int i;

	for(i = 0; i < p; i++) zcop[i] = qnorm(x[i], 0.0, 1.0, 1, 0);
	// Rprintvec("z = ", "%0.2f ", zcop, n);

	mvprod(Gphigrid[phi], zcop, zprod, p, p, 1);

	return (-wtsumsquares(zprod, phi, rho, p)/2.0);
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
