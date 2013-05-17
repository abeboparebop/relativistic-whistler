#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_bessel.h>

/* These variables are stored as globals, so that the root-finder
   always has a good place to start the next root-finding. (This is only 
   if you aren't moving very far in parameter space.) */
double Or_guess;
double Oi_guess;
double kap_guess;

/* Datatypes to pass to various GSL functions. */
typedef struct {
  double complex O;
  double complex Upar_min;
  double complex Upar_max;
  double Uperp;
  double kap;
  double eta;
  double eperp;
  double bessel;
} Params_inner;

typedef struct {
  double complex O;
  double eps2;
  double kap;
  double eta;
  double eperp;
  double bessel;
} Params_outer;

typedef struct {
  double kap;
  double eta;
  double eperp;
  double bperp;
  double eps2;
  double bessel;
} Params_disp;

typedef struct {
  double eta;
  double eperp;
  double bperp;
  double eps2;
  double bessel;
} Params_gam;

typedef struct{
  double eta;
  double eperp;
  double bperp;
  double eps2;
  double Uperp;
  double bessel;
} Params_calc;

typedef struct {
  double gam0;
  double eperp;
  double bperp;
  double eps2;
  double bessel;
} Params_max;

/* Evaluation of the integrand of the integral over Upar and Uperp: */
double complex inner_calc(double complex Upar, Params_inner p) {
  double complex gam = csqrt(1+p.Uperp*p.Uperp+Upar*Upar);
  double complex gam_aniso = csqrt(1+p.Uperp*p.Uperp+Upar*Upar*p.eta);
  double frac = sqrt(p.eta/(2*p.eperp));
  double Uperp3 = p.Uperp*p.Uperp*p.Uperp;

  double complex exp_bessel;
  if (p.eperp < 0.05) {
    exp_bessel = cexp(1./p.eperp - gam_aniso/p.eperp) *
      (sqrt(2./M_PI/p.eperp) -			
       15./4.*sqrt(p.eperp/2./M_PI));
  }
  else
    exp_bessel = exp(-gam_aniso/p.eperp)/p.bessel;

  double complex result = (-exp_bessel * sqrt(p.eta) * Uperp3 *		
		   (Upar*(p.eta-1)*2*frac*p.kap + 2*gam*p.O) / 
		   (8.*gam*p.eperp*p.eperp*gam_aniso * 
		    (1 + Upar*frac*p.kap - gam*p.O) ));

  /* printf("%e\n", result); */

  return result;
}

/* Wrapper around the real part of integrand */
double real_integrand(double Upar, void * params) {
  Params_inner * p = (Params_inner *)params;
  double result = creal(inner_calc(Upar + 1j*0.0, *p));
  return result;
}

/* Wrapper around the real part of integrand, assuming that we're 
   off the real axis in Upar. (GSL integration can't handle this 
   case automatically.) */
double real_integrand_subs(double x, void * params) {
  Params_inner * p = (Params_inner *)params;
  double complex Upar = (p->Upar_max-p->Upar_min)*x + p->Upar_min;
  double complex gprime = p->Upar_max-p->Upar_min;
  double result = creal(gprime*inner_calc(Upar, *p));
  return result;
}

/* Wrapper around the imaginary part of integrand */
double imag_integrand(double Upar, void * params) {
  Params_inner * p = (Params_inner *)params;
  double result = cimag(inner_calc(Upar, *p));
  return result;
}

/* Wrapper around the imaginary part of integrand, assuming that we're 
   off the real axis in Upar. (GSL integration can't handle this 
   case automatically.) */
double imag_integrand_subs(double x, void * params) {
  Params_inner * p = (Params_inner *)params;
  double complex Upar = (p->Upar_max-p->Upar_min)*x + p->Upar_min;
  double complex gprime = p->Upar_max-p->Upar_min;
  double result = cimag(gprime*inner_calc(Upar, *p));
  return result;
}

/* At fixed Uperp, perform integration from negative to positive 
   infinity in Upar (real part) */
double real_integrand_out(double x, void * params) {
  /* handle incoming parameters */
  Params_outer * p = (Params_outer *)params;
 
  /* variables for handling inner integrations */
  int size = 10000;
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (size);
  double result, error;
  gsl_function F;
  double total = 0.0;

  /* variables for finding locations of poles */
  double complex term1, term2, denom;
  double complex pol1, pol2, pol_tmp;
  double rep1, rep2, imp1, imp2;

  /* How close do we integrate to the poles (in Upar)? */
  double eps = 1e-3;

  /* Find locations of the two poles in Upar */
  term1 = 2.0*p->eperp*sqrt(2.0*p->eta/p->eperp)*p->kap;
  term2 = csqrt(8*p->eperp*p->eta*p->kap*p->kap +
		8*p->eperp*(1-p->O*p->O - x*x*p->O*p->O)*
		(-p->eta*p->kap*p->kap + 2*p->eperp*p->O*p->O));
  denom = (-2*p->eta*p->kap*p->kap + 4*p->eperp*p->O*p->O);
  pol1 = (term1-term2)/denom;
  pol2 = (term1+term2)/denom;
  if (creal(pol1) > creal(pol2)) {
    pol_tmp = pol1;
    pol1 = pol2;
    pol2 = pol_tmp;
  }
  rep1 = creal(pol1);
  rep2 = creal(pol2);
  imp1 = cimag(pol1);
  imp2 = cimag(pol2);

  /* Accumulate sum of 9 integrals: */
  /* Full integral is sum of three segments along real axis, plus two notches */
  /* into the complex plane under the poles. quad doesn't natively handle */
  /* integration along contours in complex plane, so real_integrand_subs */
  /* uses a change of variables to turn each of the "edges" of the rectangular */
  /* notches into an integral along the real axis from 0 to 1. Each of */
  /* nND_args define the path taken. */

  Params_inner pass_params = {
    .Upar_min = 0.0, //will be modified upon subsequent integrations
    .Upar_max = rep1-eps,
    .Uperp = x,
    .O = p->O,
    .kap = p->kap,
    .eta = p->eta,
    .eperp = p->eperp,
    .bessel = p->bessel
  };

  double epsabs = p->eps2;
  double epsrel = 0.0;

  /* 1: negative infinity to left pole */
  F.function = &real_integrand;
  F.params = &pass_params;
  gsl_integration_qagil(&F, rep1-eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 2: left notch of left pole */
  pass_params.Upar_min = rep1-eps;
  pass_params.Upar_max = rep1-eps + I*(imp1-eps);
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
		       w, &result, &error); 
  total += result;

  /* 3: bottom notch of left pole */
  pass_params.Upar_min = rep1-eps + I*(imp1-eps);
  pass_params.Upar_max = rep1+eps + I*(imp1-eps);
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 4: right notch of left pole (return to real axis) */
  pass_params.Upar_min = rep1+eps + I*(imp1-eps);
  pass_params.Upar_max = rep1+eps;
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 5: from left to right pole */
  F.function = &real_integrand;
  F.params = &pass_params;
  gsl_integration_qags(&F, rep1+eps, rep2-eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 6: left notch of right pole */
  pass_params.Upar_min = rep2-eps;
  pass_params.Upar_max = rep2-eps + I*(imp2-eps);
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 7: bottom notch of right pole */
  pass_params.Upar_min = rep2-eps + I*(imp2-eps);
  pass_params.Upar_max = rep2+eps + I*(imp2-eps);
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 8: right notch of right pole (return to real axis) */
  pass_params.Upar_min = rep2+eps + I*(imp2-eps);
  pass_params.Upar_max = rep2+eps;
  F.function = &real_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 9: from right pole to positive infinity */
  F.function = &real_integrand;
  F.params = &pass_params;
  gsl_integration_qagiu(&F, rep2+eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  gsl_integration_workspace_free(w);
  
  return total;
}

/* At fixed Uperp, perform integration from negative to positive 
   infinity in Upar (imaginary part) */

double imag_integrand_out(double x, void * params) {
  /* handle incoming parameters */
  Params_outer * p = (Params_outer *)params;
 
  /* variables for handling inner integrations */
  int size = 10000;
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (size);
  double result, error;
  gsl_function F;
  double total = 0.0;

  /* How close do we integrate to the poles (in Upar)? */
  double eps = 1e-3;

  /* variables for finding locations of poles */
  double complex term1, term2, denom;
  double complex pol1, pol2, pol_tmp;
  double rep1, rep2, imp1, imp2;

  /* Find locations of the two poles */
  term1 = 2.0*p->eperp*sqrt(2.0*p->eta/p->eperp)*p->kap;
  term2 = csqrt(8*p->eperp*p->eta*p->kap*p->kap +
		8*p->eperp*(1-p->O*p->O - x*x*p->O*p->O)*
		(-p->eta*p->kap*p->kap + 2*p->eperp*p->O*p->O));
  denom = (-2*p->eta*p->kap*p->kap + 4*p->eperp*p->O*p->O);
  pol1 = (term1-term2)/denom;
  pol2 = (term1+term2)/denom;
  if (creal(pol1) > creal(pol2)) {
    pol_tmp = pol1;
    pol1 = pol2;
    pol2 = pol_tmp;
  }
  rep1 = creal(pol1);
  rep2 = creal(pol2);
  imp1 = cimag(pol1);
  imp2 = cimag(pol2);

  F.function = &imag_integrand_subs;

  /* Accumulate sum of 9 integrals: */
  /* Full integral is sum of three segments along real axis, plus two notches */
  /* into the complex plane under the poles. quad doesn't natively handle */
  /* integration along contours in complex plane, so real_integrand_subs */
  /* uses a change of variables to turn each of the "edges" of the rectangular */
  /* notches into an integral along the real axis from 0 to 1. Each of */
  /* nND_args define the path taken. */

  Params_inner pass_params = {
    .Upar_min = 0.0, //will be modified upon subsequent integrations
    .Upar_max = rep1-eps,
    .Uperp = x,
    .O = p->O,
    .kap = p->kap,
    .eta = p->eta,
    .eperp = p->eperp,
    .bessel = p->bessel
  };
  double epsabs = p->eps2;
  double epsrel = 0.0;

  /* 1: negative infinity to left pole */
  F.function = &imag_integrand;
  F.params = &pass_params;
  gsl_integration_qagil(&F, rep1-eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 2: left notch of left pole */
  pass_params.Upar_min = rep1-eps;
  pass_params.Upar_max = rep1-eps + I*(imp1-eps);
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
		       w, &result, &error); 
  total += result;

  /* 3: bottom notch of left pole */
  pass_params.Upar_min = rep1-eps + I*(imp1-eps);
  pass_params.Upar_max = rep1+eps + I*(imp1-eps);
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 4: right notch of left pole (return to real axis) */
  pass_params.Upar_min = rep1+eps + I*(imp1-eps);
  pass_params.Upar_max = rep1+eps;
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 5: from left to right pole */
  F.function = &imag_integrand;
  F.params = &pass_params;
  gsl_integration_qags(&F, rep1+eps, rep2-eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 6: left notch of right pole */
  pass_params.Upar_min = rep2-eps;
  pass_params.Upar_max = rep2-eps + I*(imp2-eps);
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 7: bottom notch of right pole */
  pass_params.Upar_min = rep2-eps + I*(imp2-eps);
  pass_params.Upar_max = rep2+eps + I*(imp2-eps);
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 8: right notch of right pole (return to real axis) */
  pass_params.Upar_min = rep2+eps + I*(imp2-eps);
  pass_params.Upar_max = rep2+eps;
  F.function = &imag_integrand_subs;
  F.params = &pass_params;
  gsl_integration_qags(&F, 0, 1, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  /* 9: from right pole to positive infinity */
  F.function = &imag_integrand;
  F.params = &pass_params;
  gsl_integration_qagiu(&F, rep2+eps, epsabs, epsrel, size,
			w, &result, &error); 
  total += result;

  gsl_integration_workspace_free(w);

  return total;
}

/* evaluate the dispersion relation (most of the computational
   effor is the integral over Uperp and internally Upar)*/

int dispersion(const gsl_vector * x, void * params, 
		  gsl_vector * f) {
  /* handle incoming parameters */
  Params_disp * p = (Params_disp *)params;
  const double Or = gsl_vector_get(x,0);
  const double Oi = gsl_vector_get(x,1);

  printf("k = %f, O_r = %e, O_i = %e\n", p->kap, Or, Oi);

  /* Variables to handle integration */
  int size = 10000;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (size);
  gsl_function F;

  /* Results variables */
  double result_r, error_r, result_i, error_i;
  double term1_r, term2_r, term2_i;


  Params_outer pass_params = {
    .O = Or + 1j*Oi,
    .kap = p->kap,
    .eta = p->eta,
    .eperp = p->eperp,
    .eps2 = p->eps2,
    .bessel = p->bessel
  };
  
  double epsabs = p->eps2;
  double epsrel = 0.0;

  /* Integrate real part over infinite domain */ 
  F.function = &real_integrand_out;
  F.params = &pass_params;

  gsl_integration_qagiu(&F, 0, epsabs, epsrel, size,
		       w, &result_r, &error_r); 

  /* Integrate imaginary part over infinite domain */ 
  F.function = &imag_integrand_out;

  gsl_integration_workspace_free(w);
  w = gsl_integration_workspace_alloc (size);

  gsl_integration_qagiu(&F, 0, epsabs, epsrel, size,
		       w, &result_i, &error_i); 

  term1_r = p->eta*p->kap*p->kap/p->bperp;
  term2_r = -2.0*p->eperp*(Or*Or-Oi*Oi)/p->bperp;
  term2_i = -2.0*p->eperp*2*Or*Oi/p->bperp;

  gsl_integration_workspace_free(w);

  /* Strictly, the latter terms here should be imaginary, but */
  /* we're just finding the root of this function, so everything will */
  /* need to be zero anyway: */
  gsl_vector_set(f, 0, result_r+term1_r+term2_r);
  gsl_vector_set(f, 1, result_i+term2_i);

  return GSL_SUCCESS;
}


double neg_gamma (const gsl_vector *v, void * params)
{
  Params_gam * p = (Params_gam *)params;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  double kap = gsl_vector_get(v,0);

  int status;
  size_t iter=0;

  Params_disp test_params = {
    .kap = fabs(kap),
    .eta = p->eta,
    .bperp = p->bperp,
    .eperp = p->eperp,
    .eps2 = p->eps2,
    .bessel = p->bessel
  };

  /* number of dimensions */
  const size_t n=2;

  /* We're finding a root of the dispersion function: */
  gsl_multiroot_function f = {&dispersion, n, &test_params};

  /* Initial guess */
  gsl_vector *O = gsl_vector_alloc(n);
  gsl_vector_set(O, 0, Or_guess);
  gsl_vector_set(O, 1, Oi_guess);

  /* Initializing root finder */
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, n);
  gsl_multiroot_fsolver_set(s, &f, O);

  /* Do actual root finding */
  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);

    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, p->eps2);
  }
  while (status == GSL_CONTINUE && iter < 1000);
  
  printf("status = %s\n", gsl_strerror(status));
  
  /* Set guess variables to have useful starting points for following 
     iterations */
  Or_guess = gsl_vector_get(s->x, 0);
  Oi_guess = gsl_vector_get(s->x, 1);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(O);

  /* Return (negative) growth rate, in order to minimize in gamma_max */
  return -Oi_guess;
}

/* Minimizes -gamma to find maximum growth rate over k for given set 
   of plasma properties */
int gamma_max (const gsl_vector *v, void * params, gsl_vector *f)
{
  /* Handle incoming parameters */
  Params_max * p = (Params_max *)params;
  int status;
  int iter = 0, max_iter = 100;

  /* Set minimizer */
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  gsl_vector *ss, *kap_temp;

  double Oi;
  double size;

  /* set initial guess for k */
  kap_temp = gsl_vector_alloc(1);
  gsl_vector_set(kap_temp, 0, kap_guess);
  ss = gsl_vector_alloc(1);
  gsl_vector_set_all(ss,0.001);

  Params_gam test_params = {
    .eta = gsl_vector_get(v, 0),
    .bperp = p->bperp,
    .eperp = p->eperp,
    .eps2 = p->eps2,
    .bessel = p->bessel
  };

  gsl_set_error_handler_off();

  /* Minimizing neg_gamma subject to test_params */
  minex_func.n = 1;
  minex_func.f = neg_gamma;
  minex_func.params = &test_params;

  /* Initialize minimizer */
  s = gsl_multimin_fminimizer_alloc(T,1);
  gsl_multimin_fminimizer_set(s, &minex_func, kap_temp, ss);

  printf("finding gamma_max for eta=%e\n",gsl_vector_get(v,0));

  /* Iteratively minimize, varying k */
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-6);

    if (status == GSL_SUCCESS)
      printf("Converged:\n");
    Oi = -neg_gamma(s->x, &test_params);
    printf("%d %e %e\n",iter, gsl_vector_get(s->x,0), Oi);

    /* Dispersion relation is insufficiently accurate to get
       better accuracy than ~1e-8, so finding such a growth rate
       is essentially impossible */
    if (fabs(Oi) < 1e-8) {
      status = GSL_SUCCESS;
    }
    
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  printf("found gamma_max=%e for eta=%e\n",Oi, gsl_vector_get(v,0));

  /* Store location of maximum gamma, for use as an initial guess
     in future iterations */
  kap_guess = gsl_vector_get(s->x,0);

  gsl_vector_free(kap_temp);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  /* Return growth rate minus target growth rate. find_eta will rootfind
     on this parameter to find parameters corresponding to target growth
     rate. */
  gsl_vector_set(f, 0, Oi - p->gam0);

  return status;
}

/* Given beta_perp (and other parameters), find eta=Tperp/Tpar that 
   corresponds to target growth rate */
double find_eta (double eta_guess, void * params) {
  /* Handle incoming parameters */
  Params_max * p = (Params_max *)params;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  double eta_fin;

  /* number of dimensions */
  const size_t n=1;

  /* Find zero in gamma_max, corresponding to when growth rate = target */
  gsl_multiroot_function f = {&gamma_max, n, p};
  gsl_vector *eta = gsl_vector_alloc(n);
  gsl_vector_set(eta, 0, eta_guess);

  int status, iter=0;

  /* Initialize rootfinder */
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, n);
  gsl_multiroot_fsolver_set(s, &f, eta);

  /* Do actual rootfinding */
  do {
    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, 1e-8);
    if (status)
      printf("Found eta:\n");
    printf("\titer=%d, eta=%e\n",iter, gsl_vector_get(s->x,0));
  }
  while (status == GSL_CONTINUE && iter < 1000);

  eta_fin = gsl_vector_get(s->x,0);

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(eta);

  return eta_fin;
}

/* Sweeps out a curve at the given gamma threshold in 
   eta vs. beta_perp space */

/* Starts at bperp_init. Each successive bperp is given by bperp_next =
   bperp_prev * factor -- thus, factor should typically be <1 as this is 
   the easier direction to move in. Loop ends at bperp_end. */

int curve_sweep(double bperp_init, double bperp_fin, double factor, 
		double eta_guess, char * filename, void * params) {

  double bperp;
  double eta_fin;

  Params_max * p = (Params_max *)params;
  FILE *outfile;

  /* First output to file */
  outfile = fopen(filename,"w");
  fprintf(outfile, "#bperp\teta\teperp\tk\tOr\tOi\n");
  fclose(outfile);

  bperp = bperp_init;
  eta_fin = eta_guess;
  if (bperp_init < bperp_fin) {
    if (factor <= 1.0) {
      /* doesn't make sense -- loop would never end */
      return 1;
    }
    else {
      do {
	p->bperp = bperp;
	printf("\nWorking on bperp = %e...\n",bperp);
	eta_fin = find_eta(eta_fin, p);

	outfile = fopen(filename,"a");
	fprintf(outfile, "%e\t%e\t%e\t%e\t%e\t%e\n", 
		bperp, eta_fin, p->eperp, kap_guess, Or_guess, Oi_guess);
	fclose(outfile);

	printf("Found eta = %e for bperp = %e!\n",eta_fin, bperp);
  
	bperp*=factor;
      } while (bperp < bperp_fin);
    }
  }
  else if (bperp_init >= bperp_fin) {
    if (factor >= 1.0) {
      /* doesn't make sense -- loop would never end */
      return 1;
    }
    else {
      do {
	p->bperp = bperp;

	printf("\nWorking on bperp = %e...\n",bperp);

	eta_fin = find_eta(eta_fin, p);
	
	outfile = fopen(filename,"a");
	fprintf(outfile, "%e\t%e\t%e\t%e\t%e\t%e\n", 
		bperp, eta_fin, p->eperp, kap_guess, Or_guess, Oi_guess);
	fclose(outfile);

	printf("Found eta = %e for bperp = %e!\n",eta_fin, bperp);

	bperp*=factor;
      } while (bperp >= bperp_fin);
    }
  }

  return 0;
  
}

int main(int argc, char *argv[]) {
  if (argc != 11) {
    printf("usage: %s <gamma> <eperp> <bperp_i> <bperp_f> <factor> <eta_guess> <k_guess> <Or_guess> <Oi_guess> <filename>", argv[0]);
  }
  else {
    double gam0, eperp, bperp_i, test_eps2;
    double factor, bperp_f, eta_guess;
    char filename[100];

    /* Parse command line: */
    gam0 = atof(argv[1]);
    eperp = atof(argv[2]);
    bperp_i = atof(argv[3]);
    bperp_f = atof(argv[4]);
    factor = atof(argv[5]);

    eta_guess = atof(argv[6]);
    kap_guess = atof(argv[7]);
    Or_guess = atof(argv[8]);
    Oi_guess = atof(argv[9]);
    strcpy(filename, argv[10]);

    test_eps2 = 1e-9;
  
    Params_max test_params = {
      .gam0 = gam0,
      .eperp = eperp,
      .bperp = bperp_i,
      .eps2 = test_eps2,
      .bessel = gsl_sf_bessel_Kn(2, 1./eperp)
    };

    curve_sweep(bperp_i, bperp_f, factor, eta_guess, filename, &test_params);
  }
  return 0;
}
