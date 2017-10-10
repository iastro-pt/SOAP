/******************************************************************************/
/*** file:   starspot.c                                                     ***/
/*** first version: 1.0 2006/08/29     X.Bonfils                            ***/
/***         Centro de Astronomia e Astrofisica da Universidade de Lisboa   ***/
/*** this version: 2.0 2014/12/05      X. Dumusque                          ***/
/***         Harvard-Smithsonian Center for Astrophysics                    ***/
/******************************************************************************/

#include <stdlib.h>
#include <starspot.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

FILE *file; /* Declare the file pointer */

struct gauss_params
{
    /*
     * Structure of the parameters of a Gaussian
     */
    double *x, *y, *err, *fit;
    int n, m;
};

int gauss_f(const gsl_vector *v, void *p, gsl_vector *f)
{
    /*
     * Function to calculate the residuals of a Gaussian fit
     */
    double c = gsl_vector_get(v, 0);
    double k = gsl_vector_get(v, 1);
    double x0 = gsl_vector_get(v, 2);
    double fwhm = gsl_vector_get(v, 3);
    
    double *x = ((struct gauss_params *) p)->x;
    double *y = ((struct gauss_params *) p)->y;
    double *err = ((struct gauss_params *) p)->err;
    double *fit = ((struct gauss_params *) p)->fit;
    int n = ((struct gauss_params *) p)->n;
    
    int i;
    double sigma = fwhm/2/sqrt(2*log(2));
    
    for (i = 0; i < n; i++) {
        fit[i] = c + k * exp(-(x[i]-x0)*(x[i]-x0)/2/sigma/sigma);
        gsl_vector_set (f, i, (fit[i]-y[i])/err[i]);
    }
    
    return GSL_SUCCESS;
}

int gauss_df(const gsl_vector *v, void *p, gsl_matrix *J)
{
    /*
     * Function to calculate the derivative of a Gaussian
     */
    double k = gsl_vector_get(v, 1);
    double x0 = gsl_vector_get(v, 2);
    double fwhm = gsl_vector_get(v, 3);
    
    double *x = ((struct gauss_params *) p)->x;
    double *err = ((struct gauss_params *) p)->err;
    int n = ((struct gauss_params *) p)->n;
    
    int i;
    double sigma = fwhm/2/sqrt(2*log(2));
    
    for (i = 0; i < n; i++) {
        double e = exp(-(x[i]-x0)*(x[i]-x0)/2/sigma/sigma);
        gsl_matrix_set (J, i, 0, 1/err[i]);
        gsl_matrix_set (J, i, 1, e/err[i]);
        gsl_matrix_set (J, i, 2, k/err[i]*e*(x[i]-x0)/sigma/sigma);
        gsl_matrix_set (J, i, 3, k/err[i]*e*(x[i]-x0)*(x[i]-x0)/sigma/sigma/sigma/2/sqrt(2*log(2)));
    }
    
    return GSL_SUCCESS;
}

int gauss_fdf(const gsl_vector *v, void *p, gsl_vector *f, gsl_matrix *J)
{
    /*
     * Function that calls the functions "gauss_f" and "gauss_df"
     */
    gauss_f (v, p, f);
    gauss_df (v, p, J);
    
    return GSL_SUCCESS;
}

void gauss_bis(double *vrad_ccf, double *ccf, double *err, int n1, double *mod, double *para, double *depth, double *bis, int len_depth)
{
    /*
     * Function to fit a Gaussian to the CCF and calculate the bisector of the CCF
     */
    struct gauss_params p;
    
    int i, j = 0, status, npar = 4;
    double c=0, k=0, vrad=0, fwhm=0;
    double sig_c=0, sig_k=0, sig_vrad=0, sig_fwhm=0;
    double *x;
    
    gsl_vector *v = gsl_vector_alloc (npar);
    gsl_matrix *covar = gsl_matrix_alloc (npar, npar);
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    gsl_multifit_function_fdf F;
    
    x = vrad_ccf;
    p.y = ccf;
    p.err = err;
    p.n = n1;
    
    p.x = (double *) malloc(n1*sizeof(double));
    p.fit = mod;
    
    for (i = 0; i < p.n; i++)
        p.x[i] = x[i]-x[0];
    
    c = (p.y[0]+p.y[p.n-1])/2;
    fwhm = (p.x[p.n-1]-p.x[0])/5.;
    k = 0.; vrad = p.x[p.n/2];
    
    for (i = 0; i < p.n; i++)
        if (abs(p.y[i]-c) > abs(k)) { k = p.y[i]-c; vrad = p.x[i]; }
    
    F.f = &gauss_f;
    F.df = &gauss_df;
    F.fdf = &gauss_fdf;
    F.n = p.n;
    F.p = npar;
    F.params = &p;
    
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, p.n, npar);
    gsl_vector_set(v, 0, c);
    gsl_vector_set(v, 1, k);
    gsl_vector_set(v, 2, vrad);
    gsl_vector_set(v, 3, fwhm);
    gsl_multifit_fdfsolver_set (s, &F, v);
    
    do {
        j++;
        status = gsl_multifit_fdfsolver_iterate (s);
        status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-5);
    }
    
    while (status == GSL_CONTINUE && j < 10000);
    
    gsl_multifit_covar (s->J, 0.0, covar);
    
    c = gsl_vector_get(s->x, 0); sig_c = sqrt(gsl_matrix_get(covar,0,0));
    k = gsl_vector_get(s->x, 1); sig_k = sqrt(gsl_matrix_get(covar,1,1));
    vrad = gsl_vector_get(s->x, 2); sig_vrad = sqrt(gsl_matrix_get(covar,2,2));
    fwhm = gsl_vector_get(s->x, 3); sig_fwhm = sqrt(gsl_matrix_get(covar,3,3));
    
    /* Bisector part : */
    // Declarations...
    double sigma = fwhm/2./pow(2.*log(2.),.5), vr, v0=vrad+x[0];
    double dCCFdRV, d2CCFdRV2, d2RVdCCF2;
    double *norm_CCF, *p0, *p1, *p2;
    
    // Allocations...
    norm_CCF = (double *)malloc(sizeof(double)*p.n);
    p0       = (double *)malloc(sizeof(double)*p.n);
    p1       = (double *)malloc(sizeof(double)*p.n);
    p2       = (double *)malloc(sizeof(double)*p.n);
    
    // Initialization...
    for (i=0; i<p.n; i++) norm_CCF[i] = -c/k*(1.-p.y[i]/c);
    
    for (i=0; i<p.n-1; i++) {
        vr = (x[i]+x[i+1])/2.;
        dCCFdRV = -(vr-v0)/pow(sigma,2)*exp(-pow((vr-v0),2)/2./pow(sigma,2));
        d2CCFdRV2 = (pow((vr-v0),2)/pow(sigma,2)-1)/pow(sigma,2)*exp(-pow((vr-v0),2)/2/pow(sigma,2));
        d2RVdCCF2 = -d2CCFdRV2/pow(dCCFdRV,3);
        p2[i] = d2RVdCCF2/2;
        p1[i] = (x[i+1]-x[i]-p2[i]*(pow(norm_CCF[i+1],2)-pow(norm_CCF[i],2)))/(norm_CCF[i+1]-norm_CCF[i]);
        p0[i] = x[i]-p1[i]*norm_CCF[i]-p2[i]*pow(norm_CCF[i],2);
    };//};
    
    int ind_max = 0, i_b, i_r;
    for (i=0; i<p.n; i++) if (norm_CCF[i]>norm_CCF[ind_max]) ind_max=i;
    for (i=0; i<len_depth; i++) {
        i_b = ind_max; i_r = ind_max;
        while ((norm_CCF[i_b] > depth[i]) && (i_b > 1)) i_b--;
        while ((norm_CCF[i_r+1] > depth[i]) && (i_r < (p.n-2))) i_r++;
        bis[i] = (p0[i_b]+p0[i_r]) + (p1[i_b]+p1[i_r])*depth[i] + (p2[i_b]+p2[i_r])*pow(depth[i],2);
        bis[i] /= 2.;
    };
    
    int n2=0, n3=0;
    double RV_top=0., RV_bottom=0., span=0.;
    for (i=0; i<len_depth; i++) {
        if ((depth[i]>=0.1) && (depth[i] <= 0.4)) {
            n2++;
			RV_top += bis[i];};
        if ((depth[i]>=0.6) && (depth[i] <= 0.9)) {
            n3++;
			RV_bottom += bis[i];};
    };
    RV_top    /= n2;
    RV_bottom /= n3;
    span = RV_top-RV_bottom;
	
    vrad = vrad+vrad_ccf[0];
    
    para[0] = c;
    para[1] = k;
    para[2] = vrad;
    para[3] = fwhm;
    para[4] = sig_c;
    para[5] = sig_k;
    para[6] = sig_vrad;
    para[7] = sig_fwhm;
    para[8] = span;
    
    gsl_vector_free(v);
    gsl_matrix_free(covar);
    gsl_multifit_fdfsolver_free(s);
    free(p.x); free(norm_CCF); free(p0); free(p1); free(p2);
}

double rndup(double n,int nb_decimal)
{
    /*
     * Function to round up a double type at nb_decimal
     */
    double t;
    t=n*pow(10,nb_decimal) - floor(n*pow(10,nb_decimal));
    if (t>=0.5)    
    {
        n*=pow(10,nb_decimal);//where n is the multi-decimal double
        n = ceil(n);
        n/=pow(10,nb_decimal);
    }
    else 
    {
        n*=pow(10,nb_decimal);//where n is the multi-decimal double
        n = floor(n);
        n/=pow(10,nb_decimal);
    }
    return n;
}    


double Delta_lambda(double line_width, double lambda_line0)
{
    /*
     * Calculates the broadening of a spectral line (or shift) taking into account 
     * the velocity and the wavelength. line_width in [m/s] and lambda_line0 in [Angstrom]
     */
    double c=299792458.;
    double line_spread=0.;
	double beta;
	
	// relativist case
	beta = line_width/c;
	line_spread = -1 * lambda_line0 * (1 - sqrt((1+beta)/(1-beta)));
	
    return line_spread;
}


void shifting_CCF(double *lambda, double *CCF, double *CCF_blueshifted, double v_shift,int n1)
{
    /*
     * Shift the CCF with a given velocity v_shift, doing a linear interpolation
     */
    int i;
    double *CCF_derivative;
	
    CCF_derivative = (double *) malloc(n1*sizeof(double));
    for (i=0;i<n1;i++) 
    {
		CCF_derivative[i] = 0.0;
        CCF_blueshifted[i] = 0.0;
    }
	
    // Calculates the CCF derivative
	if (lambda[1]-lambda[0] == 0) CCF_derivative[0] = 0;
	else CCF_derivative[0] = (CCF[1] - CCF[0])/(lambda[1]-lambda[0]);
	
	for (i=1;i<n1-1;i++)
	{
		if (lambda[i+1]-lambda[i] == 0) CCF_derivative[i] = CCF_derivative[i-1];
		else
		{
			CCF_derivative[i] = (CCF[i+1] - CCF[i])/(lambda[i+1]-lambda[i]);
		}
	}
	CCF_derivative[i] = 0.;
	
    // Calculates the new values of the shifted CCF doing a linear interpolation
    for (i=0;i<n1;i++)
    {
        if ((v_shift >= 0) || (i == 0)) 
        {
            CCF_blueshifted[i] = CCF[i] - v_shift * CCF_derivative[i]; // dy/dx = derivative => dy = derivative*dx
        }
        else if (v_shift < 0)
        {
            CCF_blueshifted[i] = CCF[i] - v_shift * CCF_derivative[i-1]; // dy/dx = derivative => dy = derivative*dx
        }
    }
	free(CCF_derivative);
}

// Calculates the Planck function
double loi_Planck(double lambda0, int Temp)
{
    double c   = 299792458.;     // speed of light m/s
    double h   = 6.62606896e-34; // Planck constant
    double k_b = 1.380e-23;      // Boltzmann constant
    return 2*h*pow(c,2)*1./pow(lambda0,5)*1./(exp((h*c)/(lambda0*k_b*Temp))-1);
}

void itot(double v, double i, double limba1, double limba2, double modif_bis_quad, double modif_bis_lin, double modif_bis_cte, int grid,
	  double *vrad_ccf, double *intensity_ccf, double v_interval, int n_v, int n,
	  double *f_star, double *sum_star)
{
    /*
     * Calculates the flux and CCF in each cell of the grid and integrate
     * over the entire stellar disc to have the integrated flux and CCF
     */

  double omega, delta_grid, delta_v;
  double y, z, delta, r_cos, limb;
  int iy, j, iz, diff_CCF_non_v_and_v,n_v_shifted_quotient;
  double n_v_shifted, n_v_shifted_remainder;
  double *intensity_ccf_shift;

  intensity_ccf_shift = (double *)malloc(sizeof(double)*n);

  /* Conversions */
  i = i * pi/180. ; // [degree]       --> [radian]

  omega = v;
  delta_grid = 2./grid; // step of the grid. grid goes from -1 to 1, therefore 2 in total
  // v_interval is from the velocity 0 to the edge of the spectrum taking into account minimal or maximal rotation (width - v to 0 or 0 to width + v). 
  // n_v is the number of points for all the CCF from minimum rotation to maximum one (from width - v to width + v).
  // n_v represent therefore the double than v_interval, we therefore have to multiply v_interval by 2.
  delta_v = 2.*v_interval/(n_v-1); //step in speed of the CCF. There is (n_v-1) intervals

  /* Calculates the total stellar intensity (without spots) */
  #if VERBOSE
  printf("Total stellar intensity \n");
  #endif

  // Scan of each cell on the grid
  for (iy=0; iy<=grid; iy++) // y-scan
  {
    y = -1. + iy*delta_grid; // y between -1 et 1
    delta = y * omega * sin(i);  // Give the velocity of the rotation for the given grid cell as a function of y.
                                 // For y=-1 => Vrot min, for y=1 => Vrot max, and for y=0 => Vrot = 0.

      for (iz=0; iz<=grid; iz++) // z-scan
    {
      z = -1. + iz*delta_grid;// z between -1 et 1
      if ((y*y+z*z)<=1.) //projected radius on the sky smaller than 1, which means that we are on the stellar disc
      {
          //limb-darkening
          r_cos = pow(1.-(y*y+z*z),.5); //sqrt(r^2-(y^2+z^2)) where r=1. This is equal to 1 at the disc center, and 0 on the limb.
                                        //This is often referred in the literature as cos(theta)
          limb =  1. - limba1*(1-r_cos) - limba2*(1-r_cos)*(1-r_cos); //intensity due to limb-darkening (Mandel & Agol 2002)
          
          //The CCF is defined between -20 and +20 km/s with a sampling "delta_v" of 0.1 km/s, and with boundaries normalized at 1. To account for stellar rotation
          //in each cell of the simulation, the CCF will be shifted depending on the position of the cell on the stellar disc. To shift the CCF by a RV of x km/s, we
          //will first shift the CCF by "n_v_shifted_quotient" = int(x/sampling) steps, and then interpolate the CCF with the remaining velocity "n_v_shifted_remainder"
          //to have the correct shift of x km/s at the end.
          
          //This will be done by shifting the CCF by
          // an integer of steps and then interpolating the CCF to match the correct velocity
          // The sampling of the CCF is every 100 m/s, so 2 km/s corresponds to a shift of 20 steps.
          n_v_shifted = delta/delta_v; // by how much steps the CCF is shifted due to rotation
          n_v_shifted_quotient = rndup(n_v_shifted,0);
          n_v_shifted_remainder = (delta - n_v_shifted_quotient*delta_v);
          
          double v_shift = n_v_shifted_remainder;
          //shifting the CCF with the remainder of n_v_shifted, the quotient will be taken into account by shifting all the points of the spectrum
          shifting_CCF(vrad_ccf, intensity_ccf, intensity_ccf_shift,v_shift,n);
          
          diff_CCF_non_v_and_v = (n_v - n) / 2.; //difference in number of steps between the CCF without any rotation and the one with rotation

          // To take into account the rotation, we increase the width of the CCF. Because the original CCF is only defined between -20 and +20 km/s,
          // we have to extrapolate on each side of the boundaries of the CCF
          for (j=0;j<diff_CCF_non_v_and_v+n_v_shifted_quotient;j++)
          {
              // extrapolation on the left of the CCF with the value of the left boundary, weighted by the limb-darkening
              f_star[j] += intensity_ccf_shift[0]*limb;
          }
          for (j=diff_CCF_non_v_and_v+n_v_shifted_quotient;j<n+(diff_CCF_non_v_and_v+n_v_shifted_quotient);j++)
          {
              // value of the CCF, weighted by the limb-darkening
              f_star[j] += (intensity_ccf_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)]) * limb;
          }
          for (j=n+(diff_CCF_non_v_and_v+n_v_shifted_quotient);j<n_v;j++)
          {
              // extrapolation on the right of the CCF with the value of the right boundary, weighted by the limb-darkening
              f_star[j] += intensity_ccf_shift[n-1]*limb;
          }
          // calculates the total flux taking into account the limb-darkening (intensity_ccf_shift[0] is very close to 1)
          *sum_star += intensity_ccf_shift[0]*limb;
      }
    }
  }
    // Free memory
    free(intensity_ccf_shift);
}

void starmap(double v, double i, double limba1, double limba2, int grid,
	double **Fmap, double **Vmap)
{
  /*
   * Calculates the flux F and velocity V maps of the star taking into account the rotational velocity, the stellar inclination and the limb-darkening
   * v     [km/s]          stellar rotation
   * i     [degree]        inclination of the rotational axis / sky plane
   * limba1 [0-1.]         linear coefficient of the quadratic limb-darkening law
   * limba1 [0-1.]         quadratic coefficient of the quadratic limb-darkening law
   * grid
   * Fmap  [arb. unit]     Flux map (2D array)
   * Vmap  [km/s]          Velocity map (2D array)
   */
	
  double r_cos, limb, y, z, delta_grid, delta_v;
  int iy, iz;
  v = v*sin(i/180.*pi);   //Projected rotational velocity
  delta_grid = 2./grid;   //The stellar disc goes from -1 to 1, therefore 2
  delta_v    = 2.*v/grid; //The stellar disc goes from -1 to 1, therefore 2
	
    for (iy=0; iy<grid; iy++) // y-axis scan...
    {
        y = -1 + iy*delta_grid;
    
        for (iz=0; iz<grid; iz++) // z-axis scan
        {
            z = -1 + iz*delta_grid;
		
            if ((y*y+z*z)<1) // If the
            {
                //limb-darkening
                r_cos = pow(1.-(y*y+z*z),.5); //sqrt(r^2-(y^2+z^2)) where r=1. This is equal to 1 at the disc center, and 0 on the limb.
                                              //This is often referred in the literature as cos(theta)
                limb =  1. - limba1*(1-r_cos) - limba2*(1-r_cos)*(1-r_cos); //intensity due to limb-darkening (Mandel & Agol 2002)
                Fmap[iy][iz]= limb;
                Vmap[iy][iz]= -v+iy*delta_v;
            }
            else
            {
                Fmap[iy][iz]=0;
                Vmap[iy][iz]=-9999;
            }
        }
    }
}

void spot_init(double s, double longitude, double latitude, double inclination, 
	       int nrho, double **xyz)
{
  /* Position of the spot initialized at the disc center
     (star rotates around the z axis) 

     s [spot radius] 
     longitude [degree]
     latitude  [degree]
     inclination [degree] i=0  -> pole-on (North)
                          i=90 -> equator-on
     nrho : Spot circumference resolution
     xyz  : real position of the spot
     xyz2 : position of the spot at the disc center */

  double *rho, rho_step, **xyz2;
  int j;

  /* Conversion [deg] -> [rad] */
  longitude   *= pi/180.;
  latitude    *= pi/180.;
  inclination *= pi/180.;

  // In this initial disc center position we calculate the coordinates (x,y,z) of
  // points of the active region's circumference
  rho = (double *)malloc(sizeof(double)*nrho);
  xyz2 = (double **)malloc(sizeof(double *)*nrho);
  for (j=0; j<nrho; j++) xyz2[j] = (double *)malloc(sizeof(double)*3);
  rho_step = 2.*pi/(nrho-1); //A circular active region has a resolution given by nrho, which implies that we
                             //will have a point on the disc circumference every 2*pi/(nrho-1). -1 because
                             //there is (nrho-1) intervals
  for (j=0; j<nrho; j++) rho[j] = -pi + j*rho_step;  //rho goes from -pi to pi
  for (j=0; j<nrho; j++) xyz2[j][0] = pow(1-s*s,.5); //x=sqrt(r^2-s^2), where r is the radius. r=1 therefore r^2=1
                                                     //The active region is on the surface, so very close to x=1.
                                                     //However with the curvature of the sphere, the circumference of
                                                     //the active region is at x=sqrt(r^2-s^2)
  for (j=0; j<nrho; j++) xyz2[j][1] = s*cos(rho[j]); //projection of the circumference on the y axis
  for (j=0; j<nrho; j++) xyz2[j][2] = s*sin(rho[j]); //projection of the circumference on the z axis

  // to account for the real projection of the spot, we rotate the star and look how the coordinates of the
  // circumference of the spot change
  // Position according to latitude, longitude and inclination. It consists of 
  // three rotations.
  //
  // Conventions :
  // -when inclination=0 the star rotates around z axis
  // -line of sight is along x-axis
  // -and sky plane = yz-plane
  //
  // Be Rx(alpha), Ry(beta), Rz(gamma) the rotations around the x, y and z axis
  // respectively, with angles alpha, beta and gamma (counter-clockwise 
  // direction when looking toward the origin).
  //
  // The rotations to apply are:
  //   Ry(inclination) x Rz(longitude) x Ry(latitude) x X(x,y,z)
  // 
  //         |  cos(b)  0  sin(b) |                       | cos(g) -sin(g) 0 |
  // Ry(b) = |    0     1    0    |               Rz(g) = | sin(g)  cos(g) 0 |
  //         | -sin(b)  0  cos(b) |                       |   0       0    1 |
  //
  //
  // |x'|   |  cos(b2)cos(g)cos(b)-sin(b)sin(b2)  -sin(g)cos(b2)  cos(b2)cos(g)sin(b)+sin(b2)cos(b) |   |x|
  // |y'| = |              sin(g)cos(b)               cos(g)                   sin(g)sin(b)         | x |y|
  // |z'|   |  -sin(b2)cos(g)cos(b)-cos(b2)sin(b)  sin(b2)sin(g) -sin(b2)cos(g)sin(b)+cos(b2)cos(b) |   |z|

  double b  = -latitude;
  double g = longitude;
  double b2 = pi/2.-inclination;

  double R[3][3] = {{cos(b2)*cos(g)*cos(b)-sin(b)*sin(b2), -sin(g)*cos(b2), cos(b2)*cos(g)*sin(b)+sin(b2)*cos(b)},
                    {sin(g)*cos(b),                          cos(g),          sin(g)*sin(b)},
                    {-sin(b2)*cos(g)*cos(b)-cos(b2)*sin(b), sin(b2)*sin(g), -sin(b2)*cos(g)*sin(b)+cos(b2)*cos(b)}};
  
  //calculates the real xyz position of the active region after rotating the star to have the initial equatorial active region
  //at the correct longitude and latitude, taking into account the stellar inclination
  for (j=0; j<nrho; j++) xyz[j][0] = R[0][0]*xyz2[j][0] + R[0][1]*xyz2[j][1] + R[0][2]*xyz2[j][2];
  for (j=0; j<nrho; j++) xyz[j][1] = R[1][0]*xyz2[j][0] + R[1][1]*xyz2[j][1] + R[1][2]*xyz2[j][2];
  for (j=0; j<nrho; j++) xyz[j][2] = R[2][0]*xyz2[j][0] + R[2][1]*xyz2[j][1] + R[2][2]*xyz2[j][2];

  // Free memory
  free(rho);
  for (j=0; j<nrho; j++) free(xyz2[j]);
  free(xyz2);
  
}


void spot_phase(double **xyz, double inclination, int nrho, double phase, 
		double **xyz2)
{
  int i;
    double psi = -phase*(2*pi); //the phase is between 0 and 1, so psi is in radian between -2pi and 0.
  inclination = inclination*pi/180.; // in radian

  double axe[3]  = {cos(inclination),0,sin(inclination)}; //projection of the rotation axis on the xyz coordinate
    
  //The rotation around the axis (axe[0],axe[1],axe[2]) of an angle psi is given by the following matrix
  //There is some sign difference here with respect to Wikipaedia for example because in this case psi is
  //defined negative (c.f. a few lines before). In this case, cos(-psi)=cos(psi) -> no sign change,
  //but sin(-psi)=-sin(psi) -> sign change.
  double R[3][3] = {{(1-cos(psi))*axe[0]*axe[0] + cos(psi), 
                     (1-cos(psi))*axe[0]*axe[1] + sin(psi)*axe[2],
                     (1-cos(psi))*axe[0]*axe[2] - sin(psi)*axe[1]},
                    {(1-cos(psi))*axe[1]*axe[0] - sin(psi)*axe[2],
                     (1-cos(psi))*axe[1]*axe[1] + cos(psi),
                     (1-cos(psi))*axe[1]*axe[2] + sin(psi)*axe[0]},
                    {(1-cos(psi))*axe[2]*axe[0] + sin(psi)*axe[1],
                     (1-cos(psi))*axe[2]*axe[1] - sin(psi)*axe[0],
                     (1-cos(psi))*axe[2]*axe[2] + cos(psi)}};

  for (i=0; i<nrho; i++) //calculates the xyz position of the active region circumference at phase psi
  {
    xyz2[i][0] = R[0][0]*xyz[i][0] + R[0][1]*xyz[i][1] + R[0][2]*xyz[i][2];
    xyz2[i][1] = R[1][0]*xyz[i][0] + R[1][1]*xyz[i][1] + R[1][2]*xyz[i][2];
    xyz2[i][2] = R[2][0]*xyz[i][0] + R[2][1]*xyz[i][1] + R[2][2]*xyz[i][2];
  }
}


int spot_area(double **xlylzl, int nrho, int grid, int *iminy, int *iminz, 
	       int *imaxy, int *imaxz)
{
  // Determine a smaller yz-area of the stellar disk, where the active region is.
  // The different cases are :
  // - the active region is completely visible (all x of the circumference >=0)
  // - the active region is completely invisible (all x of the circumference <0)
  // - the active region is on the disk edge and partially visible only
  int j, visible=0;
  double grid_step = 2./grid; //The stellar disc goes from -1 to 1, therefore 2
  double miny=1, minz=1, maxy=-1, maxz=-1; // init to 'opposite'-extreme values
  int counton=0, countoff=0; // count how many points of the circumference are
                             // visible and how many are invisible
  for (j=0; j<nrho; j++) //scan each point of the circumference
    if (xlylzl[j][0]>=0) { // if x>=0
      counton += 1;
      // select the extreme points of the circumference
      if (xlylzl[j][1]<miny) miny = xlylzl[j][1];
      if (xlylzl[j][2]<minz) minz = xlylzl[j][2]; 
      if (xlylzl[j][1]>maxy) maxy = xlylzl[j][1]; 
      if (xlylzl[j][2]>maxz) maxz = xlylzl[j][2];
    }
    else countoff = 1;

  if ((counton>0)&&(countoff>0)) { // There are both visible and invisible points
                                   // --> active region is on the edge
    // In this situation there are cases where the yz-area define above is 
    // actually smaller than the real area of the active region on the stellar disk.
    // The minima/maxima are over/under-estimated if the active region is on one of the 
    // axis (y or z). Because if on the y axis, the minimum (or maximum) won t be on the circumference of the active region. Same for z axis
    if (miny*maxy<0) {       //active region on the z-axis because one point is on the positive side of z, and the other on the negative side of z
      if (minz<0) minz=-1;   //active region on the bottom-z axis (z<0)
      else maxz=1;}          //active region on the top-z axis (z>=0)
    if (minz*maxz<0) {       //active region on the y-axis because one point is on the positive side of y, and the other on the negative side of z
      if (miny<0) miny=-1;   //active region on the left hand-y axis (y<0)
      else maxy=1;}          //active region on the right hand-y axis (y>=0)
  };
  if (counton==0) visible = 0;
  else visible = 1;
  
  // Indices of miny, minz,... on the grid
  *iminy = floor((1.+miny)/grid_step); //floor(x) returns the largest integral value that is not greater than x.
                                       //floor of 2.3 is 2.0, floor of 3.8 is 3.0, floor of -2.3 is -3.0, floor of -3.8 is -4.0
  *iminz = floor((1.+minz)/grid_step);
  *imaxy = ceil((1.+maxy)/grid_step);  //ceil(x) returns the smallest integral value that is not less than x.
                                       //ceil of 2.3 is 3, ceil of 3.8 is 4.0, ceil of -2.3 is -2.0, ceil of -3.8 is -3.0
  *imaxz = ceil((1.+maxz)/grid_step);
  
  return visible;
}

void spot_scan(double v, double i, double limba1, double limba2, double modif_bis_quad, double modif_bis_lin, double modif_bis_cte, int grid, 
          double *vrad_ccf, double *intensity_ccf, double *intensity_ccf_spot, double v_interval, int n_v, int n,
          double s, double longitude, double phase, double latitude, 
	      int iminy, int iminz, int imaxy, 
	      int imaxz, double *f_spot_flux, double *f_spot_bconv, double *f_spot_tot, double *sum_spot,
          int magn_feature_type, int T_star, int T_diff_spot)
{
    /* Scan of the yz-area where the spot is.
     * For each grid-point (y,z) we need to check whether it belongs to the spot
     * or not. Sadly, we do not know the projected geometry of the spot in its
     * actual position. Thus, we have to do an inverse rotation to replace the
     * grid point where it would be in the initial configuration. Indeed, in the
     * initial configuration, the spot has a well known geometry of a circle
     * centered on the x-axis.
     */
    int j, iy, iz,diff_CCF_non_v_and_v,n_v_shifted_quotient;
    double n_v_shifted, n_v_shifted_remainder;
    double y, z;
    double delta_grid=2./grid, delta, r_cos;
    double limb, delta_v = 2.*v_interval/(n_v-1);
    double *xayaza; // actual coordinates
    double *xiyizi; // coordinates transformed back to the initial configuration
    double *intensity_ccf_spot_shift;
    double *intensity_ccf_shift;
    int T_spot,T_plage;
	double intensity,loi_Planck_star;
    
    loi_Planck_star = loi_Planck(5293.4115e-10,T_star); //the wavelength of the Kitt peak spectrum goes from 3921.2441+6665.5789, the mean being 5293.4115
    
    xayaza             = (double *)malloc(sizeof(double)*3);
    xiyizi             = (double *)malloc(sizeof(double)*3);

    intensity_ccf_shift = (double *)malloc(sizeof(double)*n);
    intensity_ccf_spot_shift = (double *)malloc(sizeof(double)*n);
    
    // Scan of each cell on the grid
    for (iy=iminy; iy<imaxy; iy++) // y-scan
    {
        y = -1.+iy*delta_grid; // y between -1 et 1
        delta = y * v * sin(i*pi/180.);  // Give the velocity of the rotation for the given grid cell as a function of y.
                                         // For y=-1 => Vrot min, for y=1 => Vrot max, and for y=0 => Vrot = 0.
        xayaza[1] = y;

        for (iz=iminz; iz<imaxz; iz++) // z-scan
        {
            z = -1.+iz*delta_grid; // z between -1 et 1
            if (z*z+y*y<1.)  //projected radius on the sky smaller than 1, which means that we are on the stellar disc
            {
                xayaza[0] = pow(1.-(y*y+z*z),.5); //sqrt(r^2-(y^2+z^2)) where r=1. This is equal to 1 at the disc center, and 0 on the limb.
                                                  //This is often referred in the literature as cos(theta)
                xayaza[2] = z;

                // xayaza --> xiyizi: Rotate the star so that the spot is on the disc center
                spot_inverse_rotation(xayaza,longitude,latitude,i,phase,xiyizi); 
				
				// if inside the active region when scanning the grid
                if (xiyizi[0]*xiyizi[0]>=1.-s*s) // x^2 >= 1-s^2, which means that you are inside the active region
                { 
                    //limb-darkening
                    r_cos = pow(1.-(y*y+z*z),.5); //sqrt(r^2-(y^2+z^2)) where r=1. This is equal to 1 at the disc center, and 0 on the limb.
                                                  //This is often referred in the literature as cos(theta)
                    limb =  1. - limba1*(1-r_cos) - limba2*(1-r_cos)*(1-r_cos); //intensity due to limb-darkening (Mandel & Agol 2002)
                    
                    // intensity of the spot (magn_feature_type==0) or the plage (magn_feature_type==1)
                    if (magn_feature_type==0)
                    {
                        T_spot = T_star-T_diff_spot;
                        intensity = loi_Planck(5293.4115e-10,T_spot)/loi_Planck_star;  //the wavelength of the Kitt peak spectrum goes from 3921.2441+6665.5789, the mean being 5293.4115
                    }
                    else 
                    {
                        T_plage = T_star+250.9-407.7*r_cos+190.9*pow(r_cos,2); //plages are brighter on the limb Meunier 2010
                        intensity = loi_Planck(5293.4115e-10,T_plage)/loi_Planck_star; //the wavelength of the Kitt peak spectrum goes from 3921.2441+6665.5789, the mean being 5293.4115
                    }
                    
                    n_v_shifted = delta/delta_v; // by how much steps the CCF is shifted due to rotation
                    n_v_shifted_quotient = rndup(n_v_shifted,0); // integer number of steps
                    n_v_shifted_remainder = (delta - n_v_shifted_quotient*delta_v); // remainder of the division between delta and the integer number of steps
                        
                    double v_shift = n_v_shifted_remainder;
                    //shifting the CCF with the remainder of n_v_shifted, the quotient will be taken into account by shifting all the points of the spectrum
                    shifting_CCF(vrad_ccf, intensity_ccf, intensity_ccf_shift, v_shift,n);
                    shifting_CCF(vrad_ccf, intensity_ccf_spot, intensity_ccf_spot_shift, v_shift,n);
                    
                    diff_CCF_non_v_and_v = (n_v - n) / 2.; //difference in number of step between the CCF without any rotation and the one with rotation

                    // To take into account the rotation, we increase the width of the CCF. Because the original CCF is only defined between -20 and +20 km/s,
                    // we have to extrapolate on each side of the boundaries of the CCF. In this case, we calculate the "non contribution" to the CCF.
                    // For the region inside the active region, we calculate the contribution of the quiet photosphere and then suppress the contribution of the active region.
                    // It will then be easy to include the contribution of the active region by just subtracting the "non contribution" to the integrated contribution of the
                    // star without active region calculated with the "itot" function
                    for (j=0;j<diff_CCF_non_v_and_v+n_v_shifted_quotient;j++)
                    {
                        // extrapolation on the left of the CCF with the value of the left boundary, weighted by the limb-darkening and the active region intensity
                        // We also consider limb-darkening for spots because we also observe them with different stellar depth depending on their position
                        f_spot_flux[j]  += intensity_ccf_shift[0]*limb*(1 - intensity);// only flux effect
                        f_spot_bconv[j] += intensity_ccf_shift[0]*limb*(1 - 1);        // only convective blueshift effect. The convective blueshift does not affect the boundaries of the CCF
                        f_spot_tot[j]   += intensity_ccf_shift[0]*limb*(1 - intensity);// combined effect
                    }
                    for (j=diff_CCF_non_v_and_v+n_v_shifted_quotient;j<n+(diff_CCF_non_v_and_v+n_v_shifted_quotient);j++)
                    {
                        // value of the CCF, weighted by the limb-darkening and the active region intensity
                        // We also consider limb-darkening for spots because we also observe them with different stellar depth depending on their position
                        f_spot_flux[j]  += intensity_ccf_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)]*limb - intensity * (intensity_ccf_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)])*limb;     // only flux effect
                        f_spot_bconv[j] += intensity_ccf_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)]*limb - (intensity_ccf_spot_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)])*limb;           // only convective blueshift effect
                        f_spot_tot[j]   += intensity_ccf_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)]*limb - intensity * (intensity_ccf_spot_shift[j-(diff_CCF_non_v_and_v+n_v_shifted_quotient)])*limb; // combined effect
                    }
                    for (j=n+(diff_CCF_non_v_and_v+n_v_shifted_quotient);j<n_v;j++)
                    {
                        // extrapolation on the right of the CCF with the value of the right boundary, weighted by the limb-darkening and the active region intensity
                        // We also consider limb-darkening for spots because we also observe them with different stellar depth depending on their position
                        f_spot_flux[j]  += intensity_ccf_shift[n-1]*limb*(1 - intensity);// only flux effect
                        f_spot_bconv[j] += intensity_ccf_shift[n-1]*limb*(1 - 1);        // only convective blueshift effect. The convective blueshift does not affect the boundaries of the CCF
                        f_spot_tot[j]   += intensity_ccf_shift[n-1]*limb*(1 - intensity);// combined effect
                    }
                    
                    // calculates the "non contributing" total flux of the active region taking into account the limb-darkening and the active region intensity
                    *sum_spot += intensity_ccf_shift[0]*limb*(1.-intensity);
                }
            }
        }
    }
    free(xayaza); free(xiyizi);
    free(intensity_ccf_shift); free(intensity_ccf_spot_shift);
}

void spot_inverse_rotation(double *xyz, double longitude, double latitude, 
			   double inclination, double phase, double *xiyizi)
{
    /*
    * Relocate a point (x,y,z) to the 'initial' configuration
    * i.e. when the active region is on the disc center
    *
    * Thus it consists of rotating the point, according to latitude, longitude,
    *  inclination and phase, but in the reverse order.
    */

    //
    // Conventions :
    // -when inclination=0 the star rotates around z axis (i.e. rotation axis and
    // z axis are indistinct), 
    // -line of sight is along x-axis
    // -and sky plane = yz-plane
    //

    double g2 = --phase*(2*pi);  // inverse phase ([0-1] -> [rad])
    double i = inclination * pi/180.;

    double b  =  latitude  * pi/180.;
    double g  = -longitude * pi/180.;
    double b2 = -(pi/2.-i);
    
    double R[3][3]  = {{(1-cos(g2))*cos(i)*cos(i) + cos(g2), sin(g2)*sin(i),  (1-cos(g2))*cos(i)*sin(i)},
                       {-sin(g2)*sin(i),                     cos(g2),         sin(g2)*cos(i)},
                       {(1-cos(g2))*sin(i)*cos(i),           -sin(g2)*cos(i), (1-cos(g2))*sin(i)*sin(i) + cos(g2)}};
    
    double R2[3][3] =  {{cos(b)*cos(g)*cos(b2)-sin(b2)*sin(b),  -sin(g)*cos(b), cos(b)*cos(g)*sin(b2)+sin(b)*cos(b2)},
                        {sin(g)*cos(b2),                        cos(g),         sin(g)*sin(b2)},
                        {-sin(b)*cos(g)*cos(b2)-cos(b)*sin(b2), sin(b)*sin(g),  -sin(b)*cos(g)*sin(b2)+cos(b)*cos(b2)}};
    
    //rotation for the latitude, longitude, inclination and phase, which is a combination of R and R2
    double R3[3][3] = {{R2[0][0]*R[0][0]+R2[0][1]*R[1][0]+R2[0][2]*R[2][0],
                      R2[0][0]*R[0][1]+R2[0][1]*R[1][1]+R2[0][2]*R[2][1],
                      R2[0][0]*R[0][2]+R2[0][1]*R[1][2]+R2[0][2]*R[2][2]},
                     {R2[1][0]*R[0][0]+R2[1][1]*R[1][0]+R2[1][2]*R[2][0], 
                      R2[1][0]*R[0][1]+R2[1][1]*R[1][1]+R2[1][2]*R[2][1],
                      R2[1][0]*R[0][2]+R2[1][1]*R[1][2]+R2[1][2]*R[2][2]},
                     {R2[2][0]*R[0][0]+R2[2][1]*R[1][0]+R2[2][2]*R[2][0], 
                      R2[2][0]*R[0][1]+R2[2][1]*R[1][1]+R2[2][2]*R[2][1],
                      R2[2][0]*R[0][2]+R2[2][1]*R[1][2]+R2[2][2]*R[2][2]}};
    
    //Coordinates of the active region circumference when it is in the disc center
    xiyizi[0] = R3[0][0]*xyz[0] + R3[0][1]*xyz[1] + R3[0][2]*xyz[2];
    xiyizi[1] = R3[1][0]*xyz[0] + R3[1][1]*xyz[1] + R3[1][2]*xyz[2];
    xiyizi[2] = R3[2][0]*xyz[0] + R3[2][1]*xyz[1] + R3[2][2]*xyz[2];

}

void spot_scan_npsi(double **xyz, int nrho, double *psi, int npsi, double v, 
                   double inclination, double limba1, double limba2, double modif_bis_quad, double modif_bis_lin, double modif_bis_cte, int grid, 
				   double *vrad_ccf, double *intensity_ccf, double *intensity_ccf_spot, double v_interval, int n_v, int n,
				   double s, double longitude, double latitude, 
	      		   double **f_spot_flux, double **f_spot_bconv, double **f_spot_tot, double *sum_spot,
                   int magn_feature_type, int T_star, int T_diff_spot)
{
    /* 
    * Scans the yz-area where the spot is for different phases (psi) and        
    * returns the spot's "non-contribution" to the total flux and its           
    * "non-contribution" to the ccf, for each phase.
    * Thus the result is to be subtracted to the output of the itot() function.
    */

    //tbd before:
    //spot_init(s, longitude, latitude, inclination, nrho, xyz) #out: xyz
    int ipsi, j;
    int iminy, iminz, imaxy, imaxz, vis;
    double **xyz2 = (double **)malloc(sizeof(double *)*nrho);
    for (j=0; j<nrho; j++) xyz2[j] = (double *)malloc(sizeof(double)*3);

    for (ipsi=0; ipsi<npsi; ipsi++) 
    {
        spot_phase(xyz, inclination, nrho, psi[ipsi], xyz2);
        vis = spot_area(xyz2, nrho, grid, &iminy, &iminz, &imaxy, &imaxz);
        if (vis==1)
        {
            spot_scan(v, inclination, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, grid, 
                    vrad_ccf, intensity_ccf, intensity_ccf_spot,v_interval, n_v, n,
                    s, longitude, psi[ipsi], latitude, iminy, iminz, imaxy, 
                    imaxz, f_spot_flux[ipsi], f_spot_bconv[ipsi], f_spot_tot[ipsi], &sum_spot[ipsi],magn_feature_type,T_star,T_diff_spot);
        }  
    }
    for (j=0; j<nrho; j++) free(xyz2[j]);
    free(xyz2);
}
