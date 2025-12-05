/* Copyright 2025 Måns Holmberg. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "twoface.h"

const double two_pi = 2.0 * M_PI;


void integratetransit(int m, int n, int k, double *planetx, double *planety, double *z, double p, double *r, double *f, double *spotx, double *spoty, double *spotradius, double *spotcontrast, double *planetangle, double *answer) {
/* Calculate integrated flux of a star if it is transited by a planet
  of radius p*R_star, at projected position (planetx, planety)
  in R_star units.
  Flux is normalized to out-of-transit flux.
  This algorithm works by integrating over concentric rings,
  the number of which is controlled by n. Use n=1000 for fair results.
  Planetx is the coordinate perpendicular to the transit chord
  normalized to stellar radius units, and planety is the one
  parallel to the transit chord, in a fashion such that it increases
  throughout the transit.
  We assume that the one-dimensional arrays spotx, spoty, spotradius
  and spotcontrast have the same length: the number of the spots.

  Input parameters:

  m             length of time series
  n             number of concentric rings
  k             number of spots
  planet[xy]    planetary center coordinates in stellar radii in sky-projected coordinate system [m]
  z             planetary center distance from stellar disk center in stellar radii     (cached) [m]
  p             planetary radius in stellar radii, scalar
  r             radii of integration annuli in stellar radii, non-decreasing            (cached) [n]
  f             2.0 * limb darkening at r[i] * width of annulus i                       (cached) [n]
  spotx, spoty  spot center coordinates in stellar radii in sky-projected coordinate system      [k]
  spotradius    spot radius in stellar radii                                                     [k]
  spotcontrast  spot contrast                                                                    [k]
  planetangle   value of [circleangle(r, p, z[i]) for i in xrange(m)]                 (cached) [m,n]

  (cached) means the parameter is redundant, and could be calculated from other parameters,
  but storing it and passing it to this routine speeds up iterative execution (fit or MCMC).
  Note that we do not take limb darkening coefficients, all we need is f.

  Output parameters:

  answer        model lightcurve, with oot=1.0                                                [m] */
  // Running indices for m, n, k, respectively.
  int M, N, K;
  // Out of transit flux to normalize lightcurve with.
  double ootflux;
  // Temporary storage for trapeze area to save on multiplications.
  double trapeze;
  // Temporary storage for what it is called.
  double spotcenterdistancesquared;
  // Cache for spot properties.
  double *spotcenterdistance, *spotangle;
  // Cache for trapezoid integration.
  double *values;
  // Cache for planet-spot center distance and central angle.
  double d, planetspotangle;
  /* Projected distance of center of star and center of spot, in stellar radius,
  accounting for the spot boundary plane being closer to the sphere radius
  than the tangent plane at spot center. An array of length k. */
  // If we have no spot:
  if (k == 0) {
    // Evaluate the integral for ootflux using the trapezoid method.
    ootflux = 0.0;
    for (N=0; N<n; N++) {
      ootflux += M_PI * *(r+N) * *(f+N);
    }
    for (M=0; M<m; M++) {
      // Transit?
      if (*(z+M) < 1.0 + p) {
        // Integrate over rings.
        *answer = 0.0;
        for (N=0; N<n; N++) {
          *answer += *(r+N) * (M_PI - *(planetangle + n*M + N)) * *(f+N);
        }
        // Normalize by trapezoid width and by ootflux.
        *answer /= ootflux;
      } else {
        *answer = 1.0;
      }
      // Advance pointer instead of using *(answer+i) all the times, for performance.
      answer++;
    }
    // Restore pointer.
    answer -= m;
  } else {
    //spotcenterdistance = malloc(k * sizeof(double));
    //spotangle = malloc(k * n * sizeof(double));
    //values = malloc(n * sizeof(double));
    // Instead, do a single malloc for speed.
    spotcenterdistance = malloc((k + k*n + n) * sizeof(double));
    spotangle = spotcenterdistance + k;
    values = spotcenterdistance + k * (n+1);
    // Loop over spots: fill up some arrays that do not depend on z.
    for (K=0; K<k; K++) {
      spotcenterdistancesquared = (*(spotx+K) * *(spotx+K) + *(spoty+K) * *(spoty+K)) * (1.0 - *(spotradius+K) * *(spotradius+K));
      *(spotcenterdistance+K) = sqrt(spotcenterdistancesquared);
      /* Calculate the half central angles of the spot, and store it in a single row of the 2D array.
      These values do not depend on z, that's why we cache them for all spots. */
      ellipseangle(r, *(spotradius+K), *(spotcenterdistance+K), n, spotangle + K*n);
    }
    // Evaluate the integral for ootflux using the trapezoid method.
    ootflux = 0.0;
    for (N=0; N<n; N++) {
      trapeze = M_PI;
      for (K=0; K<k; K++) {
        trapeze += (*(spotcontrast+K)-1.0) * *(spotangle + K*n + N);
      }
      ootflux += trapeze * *(r+N) * *(f+N);
    }
    // Loop over observation times and calculate answer.
    for (M=0; M<m; M++) {
      // Transit?
      if (*(z+M) < 1.0 + p) {
        /* The values array first stores the half arc length integrated contrast for no spots
        that we calculate in the next three lines. Later we add the spot contributions below,
        finally multipy by r*f and sum up. */
        for (N=0; N<n; N++) {
          *(values+N) = M_PI - *(planetangle + n*M + N);
        }
        // Cycle through spots and add their contributions.
        for (K=0; K<k; K++) {
          // Calculate distance of spot center and planet center for this moment.
          d = sqrt(pow(*(planetx+M) - *(spotx+K) * sqrt(1.0 - *(spotradius+K) * *(spotradius+K)), 2.0) + pow(*(planety+M) - *(spoty+K) * sqrt(1.0 - *(spotradius+K) * *(spotradius+K)), 2.0));
          // Calculate central angle between planet and spot.
          if ((*(spotcenterdistance+K) < 1e-12) || (*(z+M) < 1e-12)) {
            planetspotangle = 0.0;
          } else {
            double c = (pow(*(z+M), 2.0) + pow(*(spotcenterdistance+K), 2.0) - d*d) / (2.0 * *(z+M) * *(spotcenterdistance+K));
            if (c >  1.0) c =  1.0;
            if (c < -1.0) c = -1.0;
            planetspotangle = acos(c);
          }
          // Cycle through annuli.
          for (N=0; N<n; N++) {
            /* Evaluate the integrand on r[N].
            The geometry is described by planetspotangle (function of z(M))
            planetangle (function of r(N) and z(M), 2D array),
            and spotangle (function of r(N), does not depend on z, calculated above).
            For each value of r, there are five cases. */
            // Case 1: planet and spot arcs are disjoint, contributions add up.
            if (planetspotangle > *(planetangle + M*n + N) + *(spotangle + K*n + N)) {
              *(values+N) += (*(spotcontrast+K)-1.0) * *(spotangle + K*n + N);
            // Case 2: planet arc inside spot arc.
            } else if (*(spotangle + K*n + N) > planetspotangle + *(planetangle + n*M + N)) {
              *(values+N) += (*(spotcontrast+K)-1.0) * (*(spotangle + K*n + N) - *(planetangle + n*M + N));
            // Case 4: triangle inequality holds, partial overlap.
            } else if (*(planetangle + n*M + N) <= planetspotangle + *(spotangle + K*n + N)) {
              // Case 4a: partial overlap on one side only.
              if (2*M_PI - planetspotangle >= *(planetangle + n*M + N) + *(spotangle + K*n + N)) {
                *(values+N) += 0.5 * (*(spotcontrast+K)-1.0) * (*(spotangle + K*n + N) + planetspotangle - *(planetangle + n*M + N));
              // Case 4b: partial overlap on two sides.
              } else {
                *(values+N) += (*(spotcontrast+K)-1.0) * (M_PI - *(planetangle + n*M + N));
              }
            }
            // Case 3: planet arc covers spot arc, spot arc invisible. No need to do anything.
            //else
              //*(values+N) += 0.0;
          }
        }
        /* Now we multiply the half arc length integrated contrast by r*f
        to get the integrand, and sum it up right away. */
        *answer = 0.0;
        for (N=0; N<n; N++) {
          *answer += *(r+N) * *(f+N) * *(values+N);
        }
        // Finally, we normalize with ootflux.
        *answer /= ootflux;
      } else {
        // If not transit:
        *answer = 1.0;
      }
      // Advance pointer instead of using *(answer+i) all the times, for performance.
      answer++;
    }
    // Restore pointer.
    answer -= m;
    // Free spotcenterdistance, spotangle, and values, that live side by side.
    free(spotcenterdistance);
  }
  return;
}


void circleangle(double *r, double p, double z, int n, double *answer) {
/* Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside a circle of radius p (planet)
  with separation of centers z.
  This is a zeroth order homogeneous function, that is,
  circleangle(alpha*r, alpha*p, alpha*z) = circleangle(r, p, z).

  This version uses a binary search. It is only marginally faster
  than using direct comparisons: in the loop, we need to compare 
  i to a and b and n (or 0) n times. A direct loop with comparisons
  would be more expensive by one indirect addressing and by a double
  comparison instead of the integer one (assuming, of course, that the
  input array is sorted, and we only compare to the value that delimits
  the next case, not testing for all cases in each iteration). This is
  barely worth the overhead of the binary search, but is so much cooler.

  Input:
    r       array, non-negative, non-decreasing [n]
    p       scalar, non-negative
    z       scalar, non-negative
    n       number of elements

  Output:
    answer  one dimensional array [n] */
  /* For speed, we increment answer and r (possibly faster than
  adding i in each iteration), and restore it at the end. The compiler
  should not actually compile the restore operation if the value of
  the variable is not used any more. */
  /* If the circle arc of radius r is disjoint from the circular disk 
     of radius p, then the angle is zero. */
  int i, a, b;
  double zsquared = z*z;
  double psquared = p*p;
  if (p > z) {
    // Planet covers center of star.
    a = mybsearch(r, p-z, n);
    b = mybsearch(r, p+z, n);
    for(i=0; i<a; i++) {
      *answer = M_PI;
      answer++;
    }
    r += i;
    for(; i<b; i++) {
      *answer = acos(((*r)*(*r)+zsquared-psquared)/(2*z*(*r)));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else {
    // Planet does not cover center of star.
    a = mybsearch(r, z-p, n);
    b = mybsearch(r, z+p, n);
    for(i=0; i<a; i++) {
      *answer = 0.0;
      answer++;
    }
    r += i;
    for(; i<b; i++) {
      *answer = acos(((*r)*(*r)+zsquared-psquared)/(2*z*(*r)));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  }
  return;
}

void ellipseangle(double *r, double a, double z, int n, double *answer) {
/*  Calculate half central angle of the arc of circle of radius r
  (which concentrically spans the inside of the star during integration)
  that is inside an ellipse of semi-major axis a with separation of centers z.
  The orientation of the ellipse is so that the center of the circle lies on 
  the continuation of the minor axis. This is the orientation if the ellipse
  is a circle on the surface of a sphere viewed in projection, and the circle
  is concentric with the projection of the sphere.
  b is calculated from a and z, assuming projection of a circle of radius a
  on the surface of a unit sphere. If a and z are not compatible, a is clipped.
  This is not zeroth order homogeneous function, because it calculates b based
  on a circle of radius a living on the surface of the unit sphere.
  r is an array, a, and z are scalars. They should all be non-negative.
  We store the result on the n double positions starting with *answer.
  
  Input:

  r        radius of circle [n]
  a        semi-major axis of ellipse, non-negative
  z        distance between centers of circle and ellipse,
           non-negative and at most 1
  n        size of array a

  Output:

  answer   half central angle of arc of circle that lies inside ellipes [n]. */
  int i;
  // Degenerate case.
  if ((a<=0.0) || (z>=1.0)) {
    for (i=0; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else if (z<=0) {
  // Concentric case.
    int bound = mybsearch(r, a, n);
    for(i=0; i<bound; i++) {
      *answer = M_PI;
      answer++;
    }
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  } else if (a*a + z*z >= 1.0) {
  // Unphysical case: clip a. Now b=0.
    a = sqrt(1-z*z);
    int bound = mybsearch(r, z, n);
    for(i=0; i<bound; i++) {
      *answer = 0.0;
      answer++;
    }
    r += i;
    for(; i<n; i++) {
      *answer = acos(z/(*r));
      answer++;
      r++;
    }
    answer -= i;
    r -= i;
  } else {
  // Case of ellipse.
    double b = a * sqrt(1.0-z*z/(1-a*a));
    double zsquared = z*z;
    double asquared = a*a;
    // Calculate A based on a and z to mitigate rounding errors.
    //double A = pow(a/b,2.0) - 1.0;
    double A = zsquared / (1.0 - asquared - zsquared);
    /* First, go through all the cases where the ellipse covers C_r.
    If there is no such r, then bound=0, nothing happens. */
    int bound = mybsearch(r, b-z, n);
    for(i=0; i<bound; i++) {
      *answer = M_PI;
      answer++;
    }
    /* Now go through all the cases when C_r does not reach out to the ellipse.
    Again, the loop body might not get executed. */
    bound = mybsearch(r, z-b, n);
    for(; i<bound; i++) {
      *answer = 0.0;
      answer++;
    }
    /* Now take care of the cases where C_r and the ellipse intersect. z_crit = 1-a^2. */
    if (z < 1.0 - a*a)
      bound = mybsearch(r, z+b, n);
    else
      bound = n;
    r += i;
    for(; i<bound; i++) {
      // If not disjoint from the outside.
      // We must have y_+ > -b now.
      //yp = (-z + sqrt(z*z - A*((*r)*(*r) - zsquared - asquared)))/A;
      //*answer = acos((z - yp)/(*r));
      *answer = acos((z - (-z + sqrt(z*z - A*((*r)*(*r) - zsquared - asquared)))/A)/(*r));
      answer++;
      r++;
    }
    r -= i;
    for(; i<n; i++) {
      *answer = 0.0;
      answer++;
    }
    answer -= i;
  }
  return;
}

int mybsearch(double *array, double val, int n) {
/* Return the smallest of i = 0, 1, ..., n such that val < *(array+i),
with the convention that *(array+n) = inf. */
  if (val < *array)
    return 0;
  if (*(array+n-1) <= val)
    return n;
  int a = 0, b = n-1, c;
  /* From this point on, we always have
  0 <= a <= c <= b <= n-1,
  *(array+a) <= val, and val < *(array+b). */
  while (a+1 < b) {
    c = (a+b)/2;
    if (val < *(array+c))
      b = c;
    else
      a = c;
  }
  return b;
}


static double wrap_angle(double a) {
    a = fmod(a, two_pi);
    if (a < 0.0) a += two_pi;
    return a;
}

static void sort_angles(double *a, int n) {
    int i, j;
    for (i = 1; i < n; ++i) {
        double x = a[i];
        j = i - 1;
        while (j >= 0 && a[j] > x) {
            a[j + 1] = a[j];
            --j;
        }
        a[j + 1] = x;
    }
}

static int in_planet_twohalves(double dx, double dy, double p1_squared, double p2_squared, double cosphi, double sinphi) {
  double dist2 = dx*dx + dy*dy;
    if (dist2 <= p1_squared) {
        double dot1 = dx * cosphi + dy * sinphi;
        if (dot1 >= 0.0) return 1;
    }
    if (dist2 <= p2_squared) {
        double dot2 = dx * (-cosphi) + dy * (-sinphi);
        if (dot2 >= 0.0) return 1;
    }
    return 0;
}

static int in_spot_on_ring(double theta,
                           double spottheta,
                           double spot_halfangle)
{
    if (spot_halfangle <= 0.0) return 0;
    double dtheta = theta - spottheta;
    while (dtheta >  M_PI) dtheta -= 2.0 * M_PI;
    while (dtheta < -M_PI) dtheta += 2.0 * M_PI;
    return (fabs(dtheta) <= spot_halfangle);
}


static void get_halfdisk_bounds(double r, double z, double psi_p, double p, double phi, int *nbounds, double *bounds) {

  /* 1) Planet *circular* disk boundaries only for partial overlap */
  if (z > 0.0) {
    if ( (z > fabs(r - p)) && (z < r + p) && r > 0.0 ) {
      double k_c = (r*r + z*z - p*p) / (2.0 * r * z);
      if (k_c < -1.0) k_c = -1.0;
      if (k_c >  1.0) k_c =  1.0;
        double alpha_c = acos(k_c);
        *(bounds + (*nbounds)++) = psi_p - alpha_c;
        *(bounds + (*nbounds)++) = psi_p + alpha_c;
    }
  }
  /* 2) Planet half-plane boundaries: these must be handled
                  even if the full disk completely covers the ring */
  double m;
  if (z == 0.0) {
    /* Planet at star center: line through origin cuts every ring
                      at cos(theta - phi) = 0 -> m = 0 */
    m = 0.0;
  } else {
    double cos_delta = cos(psi_p - phi);
    m = (z / r) * cos_delta;
  }
  if (m > -1.0 && m < 1.0) {
    double beta = acos(m);
    *(bounds + (*nbounds)++) = phi - beta;
    *(bounds + (*nbounds)++) = phi + beta;
  }
}

static inline int angles_very_close(double a, double b, double eps) {
    double d = fabs(a - b);
    if (d > M_PI) d = fabs(d - 2.0*M_PI);
    return d <= eps;
}

static void get_limits(int n, double z, double dr, double r0, double pmax, int *Nmin, int *Nmax) {
  double rmin = z - pmax;
  if (rmin < 0.0) rmin = 0.0;
  double rmax = z + pmax;
  if (rmax > 1.0) rmax = 1.0;  // stellar radius is 1
  *Nmin = (int)ceil((rmin - r0) / dr);
  if (*Nmin < 0) *Nmin = 0;
  *Nmax = (int)floor((rmax - r0) / dr);
  if (*Nmax >= n) *Nmax = n - 1;
}


static double getE(double M, double e)	//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;
	double fe, fs;

	// modification from LK 05/07/2017:
	// add fmod to ensure convergence for diabolical inputs (following Eastman et al. 2013; Section 3.1)
	while(fmod(fabs(E - e*sin(E) - M), two_pi) > eps)
	{
		fe = fmod(E - e*sin(E) - M, two_pi);
		fs = fmod(1 - e*cos(E), two_pi);
		E = E - fe/fs;
	}
	return E;
}

void sky(int m, double *restrict t_array, double tc, double per, double a, double inc, double ecc, double omega, double *restrict z, double *restrict x, double *restrict y, double *restrict psi) {
  
  double n = two_pi / per;
  for(int i = 0; i < m; i++) {
	  double t = t_array[i];
		
		//calculates time of periastron passage from time of inferior conjunction
		double f = M_PI/2. - omega;								//true anomaly corresponding to time of primary transit center
		double E = 2.*atan(sqrt((1. - ecc)/(1. + ecc))*tan(f/2.));				//corresponding eccentric anomaly
		double M = E - ecc*sin(E);
		double tp = tc - per*M/2./M_PI;							//time of periastron

		if(ecc < 1.0e-5)
		{
			f = ((t - tp)/per - (int)((t - tp)/per))*two_pi;			//calculates f for a circular orbit
		}
		else
		{
			M = n*(t - tp);
			E = getE(M, ecc);
			f = 2.*atan(sqrt((1.+ecc)/(1.-ecc))*tan(E/2.));
		}

		double r = a * (1-ecc*cos(E));
		*(psi + i) = atan((-1/tan(omega+f))*cos(inc));
		*(y + i) = -r*sin(omega+f)*cos(inc);
		*(x + i) = -r*cos(omega+f);

		if (sin(f + omega)*sin(inc) <= 0.) *(z + i) = 100; //z < 0, so d is set to large value in order to not model primary transit during secondary eclipse
		else *(z + i) = a*(1.0 - ecc*ecc)/(1.0 + ecc*cos(f))*sqrt(1.0 - sin(omega + f)*sin(omega + f)*sin(inc)*sin(inc));	//calculates separation of centers
  }
}


void integratetransit_asymmetric(int m, int n, int k, double *restrict planetx, double *restrict planety, double *restrict z, double p, double p2, double *restrict r, double *restrict f, double *restrict spotx, double *restrict spoty, double *restrict spotradius, double *restrict spotcontrast, double *restrict planetphi, double *restrict answer, double *restrict spotcenterdistance, double *restrict bounds) {
  // Running indices for m, n, k, respectively.
  int M, N, K;
  // Out of transit flux to normalize lightcurve with.
  double ootflux, ootflux_inv;
  // Temporary storage for trapeze area to save on multiplications.
  double trapeze;
  // Temporary storage for what it is called.
  double spotcenterdistancesquared;
  // Cache for spot properties.
  // double *spotcenterdistance, *spotangle;
  double *spotangle;
  // Cache for trapezoid integration.
  // Cache for planet-spot center distance and central angle.
  // double d, planetspotangle;
  /* Projected distance of center of star and center of spot, in stellar radius,
  accounting for the spot boundary plane being closer to the sphere radius
  than the tangent plane at spot center. An array of length k. */

  double pmax = (p >= p2) ? p : p2;
  int Nmin, Nmax;
  double r0 = r[0];
  double dr = r[1] - r[0];   // uniform grid
  // double bounds[8+2*k]; // Maximum 6 spots
  // int nbounds = 0;
  double p_squared = p * p;
  double p2_squared = p2 * p2;

  spotangle = spotcenterdistance + k;
  // Loop over spots: fill up some arrays that do not depend on z.
  for (K=0; K<k; K++) {
    spotcenterdistancesquared = (*(spotx+K) * *(spotx+K) + *(spoty+K) * *(spoty+K)) * (1.0 - *(spotradius+K) * *(spotradius+K));
    *(spotcenterdistance+K) = sqrt(spotcenterdistancesquared);
    /* Calculate the half central angles of the spot, and store it in a single row of the 2D array.
    These values do not depend on z, that's why we cache them for all spots. */
    ellipseangle(r, *(spotradius+K), *(spotcenterdistance+K), n, spotangle + K*n);
  }
  // Evaluate the integral for ootflux using the trapezoid method.
  ootflux = 0.0;
  ootflux_inv = 0.0;
  for (N=0; N<n; N++) {
    trapeze = M_PI;
    for (K=0; K<k; K++) {
      trapeze += (*(spotcontrast+K)-1.0) * *(spotangle + K*n + N);
    }
    ootflux += trapeze * *(r+N) * *(f+N);
  }
  ootflux_inv = 1.0 / ootflux;

  // double *spottheta = malloc(k * sizeof(double));
  double delta_full[k];
  double spottheta[k];
  for (K = 0; K < k; ++K) {
    *(spottheta+K) = atan2(spoty[K], spotx[K]);
  }

  // Loop over observation times and calculate answer.
  for (M=0; M<m; M++) {
    double zM  = *(z + M);       /* distance star–planet center */
    // Transit?
    if (zM < 1.0 + pmax) {
      get_limits(n, zM, dr, r0, pmax, &Nmin, &Nmax);

      double x_p = *(planetx + M);
      double y_p = *(planety + M);
      double phi = *(planetphi + M);
      double cosphi = cos(phi);
      double sinphi = sin(phi);
      double psi_p = (zM > 0.0) ? atan2(y_p, x_p) : 0.0;

      double delta = 0.0;
      double value = 0.0;
      for (N = Nmin; N <= Nmax; N++) {

        double rN = *(r + N);
        int nbounds = 0;
        get_halfdisk_bounds(rN, zM, psi_p, p, phi, &nbounds, &bounds[0]);
        get_halfdisk_bounds(rN, zM, psi_p, p2, phi + M_PI, &nbounds, &bounds[0]);

        double total_angle = 0.0;
        double total_angle_planet = 0.0;

        // Spot bounds: add for all spots that intersect this ring and are not full-ring
        for (K = 0; K < k; ++K) {
            delta_full[K] = 0.0;
            double alpha_spot = spotangle[K*n + N];
            if (alpha_spot <= 0.0) continue;                    // no intersection
            if (alpha_spot >= M_PI - 1e-12) continue;           // full ring: no angular boundary

            // double th = spottheta[K];
             bounds[nbounds++] = spottheta[K] - alpha_spot;
             bounds[nbounds++] = spottheta[K] + alpha_spot;
        }

        if (nbounds == 0) {
            // No boundaries at all -> everything is uniform on the ring.
            double theta = 0.0;
            int inp = in_planet_twohalves(rN - x_p, - y_p,
                                              p_squared, p2_squared,
                                              cosphi, sinphi);

            if (inp) {
                  // Planet covers the entire ring
                  total_angle_planet = two_pi;
                  // No spot-visible-without-planet anywhere
            } else {
                // Planet covers none of the ring; spots are either full-ring or empty
                for (K = 0; K < k; ++K) {
                    double alpha_spot = spotangle[K*n + N];
                    if (alpha_spot <= 0.0) continue;
                    if (alpha_spot >= M_PI - 1e-12) {
                        // full ring in spot
                        delta_full[K] = two_pi;
                    } else {
                        int ins = in_spot_on_ring(theta, spottheta[K], alpha_spot);
                        if (ins) delta_full[K] = two_pi;
                    }
                }
            }
        } else {
            // 3) Wrap and sort unified boundaries
            for (int i = 0; i < nbounds; ++i) {
                bounds[i] = wrap_angle(bounds[i]);
            }
            sort_angles(bounds, nbounds);

            int nseg = nbounds;
            for (int i = 0; i < nseg; ++i) {
                double a1 = bounds[i];
                double a2 = bounds[(i + 1) % nseg];
                if (angles_very_close(a1, a2, 1e-12)) continue;

                double seg = a2 - a1;
                if (seg <= 0.0) seg += two_pi;

                double mid = a1 + 0.5 * seg;
                if (mid >= two_pi) mid -= two_pi;

                // Planet membership on this segment
                int inp = in_planet_twohalves(rN * cos(mid) - x_p, rN * sin(mid) - y_p,
                                                  p_squared, p2_squared,
                                                  cosphi, sinphi);
  
                if (inp) {
                    total_angle_planet += seg;
                } else {
                    // Planet NOT present: spots can contribute spot-visible-no-planet
                    for (K = 0; K < k; ++K) {
                        double alpha_spot = spotangle[K*n + N];
                        if (alpha_spot <= 0.0) continue;

                        int ins;
                        if (alpha_spot >= M_PI - 1e-12) {
                            ins = 1;  // full ring in spot
                        } else {
                            ins = in_spot_on_ring(mid, spottheta[K], alpha_spot);
                        }

                        if (ins) {
                            delta_full[K] += seg;
                        }
                    }
                }
            }
        }

        // 4) Combine into per-ring blocked value:
        //    value_N = -1/2 * total_angle_planet
        //              - Σ_K (c_K - 1) * spotangle[K,N]
        //              + Σ_K 1/2 (c_K - 1) * delta_full[K]
        double value = -0.5 * total_angle_planet;

        for (K = 0; K < k; ++K) {
            double alpha_spot = spotangle[K*n + N];
            if (alpha_spot <= 0.0) continue;

            double ck = spotcontrast[K];
            value -= (ck - 1.0) * alpha_spot;
            value += 0.5 * (ck - 1.0) * delta_full[K];
        }

        // 5) Accumulate flux contribution from this ring
        delta += r[N] * f[N] * value;
      }
      *answer = 1.0 + delta * ootflux_inv;
    } else {
      // If not transit:
      *answer = 1.0;
    }
    // Advance pointer instead of using *(answer+i) all the times, for performance.
    answer++;
  }
  // Restore pointer.
  answer -= m;
return;
}


