/*! \file
    \brief Lagrangian of the Curved N Body Problem (CNBP).

    $Author: roldan $
    $Date: $
*/

#include <math.h>       // M_PI
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // memcpy
#include <assert.h>
#include <gsl/gsl_errno.h>  // GSL_SUCCESS

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>		// To compute determinant of matrix
#include <gsl/gsl_eigen.h>		// To compute eigenvalues of matrix
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// number of masses in the configuration
#define N 4	

/** Lagrangian of the CNBP */
double Lagr(double z[2*N], double kappa, double m[N]);

/** Kinetic energy */
double kinetic(double z[2*N], double kappa, double m[N]);

/** Conformal factor */
double conf(double z[2], double kappa);

/** Norm square */
double normsq(double z[2]);

/** Potential energy */
double potential(double z[2*N], double kappa, double m[N]);
double V(double *zj, double *zk, double kappa);

/** Norm */
double norm(double z[2]);

/** Gradient of Lagrangian L w.r.t. positions u */
void gradL(double u[2*N], double K, double m[N], double grad[2*N]);

/** Partial derivative of T (kinetic energy) w.r.t. x_k */
double pd_T_xk(double u[2*N], int k, double K, double m[N]);
double pd_T_yk(double u[2*N], int k, double K, double m[N]);

/** Partial derivative of U (potential energy) w.r.t. x_k */
double pd_U_xk(double u[2*N], int k, double K, double m[N]);
double pd_U_yk(double u[2*N], int k, double K, double m[N]);

double pd_V_xk(double xj, double yj, double xk, double yk, double K);
double pd_V_yk(double xj, double yj, double xk, double yk, double K);

/** Dot (inner) product of u and v */
double dot(double u[2], double v[2]);

/** Print eigenvalues of matrix A **/
void eigenvalues(gsl_matrix *A, int inPlace);

/** Determinant of matrix A **/
double det_get(gsl_matrix *A, int inPlace);

struct grad_params
{
	double K;		// curvature \kappa
	double m[N];	// masses m_1, m_2, ..., m_N
	double a[2*N];	// RE
};

int grad_f(const gsl_vector *x, void *params, gsl_vector *f)
{
	double K = ((struct grad_params *)params)->K;
	double m[N];
	double a[2*N];

	for(int i=0; i<N; i++)
		m[i] = ((struct grad_params *)params)->m[i];

	for(int i=0; i<2*N; i++)
		a[i] = ((struct grad_params *)params)->a[i];

	double u[2*N];

	for(int i=0; i<2*N; i++)
	{
		u[i] = gsl_vector_get(x, i);
	}

	double alpha = gsl_vector_get(x,2*N);	// Lagrange multiplier

	double grad[2*N];
	gradL(u, K, m, grad);

	for(int i=0; i<N; i++)
	{
		int k0 = 2*i;
		int k1 = 2*i+1;
		gsl_vector_set(f, k0, grad[k0] - alpha*u[k1]);
		gsl_vector_set(f, k1, grad[k1] + alpha*u[k0]);
	}
	double uJa = 0;
	for(int i=0;i<N; i++)
	{
		int k0 = 2*i;
		int k1 = 2*i+1;
		uJa = uJa -u[k0]*a[k1] + u[k1]*a[k0];
	}
	gsl_vector_set(f, 2*N, uJa);
	return GSL_SUCCESS;
}

int grad_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double K = ((struct grad_params *)params)->K;
	double m[N];
	double a[2*N];
	double det;	// determinant of Jac(grad_L)

	int n = 2*N+1;

	for(int i=0; i<N; i++)
		m[i] = ((struct grad_params *)params)->m[i];

	for(int i=0; i<2*N; i++)
		a[i] = ((struct grad_params *)params)->a[i];

	gsl_multiroot_function f = {&grad_f, n, params};

	gsl_vector *f_vec = gsl_vector_alloc(n);
	grad_f(x,params,f_vec);
	gsl_multiroot_fdjacobian (&f, x, f_vec, GSL_SQRT_DBL_EPSILON, J);

	for(int i=0; i<N; i++)
	{
		int k0 = 2*i;
		int k1 = 2*i+1;
		gsl_matrix_set(J, k0, 2*N, -gsl_vector_get(x, k1));
		gsl_matrix_set(J, k1, 2*N,  gsl_vector_get(x, k0));
	}
	gsl_matrix_set(J, 2*N, 2*N,  0.0);

	fprintf(stderr, "JACOBIAN:\n");
	for(int i=0; i<n; i++)
	{
		for(int k=0; k<n; k++)
			fprintf(stderr, "% .14e ", gsl_matrix_get(J,i,k));
		fprintf(stderr, "\n");
	}

	// Compute determinant (in place).
	// det = det_get(J, 1);	
	// fprintf(stderr, "det(JACOBIAN): %f\n", det);

	// Compute eigenvalues of A (not in place)
	fprintf(stderr, "eigenvals(JACOBIAN):\n");
	eigenvalues(J, 0);

	return(GSL_SUCCESS);
}

int grad_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	grad_f(x, params, f);
	grad_df(x, params, J);

	return(GSL_SUCCESS);
}

int
print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
	fprintf (stderr, "iter = %3lu\n", iter);

	fprintf (stderr, "x = ");
	for(int i=0; i<=2*N; i++)
		fprintf (stderr, "% .5f ", gsl_vector_get(s->x, i));
	fprintf (stderr, "\n");

	fprintf (stderr, "f(x) = ");
	for(int i=0; i<=2*N; i++)
		fprintf (stderr, "% .5f ", gsl_vector_get(s->f, i));
	fprintf (stderr, "\n\n");
}

int main( )
{
	// Cuadrilateral RE with one smaller mass than the other three
	double m[N] = {0.5, 1, 1, 1};

	// positions (x_j, y_j) for each of the N bodies.

	double a[2*N] = {-0.00000000000000, 0.56884869024927, 
		0.00000000000000, -0.46321826526006, 
		-0.42801097159491, 0.08939696006771,
		0.42801097159491, 0.08939696006771};

	double kappa;	// curvature

	/*
	double grad[6];
	gradL(z,kappa,m,grad);

	printf("Gradient: %f %f %f %f %f %f\n", 
			grad[0], grad[1], 
			grad[2], grad[3], 
			grad[4], grad[5]);

	double deltaAlpha = 3.0/100;
	for(int i=1; i<=100; i++)
	{
		alpha = i*deltaAlpha;

		// positions (x_j, y_j) for each of the 3 bodies.
		double z[6] = {alpha, 0, -alpha/2, sqrt(3)*alpha/2, -alpha/2,
			-sqrt(3)*alpha/2};

		printf("%f %f\n", alpha, Lagr(z,kappa,m));
	}
	*/

	int num_K = 20000;
	double delta_K = -200.0/num_K;

	const size_t n = 2*N+1;
	int i;
	double x_init[n];
	for(i=0; i<2*N; i++)
		x_init[i] = a[i];

	x_init[2*N] = 0;	// alpha (Lagrange multiplier)

	for(int j=0; j<num_K; j++)
	{
		kappa = j*delta_K;

		const gsl_multiroot_fdfsolver_type *T;
		gsl_multiroot_fdfsolver *s;	

		int status;
		size_t iter = 0;

		struct grad_params p; 
		p.K = kappa;
		for(int i=0; i<N; i++)
			p.m[i] = m[i];
		for(int i=0; i<2*N; i++)
			p.a[i] = a[i];

		gsl_multiroot_function_fdf f = {&grad_f, &grad_df, &grad_fdf, n, &p};

		gsl_vector *x = gsl_vector_alloc (n);

		for(i=0; i<n; i++)
			gsl_vector_set (x, i, x_init[i]);

		T = gsl_multiroot_fdfsolver_hybridsj;
		s = gsl_multiroot_fdfsolver_alloc (T, n);
		gsl_multiroot_fdfsolver_set (s, &f, x);

		//print_state (iter, s);

		do
			{
			  iter++;
			  status = gsl_multiroot_fdfsolver_iterate (s);


			  if (status)   /* check if solver is stuck */
				break;

			  status =
				gsl_multiroot_test_residual (s->f, 1e-12);


			  print_state (iter, s);
			}
		while (status == GSL_CONTINUE && iter < 1000);

		fprintf (stderr, "status = %s\n", gsl_strerror (status));

		if (status)   /* check if solver is stuck */
		{
			print_state (iter, s);
			exit(EXIT_FAILURE);
		}

		// Output relative equilibrium
		printf ("% .3f ", kappa);
		for(i=0; i<=2*N; i++)
			printf("% .14f ", gsl_vector_get(s->x, i));
		printf("\n");

		// Use current root as initial guess for next curvature
		for(i=0; i<n; i++)
			x_init[i] = gsl_vector_get (s->x, i);

		gsl_multiroot_fdfsolver_free (s);
		gsl_vector_free (x);
	}

	exit(EXIT_SUCCESS);
}

/** 
  Lagrangian of the CNBP in rotating coordinates.

  The expression of the Lagrangian given here is valid for both 
  positive and negative curvature kappa (as well as for the 
  Euclidean case kappa = 0).

  The position (x_j, y_j) of the bodies is parametrized by the 
  stereographic projection.

  \param[in] z
  positions (x_j, y_j) for each of the 3 bodies.

  \param[in] kappa
  curvature.

  \param[in] m
  masses m_j for each of the 3 bodies

  \returns Value of the Lagrangian of the system.
  */

double Lagr(double z[2*N], double kappa, double m[N])
{
	return kinetic(z,kappa,m) + potential(z,kappa,m);
}

/** Kinetic energy */
double kinetic(double z[2*N], double kappa, double m[N])
{
	int j;

	double T = 0;
	for(j=0; j<N; j++)
	{
		T += m[j]*conf(z+2*j, kappa)*normsq(z+2*j);
	}
	return 0.5*T;
}

/** Conformal factor */
double conf(double z[2], double kappa)
{
	double temp = 1 + kappa*normsq(z);
	return 4/(temp*temp);
}

/** Norm square */
double normsq(double z[2])
{
	double x=z[0];
	double y=z[1];

	return x*x+y*y;
}

/** Potential energy */
double potential(double z[2*N], double kappa, double m[N])
{
	int j, k;

	double U = 0;
	for(k=0; k<N; k++)
	{
		for(j=0; j<k; j++)
		{
			U += m[k]*m[j]*V(z+2*j, z+2*k, kappa);
		}
	}
	return U;
}

double V(double *zj, double *zk, double kappa)
{
	double xj = *zj;
	double yj = *(zj+1);
	double xk = *zk;
	double yk = *(zk+1);

	double result;

	if(kappa != 0)
	{
		double dot = xj*xk+yj*yk;	/* dot product z_j \cdot z_k */
		double tempj = normsq(zj)*kappa - 1;
		double tempk = normsq(zk)*kappa - 1;

		double temp2j = normsq(zj)*kappa + 1;
		double temp2k = normsq(zk)*kappa + 1;

		double temp2jsq = temp2j*temp2j;
		double temp2ksq = temp2k*temp2k;

		double num = 4*dot*kappa + tempk*tempj;
		double den = sqrt(temp2ksq*temp2jsq/kappa - num*num/kappa);

		result = num/den;
	}
	else // kappa == 0
	{
		double d[2] = {xj-xk, yj-yk};
		result = 1.0/(2*norm(d));
	}
	return result;
}

/** Norm */
double norm(double z[2])
{
	double x=z[0];
	double y=z[1];

	return sqrt(x*x+y*y);
}

/** 
  Gradient of Lagrangian L w.r.t. positions u.

  \param[in] u
  positions (x_j, y_j) for each of the N bodies.

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the N bodies

  \param[out] grad[2*N]
  Value of the gradient evaluated at u.

  \returns 
  */

void gradL(double u[2*N], double K, double m[N], double grad[2*N])
{
	for(int i=0; i<N; i++)
	{
		grad[2*i  ] = pd_T_xk(u,i,K,m) + pd_U_xk(u,i,K,m);	// w.r.t. x_k
		grad[2*i+1] = pd_T_yk(u,i,K,m) + pd_U_yk(u,i,K,m);	// w.r.t. y_k
	}
}

/** 
  Partial derivative of T (kinetic energy) w.r.t. x_k

  \param[in] u
  positions (x_j, y_j) for each of the N bodies.

  \param[in] k
  differentiate w.r.t. x_k (0 <= k <= N-1)

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the N bodies

  \returns Value of the partial derivative evaluated at u.
  */

double pd_T_xk(double u[2*N], int k, double K, double m[N])
{
	assert(k >= 0 && k <= N-1);
	double *uk = u + 2*k;
	double mk = m[k];
	double xk = uk[0];

	return mk*conf(uk, K)*((-2.0*K*xk)/(1+K*normsq(uk))*normsq(uk) + xk);
}

double pd_T_yk(double u[2*N], int k, double K, double m[N])
{
	assert(k >= 0 && k <= N-1);
	double *uk = u + 2*k;
	double mk = m[k];
	double yk = uk[1];

	return mk*(conf(uk, K)*(-2.0*K*yk)/(1+K*normsq(uk))*normsq(uk) + 
				conf(uk, K)*yk);
}

/** 
  Partial derivative of U (potential energy) w.r.t. x_k

  \param[in] u
  positions (x_j, y_j) for each of the N bodies.

  \param[in] k
  differentiate w.r.t. x_k (0 <= k <= N-1)

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the N bodies

  \returns Value of the partial derivative evaluated at u.
  */

double pd_U_xk(double u[2*N], int k, double K, double m[N])
{
	double *uk = u + 2*k;
	double xk = uk[0];
	double yk = uk[1];

	double sum = 0;
	for(int j=0; j<N; j++)
	{
		if (j!=k)
		{
			double *uj = u + 2*j;
			double xj = uj[0];
			double yj = uj[1];
			sum = sum + m[k]*m[j]*pd_V_xk(xj,yj,xk,yk,K);
		}
	}
	return sum;
}

double pd_U_yk(double u[2*N], int k, double K, double m[N])
{
	double *uk = u + 2*k;
	double xk = uk[0];
	double yk = uk[1];

	double sum = 0;
	for(int j=0; j<N; j++)
	{
		if (j!=k)
		{
			double *uj = u + 2*j;
			double xj = uj[0];
			double yj = uj[1];
			sum = sum + m[k]*m[j]*pd_V_yk(xj,yj,xk,yk,K);
		}
	}
	return sum;
}

double pd_V_xk(double xj, double yj, double xk, double yk, double K)
{
	double ujuk = xj*xk + yj*yk;

	double ujsq = xj*xj + yj*yj;
	double uksq = xk*xk + yk*yk;

	double tempj = ujsq*K-1.0;
	double tempk = uksq*K-1.0;

	double dx = xj-xk;
	double dy = yj-yk;

	double A = dx*dx + dy*dy;
	double sqrtA = sqrt(A);

	double B = uksq*ujsq*K*K + 2*ujuk*K + 1.0;
	double sqrtB = sqrt(B);

	double num = 4*ujuk*K + tempk*tempj;
	double den = 2*sqrtA*sqrtB;

	double pd_num_xk = 4*xj*K + 2*xk*K*tempj;
	double pd_den_xk = -2.0*dx/sqrtA*sqrtB + 
		sqrtA*(2*xk*ujsq*K*K + 2*xj*K)/sqrtB;

	return (pd_num_xk*den - num*pd_den_xk)/(den*den);
}

double pd_V_yk(double xj, double yj, double xk, double yk, double K)
{
	double ujuk = xj*xk + yj*yk;

	double ujsq = xj*xj + yj*yj;
	double uksq = xk*xk + yk*yk;

	double tempj = ujsq*K-1.0;
	double tempk = uksq*K-1.0;

	double dx = xj-xk;
	double dy = yj-yk;

	double A = dx*dx + dy*dy;
	double sqrtA = sqrt(A);

	double B = uksq*ujsq*K*K + 2*ujuk*K + 1.0;
	double sqrtB = sqrt(B);

	double num = 4*ujuk*K + tempk*tempj;
	double den = 2*sqrtA*sqrtB;

	double pd_num_yk = 4*yj*K + 2*yk*K*tempj;
	double pd_den_yk = -2.0*dy/sqrtA*sqrtB + 
		sqrtA*(2*yk*ujsq*K*K + 2*yj*K)/sqrtB;

	return (pd_num_yk*den - num*pd_den_yk)/(den*den);
}

/** Dot (inner) product of u and v */
double dot(double u[2], double v[2])
{
	return u[0]*v[0] + u[1]*v[1];
}

void eigenvalues(gsl_matrix *A, int inPlace) {

/*
  inPlace = 1 => A is replaced
  inPlace = 0 => A is retained
*/

   gsl_matrix *tmpA;

   if (inPlace)
      tmpA = A;
   else {
     tmpA = gsl_matrix_alloc(A->size1, A->size2);
     gsl_matrix_memcpy(tmpA , A);
   }

  int n = tmpA->size1;

  gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);

  gsl_eigen_nonsymmv_workspace * w =
    gsl_eigen_nonsymmv_alloc (n);

  gsl_eigen_nonsymmv (tmpA, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec,
                           GSL_EIGEN_SORT_ABS_DESC);

  {
    int i, j;

    for (i = 0; i < n; i++)
      {
        gsl_complex eval_i
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i
           = gsl_matrix_complex_column (evec, i);

        fprintf (stderr, "eigenvalue = %g + %gi\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
/*
        printf ("eigenvector = \n");
        for (j = 0; j < n; ++j)
          {
            gsl_complex z =
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
          }
*/
      }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  if (! inPlace)
    gsl_matrix_free(tmpA);
}

// Not used anymore
double det_get(gsl_matrix *A, int inPlace) {

/*
  inPlace = 1 => A is replaced with the LU decomposed copy.
  inPlace = 0 => A is retained, and a copy is used for LU.
*/

   printf("A->size1 = %ld, A->size2 = %ld\n", A->size1, A->size2);

   double det;
   int signum;
   gsl_permutation *p = gsl_permutation_alloc(A->size1);
   gsl_matrix *tmpA;

   if (inPlace)
      tmpA = A;
   else {
     tmpA = gsl_matrix_alloc(A->size1, A->size2);
     gsl_matrix_memcpy(tmpA , A);
   }


   gsl_linalg_LU_decomp(tmpA , p , &signum);
   det = gsl_linalg_LU_det(tmpA , signum);
   gsl_permutation_free(p);
   if (! inPlace)
      gsl_matrix_free(tmpA);

   
   return det;
}
