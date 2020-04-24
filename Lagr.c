/*! \file
    \brief Lagrangian of the Curved Three Body Problem (C3BP).

    $Author: roldan $
    $Date: $
*/

#include <math.h>       // M_PI
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // memcpy
#include <assert.h>
#include <gsl/gsl_errno.h>  // GSL_SUCCESS

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

/** Lagrangian of the C3BP */
double Lagr(double z[6], double kappa, double m[3]);

/** Kinetic energy */
double kinetic(double z[6], double kappa, double m[3]);

/** Conformal factor */
double conf(double z[2], double kappa);

/** Norm square */
double normsq(double z[2]);

/** Potential energy */
double potential(double z[6], double kappa, double m[3]);
double V(double *zj, double *zk, double kappa);

/** Norm */
double norm(double z[2]);

/** Gradient of Lagrangian L w.r.t. positions u */
void gradL(double u[6], double K, double m[3], double grad[6]);

/** Partial derivative of T (kinetic energy) w.r.t. x_k */
double pd_T_xk(double u[6], int k, double K, double m[3]);
double pd_T_yk(double u[6], int k, double K, double m[3]);

/** Partial derivative of U (potential energy) w.r.t. x_k */
double pd_U_xk(double u[6], int k, double K, double m[3]);
double pd_U_yk(double u[6], int k, double K, double m[3]);

double pd_V_xk(double xj, double yj, double xk, double yk, double K);
double pd_V_yk(double xj, double yj, double xk, double yk, double K);

/** Dot (inner) product of u and v */
double dot(double u[2], double v[2]);

struct grad_params
{
	double K;		// curvature \kappa
	double m[3];	// masses m_1, m_2, m_3
	double a[6];	// RE
};

int grad_f(const gsl_vector *x, void *params, gsl_vector *f)
{
	double K = ((struct grad_params *)params)->K;
	double m[3];
	double a[6];

	for(int i=0; i<3; i++)
		m[i] = ((struct grad_params *)params)->m[i];

	for(int i=0; i<6; i++)
		a[i] = ((struct grad_params *)params)->a[i];

	double u[6];

	for(int i=0; i<6; i++)
	{
		u[i] = gsl_vector_get(x, i);
	}

	double alpha = gsl_vector_get(x,6);	// Lagrange multiplier

	double grad[6];
	gradL(u, K, m, grad);

	gsl_vector_set(f, 0, grad[0] - alpha*u[1]);
	gsl_vector_set(f, 1, grad[1] + alpha*u[0]);
	gsl_vector_set(f, 2, grad[2] - alpha*u[3]);
	gsl_vector_set(f, 3, grad[3] + alpha*u[2]);
	gsl_vector_set(f, 4, grad[4] - alpha*u[5]);
	gsl_vector_set(f, 5, grad[5] + alpha*u[4]);
	gsl_vector_set(f, 6, 
			-u[0]*a[1] + u[1]*a[0] \
			-u[2]*a[3] + u[3]*a[2] \
			-u[4]*a[5] + u[5]*a[4]);
	return GSL_SUCCESS;
}

int grad_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double K = ((struct grad_params *)params)->K;
	double m[3];
	double a[6];

	int n = 7;

	for(int i=0; i<3; i++)
		m[i] = ((struct grad_params *)params)->m[i];

	for(int i=0; i<6; i++)
		a[i] = ((struct grad_params *)params)->a[i];

	gsl_multiroot_function f = {&grad_f, n, params};

	gsl_vector *f_vec = gsl_vector_alloc(n);
	grad_f(x,params,f_vec);
	gsl_multiroot_fdjacobian (&f, x, f_vec, GSL_SQRT_DBL_EPSILON, J);

	gsl_matrix_set(J, 0, 6, -gsl_vector_get(x, 1));
	gsl_matrix_set(J, 1, 6,  gsl_vector_get(x, 0));
	gsl_matrix_set(J, 2, 6, -gsl_vector_get(x, 3));
	gsl_matrix_set(J, 3, 6,  gsl_vector_get(x, 2));
	gsl_matrix_set(J, 4, 6, -gsl_vector_get(x, 5));
	gsl_matrix_set(J, 5, 6,  gsl_vector_get(x, 4));
	gsl_matrix_set(J, 6, 6,  0.0);

	/*
	fprintf(stderr, "JACOBIAN:\n");
	for(int i=0; i<n; i++)
	{
		for(int k=0; k<n; k++)
			fprintf(stderr, "% .14e ", gsl_matrix_get(J,i,k));
		fprintf(stderr, "\n");
	}
	*/

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
  fprintf (stderr, "iter = %3lu\n"
		  "x = % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n"
          "f(x) = % .5e % .5e % .5e % .5e % .5e % .5e % .5e\n\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4),
          gsl_vector_get (s->x, 5),
          gsl_vector_get (s->x, 6),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3),
          gsl_vector_get (s->f, 4),
          gsl_vector_get (s->f, 5),
          gsl_vector_get (s->f, 6));
}

int main( )
{
	// Lagrange's equilateral RE, equal masses m1=m2=m3=sqrt(3)
	//double m[3] = {sqrt(3), sqrt(3), sqrt(3)};
	//double alpha = 0.5;

	// Lagrange's equilateral RE, equal masses m1=m2=m3=1
	//double m[3] = {1, 1, 1};
	//double alpha = 0.5*cbrt(1.0/sqrt(3.0));

	// Lagrange's equilateral RE, different masses m1=1, m2=2, m3=3
	double m[3] = {1, 2, 3};
	//double alpha = 2.0;
	//double alpha = 0.5*cbrt(2.0)*cbrt(1.0/sqrt(3.0));

	// positions (x_j, y_j) for each of the 3 bodies.

	// Lagrange's equilateral RE, equal masses
	//double a[6] = {alpha, 0, -alpha/2, sqrt(3)*alpha/2, -alpha/2,
	//	-sqrt(3)*alpha/2};

	// Lagrange's equilateral RE, different masses m1=1, m2=2, m3=3
	double tau = cbrt(6)/2.0;		// length of side of eq triangle
	double a[6] = {0.655696914638530, 0.0757133580346725, 
		-0.131139382927706, 0.529993506242707, 
		-0.131139382927706, -0.378566790173362};

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

	int num_K = 10000;
	double delta_K = -100.0/num_K;

	const size_t n = 7;
	int i;
	double x_init[n];
	for(i=0; i<6; i++)
		x_init[i] = a[i];

	x_init[6] = 0;	// alpha (Lagrange multiplier)

	for(int j=0; j<num_K; j++)
	{
		kappa = j*delta_K;

		const gsl_multiroot_fdfsolver_type *T;
		gsl_multiroot_fdfsolver *s;	

		int status;
		size_t iter = 0;

		struct grad_params p; 
		p.K = kappa;
		for(int i=0; i<3; i++)
			p.m[i] = m[i];
		for(int i=0; i<6; i++)
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
				gsl_multiroot_test_residual (s->f, 1e-13);
			}
		while (status == GSL_CONTINUE && iter < 1000);

		fprintf (stderr, "status = %s\n", gsl_strerror (status));

		if (status)   /* check if solver is stuck */
		{
			print_state (iter, s);
			exit(EXIT_FAILURE);
		}

		// Output relative equilibrium
		printf ("% .2f % .14f % .14f % .14f % .14f % .14f % .14f % .14f\n",
				kappa,
				gsl_vector_get (s->x, 0),
				gsl_vector_get (s->x, 1),
				gsl_vector_get (s->x, 2),
				gsl_vector_get (s->x, 3),
				gsl_vector_get (s->x, 4),
				gsl_vector_get (s->x, 5),
				gsl_vector_get (s->x, 6));

		// Use current root as initial guess for next curvature
		for(i=0; i<n; i++)
			x_init[i] = gsl_vector_get (s->x, i);

		gsl_multiroot_fdfsolver_free (s);
		gsl_vector_free (x);
	}

	exit(EXIT_SUCCESS);
}

/** 
  Lagrangian of the C3BP in rotating coordinates.

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

double Lagr(double z[6], double kappa, double m[3])
{
	return kinetic(z,kappa,m) + potential(z,kappa,m);
}

/** Kinetic energy */
double kinetic(double z[6], double kappa, double m[3])
{
	int j;

	double T = 0;
	for(j=0; j<3; j++)
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
double potential(double z[6], double kappa, double m[3])
{
	int j, k;

	double U = 0;
	for(k=0; k<3; k++)
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
  positions (x_j, y_j) for each of the 3 bodies.

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the 3 bodies

  \param[out] grad[6]
  Value of the gradient evaluated at u.

  \returns 
  */

void gradL(double u[6], double K, double m[3], double grad[6])
{
	for(int i=0; i<3; i++)
	{
		grad[2*i  ] = pd_T_xk(u,i,K,m) + pd_U_xk(u,i,K,m);	// w.r.t. x_k
		grad[2*i+1] = pd_T_yk(u,i,K,m) + pd_U_yk(u,i,K,m);	// w.r.t. y_k
	}
}

/** 
  Partial derivative of T (kinetic energy) w.r.t. x_k

  \param[in] u
  positions (x_j, y_j) for each of the 3 bodies.

  \param[in] k
  differentiate w.r.t. x_k (0 <= k <= 2)

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the 3 bodies

  \returns Value of the partial derivative evaluated at u.
  */

double pd_T_xk(double u[6], int k, double K, double m[3])
{
	assert(k >= 0 && k <= 2);
	double *uk = u + 2*k;
	double mk = m[k];
	double xk = uk[0];

	return mk*conf(uk, K)*((-2.0*K*xk)/(1+K*normsq(uk))*normsq(uk) + xk);
}

double pd_T_yk(double u[6], int k, double K, double m[3])
{
	assert(k >= 0 && k <= 2);
	double *uk = u + 2*k;
	double mk = m[k];
	double yk = uk[1];

	return mk*(conf(uk, K)*(-2.0*K*yk)/(1+K*normsq(uk))*normsq(uk) + 
				conf(uk, K)*yk);
}

/** 
  Partial derivative of U (potential energy) w.r.t. x_k

  \param[in] u
  positions (x_j, y_j) for each of the 3 bodies.

  \param[in] k
  differentiate w.r.t. x_k (0 <= k <= 2)

  \param[in] K
  curvature kappa.

  \param[in] m
  masses m_j for each of the 3 bodies

  \returns Value of the partial derivative evaluated at u.
  */

double pd_U_xk(double u[6], int k, double K, double m[3])
{
	double *uk = u + 2*k;
	double xk = uk[0];
	double yk = uk[1];

	double sum = 0;
	for(int j=0; j<3; j++)
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

double pd_U_yk(double u[6], int k, double K, double m[3])
{
	double *uk = u + 2*k;
	double xk = uk[0];
	double yk = uk[1];

	double sum = 0;
	for(int j=0; j<3; j++)
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

