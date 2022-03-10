#pragma once

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "../include/solver.h"
#include <boost/multi_array.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <boost/math/constants/constants.hpp>
#include <omp.h>

template <typename T>
Solver<T>::Solver(T shape, int n, double dt, int ns, double atol) {
    m_atol = atol;
    m_t = 0;
    m_iterations = 0;
    m_dt = dt;

    m_smin = shape.smin;
    m_smax = shape.smax;
    m_ds = (m_smax - m_smin) / (ns + 1);

    m_n = n;
    m_xmin = shape.xmin;
    m_xmax = shape.xmax;
    m_ymin = shape.ymin;
    m_ymax = shape.ymax;
    m_dx = (m_xmax - m_xmin) / m_n;
    m_dy = (m_ymax - m_ymin) / m_n;
    m_r = m_dt / m_dx / m_dx;

    // Allocate and fill axes
    m_x.resize(boost::extents[m_n+1]);
    m_y.resize(boost::extents[m_n+1]);
    for (int i=0; i < m_n+1; i++)
    {
	m_x[i] = m_xmin + i * m_dx;
	m_y[i] = m_ymin + i * m_dy;
    }

    // Arrays for the scheme
    m_rho_l.resize(boost::extents[m_n+1][m_n+1]);
    m_rho_r.resize(boost::extents[m_n+1][m_n+1]);
    m_rho_u.resize(boost::extents[m_n+1][m_n+1]);
    m_rho_d.resize(boost::extents[m_n+1][m_n+1]);

    m_un_rho_l.resize(boost::extents[m_n+1][m_n+1]);
    m_un_rho_r.resize(boost::extents[m_n+1][m_n+1]);
    m_un_rho_u.resize(boost::extents[m_n+1][m_n+1]);
    m_un_rho_d.resize(boost::extents[m_n+1][m_n+1]);

    m_on_shape.resize(boost::extents[m_n+1][m_n+1]);
    m_in_shape.resize(boost::extents[m_n+1][m_n+1]);

    m_s0.resize(boost::extents[m_n+1][m_n+1]);
    m_s1.resize(boost::extents[m_n+1][m_n+1]);
    m_s2.resize(boost::extents[m_n+1][m_n+1]);
    m_s3.resize(boost::extents[m_n+1][m_n+1]);
    m_s4.resize(boost::extents[m_n+1][m_n+1]);

    // State variables
    m_un.resize(boost::extents[m_n+1][m_n+1]);
    m_unp1.resize(boost::extents[m_n+1][m_n+1]);

    // Initial values
    for (int j=0; j < m_n+1; j++)
	for (int i=0; i < m_n+1; i++)
	{
	    m_rho_l[i][j] = 1.0d;
	    m_rho_r[i][j] = 1.0d;
	    m_rho_u[i][j] = 1.0d;
	    m_rho_d[i][j] = 1.0d;

	    m_un_rho_l[i][j] = 1.0d;
	    m_un_rho_r[i][j] = 1.0d;
	    m_un_rho_u[i][j] = 1.0d;
	    m_un_rho_d[i][j] = 1.0d;

	    m_s0[i][j] = 0.0d;
	    m_s1[i][j] = 0.0d;
	    m_s2[i][j] = 0.0d;
	    m_s3[i][j] = 0.0d;
	    m_s4[i][j] = 0.0d;

	    m_un[i][j] = 0.0d;
	    m_unp1[i][j] = 0.0d;

	    m_in_shape[i][j] = false;
	    m_on_shape[i][j] = false;
	}

    // Initialize surfs
    for (const auto & xi : m_x)
    {
	m_surf_x[xi] = std::vector<double>();
    }
    for (const auto & yi : m_y)
    {
	m_surf_y[yi] = std::vector<double>();
    }

    // Calculate and set all scheme parameters needed
    calc_surf_points();
    set_bools_in_shape();
    set_bools_on_shape();
    calc_dist_from_surf();
    calc_scheme_parameters();
};

template <typename T>
Solver<T>::~Solver()
{

}

template <typename T>
void Solver<T>::calc_surf_points()
{
    double sn = m_smin;
    double xn = T::f_x(sn);
    double yn = T::f_y(sn);
    std::vector<point*> points;
    point *pt = new point(xn, yn, sn);
    points.push_back(pt);

    while (sn < m_smax)
    {
	double x_ref = points.back()->x;
	double y_ref = points.back()->y;
	double s_ref = points.back()->s;
	
	double t_ds = m_ds;
	xn = T::f_x(sn + t_ds);
	yn = T::f_y(sn + t_ds);
	while (std::abs(xn - x_ref) > (m_dx - m_atol) || (std::abs(yn - y_ref) > (m_dy - m_atol)))
	{
	    t_ds = t_ds / 2;
	    xn = T::f_x(sn + t_ds);
	    yn = T::f_y(sn + t_ds);
	}
	sn = sn + t_ds;
	pt = new point(xn, yn, sn);
	points.push_back(pt);
    }

    for (int i=0; i < points.size()-1; i++)
    {
	point* pa = points[i];
	point* pb = points[i+1];
	double xa = pa->x;
	double ya = pa->y;
	double sa = pa->s;
	double xb = pb->x;
	double yb = pb->y;
	double sb = pb->s;

	if (std::abs(xa - xb) < m_atol)
	    continue;

	int idx = -1;
	for (int j=0; j < m_n+1; j++)
	{
	    bool b1 = ((xa + m_atol) > m_x[j]) & ((xb - m_atol) < m_x[j]);
	    bool b2 = ((xa - m_atol) < m_x[j]) & ((xb + m_atol) > m_x[j]);
	    if (b1 || b2)
		idx = j;
	}
	if (idx == -1)
	    continue;

	double xj = m_x[idx];
	double yj = root_finder_x(xj, ya, sa);
	if (!std::isnan(yj))
	{
	    m_surf_x[xj].push_back(yj);
	    continue;
	}

	yj = root_finder_x(xj, ya, sb);
	if (!std::isnan(yj))
	{
	    m_surf_x[xj].push_back(yj);
	    continue;
	}
	yj = root_finder_x(xj, yb, sa);
	if (!std::isnan(yj))
	{
	    m_surf_x[xj].push_back(yj);
	    continue;
	}
	yj = root_finder_x(xj, yb, sb);
	if (!std::isnan(yj))
	{
	    m_surf_x[xj].push_back(yj);
	    continue;
	}

	std::cout << "NO SOLUTION FOUND!" << std::endl;
	return;
    }

    for (int i=0; i < points.size()-1; i++)
    {
	point* pa = points[i];
	point* pb = points[i+1];
	double xa = pa->x;
	double ya = pa->y;
	double sa = pa->s;
	double xb = pb->x;
	double yb = pb->y;
	double sb = pb->s;

	if (std::abs(ya - yb) < m_atol)
	    continue;

	int idx = -1;
	for (int j=0; j < m_n+1; j++)
	{
	    bool b1 = ((ya + m_atol) > m_y[j]) & ((yb - m_atol) < m_y[j]);
	    bool b2 = ((ya - m_atol) < m_y[j]) & ((yb + m_atol) > m_y[j]);
	    if (b1 || b2)
		idx = j;
	}
	if (idx == -1)
	    continue;

	double yj = m_y[idx];
	double xj = root_finder_y(xa, yj, sa);
	if (!std::isnan(xj))
	{
	    m_surf_y[yj].push_back(xj);
	    continue;
	}

	xj = root_finder_y(xa, yj, sb);
	if (!std::isnan(xj))
	{
	    m_surf_y[yj].push_back(xj);
	    continue;
	}

	xj = root_finder_y(xb, yj, sa);
	if (!std::isnan(xj))
	{
	    m_surf_y[yj].push_back(xj);
	    continue;
	}

	xj = root_finder_y(xb, yj, sb);
	if (!std::isnan(xj))
	{
	    m_surf_y[yj].push_back(xj);
	    continue;
	}

	std::cout << "NO SOLUTION FOUND!" << std::endl;
	return;
    }
}

template <typename T>
int Solver<T>::fsolver_function_x(const gsl_vector* x, void* params, gsl_vector* f)
{
    double x0 = ((struct rparams*) params)->a;

    const double y = gsl_vector_get (x, 0);
    const double s = gsl_vector_get (x, 1);

    const double f1 = y - T::f_y(s);
    const double f2 = x0 - T::f_x(s);

    gsl_vector_set(f, 0, f1);
    gsl_vector_set(f, 1, f2);

    return GSL_SUCCESS;
}

template <typename T>
int Solver<T>::fsolver_function_y(const gsl_vector* y, void* params, gsl_vector* f)
{
    double y0 = ((struct rparams*) params)->a;

    const double x = gsl_vector_get (y, 0);
    const double s = gsl_vector_get (y, 1);

    const double f1 = y0 - T::f_y(s);
    const double f2 = x - T::f_x(s);

    gsl_vector_set(f, 0, f1);
    gsl_vector_set(f, 1, f2);

    return GSL_SUCCESS;
}

template <typename T>
double Solver<T>::root_finder(gsl_multiroot_function f, rparams p, double x_init[], size_t n)
{
    gsl_vector *x_gsl = gsl_vector_alloc (n);
    gsl_vector_set (x_gsl, 0, x_init[0]);
    gsl_vector_set (x_gsl, 1, x_init[1]);

    const gsl_multiroot_fsolver_type *fs_type = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *solver_gsl = gsl_multiroot_fsolver_alloc (fs_type, 2);
    gsl_multiroot_fsolver_set (solver_gsl, &f, x_gsl);

    int status;
    size_t iter = 0;
    do
    {
	iter++;
	status = gsl_multiroot_fsolver_iterate (solver_gsl);
	if (status)   /* check if solver is stuck */
	    return std::nan("");

	status = gsl_multiroot_test_residual (solver_gsl->f, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 1000);
    if (iter == 1000)
	return std::nan("");
    
    double root = gsl_vector_get (solver_gsl->x, 0);
    gsl_multiroot_fsolver_free (solver_gsl);
    gsl_vector_free (x_gsl);
    return root;
}

template <typename T>
double Solver<T>::root_finder_x(double x, double y, double s)
{
    struct rparams p = {x};
    const size_t n = 2; // number of variables to solve for; NOT size of rparam p
    gsl_multiroot_function f = {&fsolver_function_x, n, &p};
    double x_init[2] = {y, s};
    return root_finder(f, p, x_init, n);
}

template <typename T>
double Solver<T>::root_finder_y(double x, double y, double s)
{
    struct rparams p = {y};
    const size_t n = 2; // number of variables to solve for; NOT size of rparam p
    gsl_multiroot_function f = {&fsolver_function_y, n, &p};
    double x_init[2] = {x, s};
    return root_finder(f, p, x_init, n);
}

template <typename T>
void Solver<T>::set_bools_in_shape()
{
    bool out_x;
    bool out_y;
    double min_elem;
    double max_elem;
    double yj;
    double xi;
    for (int j=0; j < m_n+1; j++)
    {
	yj = m_y[j];
	for (int i=0; i < m_n+1; i++)
	{
	    xi = m_x[i];
	    out_x = true;
	    out_y = true;
	    if (m_surf_x[xi].size() > 0)
	    {
		min_elem = std::min_element(m_surf_x[xi].begin(), m_surf_x[xi].end())[0];
		max_elem = std::max_element(m_surf_x[xi].begin(), m_surf_x[xi].end())[0];
		out_x = (yj < min_elem + m_atol) || (yj > max_elem - m_atol);
	    }
	    if (m_surf_y[yj].size() > 0)
	    {
		min_elem = std::min_element(m_surf_y[yj].begin(), m_surf_y[yj].end())[0];
		max_elem = std::max_element(m_surf_y[yj].begin(), m_surf_y[yj].end())[0];
		out_y = (xi < min_elem + m_atol) || (xi > max_elem - m_atol);
	    }
	    m_in_shape[i][j] = !(out_x || out_y);
	}
    }
}

template <typename T>
void Solver<T>::set_bools_on_shape()
{
    double xi;
    bool b;
    int idx;
    for (int i=0; i < m_n+1; i++)
    {
	xi = m_x[i];
	for (const auto & yi : m_surf_x[xi])
	{
	    idx = -1;
	    for (int j=0; j < m_n+1; j++)
	    {
		b = (yi - m_atol < m_y[j]) & (yi + m_atol > m_y[j]);
		if (b)
		    idx = j;
	    }
	    if (idx == -1)
		continue;

	    m_on_shape[i][idx] = !m_in_shape[i][idx];
	}
    }

    double yi;
    for (int i=0; i < m_n+1; i++)
    {
	yi = m_y[i];
	for (const auto & xi : m_surf_y[yi])
	{
	    idx = -1;
	    for (int j=0; j < m_n+1; j++)
	    {
		b = (xi - m_atol < m_x[j]) & (xi + m_atol > m_x[j]);
		if (b)
		    idx = j;
	    }
	    if (idx == -1)
		continue;

	    m_on_shape[idx][i] = !m_in_shape[idx][i];
	}
    }

}

template <typename T>
void Solver<T>::calc_dist_from_surf()
{
    double distance;
    for (int i=0; i < m_n+1; i++)
    {
	double xi = m_x[i];
	for (const auto & yi : m_surf_x[xi])
	{
	    int j = -1;
	    for (int k=0; k < m_n+1; k++)
	    {
		if (m_y[k] > yi + m_atol)
		{
		    j = k - 1;
		    break;
		}
	    }
	    if (m_on_shape[i][j])
		continue;
	    if (!m_in_shape[i][j] && m_in_shape[i][j+1])
	    {
		distance = yi - m_y[j];
		m_rho_u[i][j] = distance / m_dy;
		m_un_rho_u[i][j] = T::boundary_conditions(xi, yi);
	    }
	    else if (!m_in_shape[i][j+1] && m_in_shape[i][j])
	    {
		distance = m_y[j+1] - yi;
		m_rho_d[i][j] = distance / m_dy;
		m_un_rho_d[i][j+1] = T::boundary_conditions(xi, yi);
	    }
	}
    }

    for (int i=0; i < m_n+1; i++)
    {
	double yi = m_y[i];
	for (const auto & xi : m_surf_y[yi])
	{
	    int j = -1;
	    for (int k=0; k < m_n+1; k++)
	    {
		if (m_x[k] > xi + m_atol)
		{
		    j = k - 1;
		    break;
		}
	    }
	    if (m_on_shape[j][i])
		continue;
	    if (!m_in_shape[j][i] && m_in_shape[j+1][i])
	    {
		distance = xi - m_x[j];
		m_rho_r[j][i] = distance / m_dx;
		m_un_rho_r[j][i] = T::boundary_conditions(xi, yi);
	    }
	    else if (!m_in_shape[j+1][i] && m_in_shape[j][i])
	    {
		distance = m_x[j+1] - xi;
		m_rho_l[j+1][i] = distance / m_dx;
		m_un_rho_l[j+1][i] = T::boundary_conditions(xi, yi);
	    }
	}
    }
}

template <typename T>
void Solver<T>::calc_scheme_parameters()
{
    for (int j=0; j < m_n+1; j++)
	for (int i=0; i < m_n+1; i++)
	{
	    m_s0[i][j] = 1.0d - 2.0d * m_r / m_rho_l[i][j] / m_rho_r[i][j] - 2 * m_r / m_rho_u[i][j] / m_rho_d[i][j];
	    m_s1[i][j] = 4 * m_r / (1 + m_rho_l[i][j]) / m_rho_r[i][j] / (1 + m_rho_r[i][j]);
	    m_s2[i][j] = 4 * m_r / (1 + m_rho_r[i][j]) / m_rho_l[i][j] / (1 + m_rho_l[i][j]);
	    m_s3[i][j] = 4 * m_r / (1 + m_rho_d[i][j]) / m_rho_u[i][j] / (1 + m_rho_u[i][j]);
	    m_s4[i][j] = 4 * m_r / (1 + m_rho_u[i][j]) / m_rho_d[i][j] / (1 + m_rho_d[i][j]);

	    m_un[i][j] = initial_conditions(m_x[i], m_y[j]);
	    if (m_on_shape[i][j])
		m_un[i][j] = T::boundary_conditions(m_x[i], m_y[j]);
	    else if (m_in_shape[i][j])
		m_un[i][j] = 1.0d;
	}
}

template <typename T>
void Solver<T>::iteration(int n)
{
    int i, j, k, m;
    for (m=0; m < n; m++)
    {
#pragma omp parallel for private(k,i,j)
	for (k=0; k < (m_n+1)*(m_n+1); k++)
	{
	    i = k / (m_n+1);
	    j = k % (m_n+1);

	    // BCs on shape surface
	    if (m_on_shape[i][j] || m_in_shape[i][j])
		m_unp1[i][j] = m_un[i][j];

	    // BCs on edges
	    else if (i == 0)
		m_unp1[i][j] = 0;
	    else if (i == m_n)
		m_unp1[i][j] = 0;
	    else if (j == 0)
		m_unp1[i][j] = 0;
	    else if (j == m_n)
		m_unp1[i][j] = 0;

	    // Numerical scheme
	    else
	    {
		m_unp1[i][j] = m_s0[i][j] * m_un[i][j]
		    + m_s1[i][j] * m_un[i+1][j] * m_un_rho_r[i][j]
		    + m_s2[i][j] * m_un[i-1][j] * m_un_rho_l[i][j]
		    + m_s3[i][j] * m_un[i][j+1] * m_un_rho_u[i][j]
		    + m_s4[i][j] * m_un[i][j-1] * m_un_rho_d[i][j];
	    }
	}

	// Swap pointers
	std::swap(m_un, m_unp1);

	m_iterations = m_iterations + 1;
	if (m_iterations % 1000 == 0)
	{
	    std::cout << "Iteration " << m_iterations << std::endl;
	}
    }
}

template <typename T>
double Solver<T>::initial_conditions(double x, double y)
{
    return 0;
}

template <typename T>
void Solver<T>::dump_solution(std::string filename)
{
    FILE *myfile;
    myfile = fopen(filename.c_str(), "w");
    fprintf(myfile, "grid %u %u\n", m_n+1, m_n+1);
    for (int j=0; j < m_n+1; j++)
	for (int i=0; i < m_n+1; i++)
	{
	    if (m_on_shape[i][j] || m_in_shape[i][j])
		fprintf(myfile, "%.6e %.6e NAN\n", m_x[i], m_y[j]);
	    else
		fprintf(myfile, "%.6e %.6e %.6e\n", m_x[i], m_y[j], m_un[i][j]);
	}
    fclose(myfile);
}
