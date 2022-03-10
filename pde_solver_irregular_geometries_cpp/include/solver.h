#pragma once

#include <string>
#include <map>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/math/constants/constants.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

typedef boost::multi_array<double, 2> array_double_2d;
typedef boost::multi_array<double, 1> array_double_1d;
typedef boost::multi_array<bool, 2> array_bool_2d;

const double PI = boost::math::constants::pi<double>();

typedef struct {
    double smin = 0.0d;
    double smax = 2*PI;
    double xmin = -4.0d;
    double xmax = 4.0d;
    double ymin = -4.0d;
    double ymax = 4.0d;
    static double f_x(double s) { return (-9.0*std::sin(2*s) - 5.0*std::sin(3*s)) / 8.0; };
    static double f_y(double s) { return ( 9.0*std::cos(2*s) - 5.0*std::cos(3*s)) / 8.0; };
    static double boundary_conditions(double x, double y) { return std::log10(x*x+y*y); };
} Star;

typedef struct {
    double smin = 0.0d;
    double smax = 4.0d;
    double xmin = 0.0d;
    double xmax = 1.0d;
    double ymin = 0.0d;
    double ymax = 1.0d;
    static double f_x(double s) {
	if (s < 1)
	    return 0.3;
	else if (s < 2)
	    return 0.3 + (s-1) * 0.4;
	else if (s < 3)
	    return 0.7;
	else
	    return 0.7 - (s-3) * 0.4;
    };
    static double f_y(double s) {
	if (s < 1)
	    return 0.3 + s * 0.4;
	else if (s < 2)
	    return 0.7;
	else if (s < 3)
	    return 0.7 - (s - 2) * 0.4;
	else
	    return 0.3;
    };
    static double boundary_conditions(double x, double y) {
	return std::sin(PI * (x-0.5) / 0.4) * std::sin(PI * (y-0.5) / 0.4);
    };
} Square;


typedef struct {
    double smin = 0.0d;
    double smax = 6.28d;
    double xmin = 0.0d;
    double xmax = 1.0d;
    double ymin = 0.0d;
    double ymax = 1.0d;
    static double f_x(double s) { return 0.5 + std::cos(s) / 5.0; };
    static double f_y(double s) { return 0.5 + std::sin(s) / 5.0; };
    static double boundary_conditions(double x, double y) { return std::sin(4.0d * PI * (x - 0.5d)); };
} Circle;


struct rparams
{
    double a;
};

template <typename T>
class Solver {
private:
    int m_n;
    double m_dx, m_dy;
    double m_xmin, m_xmax, m_ymin, m_ymax;
    array_double_1d m_x, m_y;

    array_double_2d m_un, m_unp1;
    array_double_2d m_rho_l, m_rho_r, m_rho_u, m_rho_d;
    array_double_2d m_un_rho_l, m_un_rho_r, m_un_rho_u, m_un_rho_d;
    array_double_2d m_s0, m_s1, m_s2, m_s3, m_s4;
    array_bool_2d m_on_shape, m_in_shape;

    std::map<double, std::vector<double>> m_surf_x, m_surf_y;

    double m_smin, m_smax, m_ds;
    int m_iterations;
    double m_t, m_dt;

    double m_r;
    double m_atol;

    void calc_surf_points();
    void set_bools_in_shape();
    void set_bools_on_shape();
    void calc_dist_from_surf();
    void calc_scheme_parameters();

    // double boundary_conditions(double x, double y);
    double initial_conditions(double x, double y);

public:
    Solver(T shape, int n, double dt, int ns, double atol);

    ~Solver();
    
    void iteration(int n);
    void dump_solution(std::string filename);

    // methods setting up for and running gsl multiroot
    double root_finder(gsl_multiroot_function func, rparams p, double x_init[], size_t n);
    double root_finder_x(double x, double y, double s);
    double root_finder_y(double x, double y, double s);
    static int fsolver_function_x(const gsl_vector* x, void* params, gsl_vector* f);
    static int fsolver_function_y(const gsl_vector* y, void* params, gsl_vector* f);   
};

struct point {
    double x;
    double y;
    double s;
    point(double xi, double yi, double si) : x(xi), y(yi), s(si) {};
};

template class Solver<Circle>;
template class Solver<Star>;
template class Solver<Square>;
