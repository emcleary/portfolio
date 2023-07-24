#include "rootfinder.hpp"


/// @brief Additional arguments for the rootfinder
struct arguments {
    double value_fixed;
    Shape *shape;
};


/// @brief Itertaive multiroot solver using the gls library
/// @param f Struct including the function used for finding roots
/// @param x_init Array of parameters (x or y) and s
/// @param n Number of parameters
/// @return Returns first root (x or y)
double _root_finder(gsl_multiroot_function f, double x_init[], size_t n);


/// @brief Defines residuals for calculating roots along the x-axis
/// @param[in] parameters Tunable parameters: dimension y, parameter s
/// @param[in] args Additional arguments, including fixed x and parametric function
/// @param[out] residuals Differences between x, y and evaluations of the parametric function
int _fsolverfunction_x(const gsl_vector *parameters, void *args, gsl_vector *residuals);


/// @brief Defines residuals for calculating roots along the y-axis
/// @param[in] parameters Tunable parameters: dimension x, parameter s
/// @param[in] args Additional arguments, including fixed y and parametric function
/// @param[out] residuals Differences between x, y and evaluations of the parametric function
int _fsolverfunction_y(const gsl_vector *parameters, void *args, gsl_vector *residuals);


double _root_finder(gsl_multiroot_function f, double x_init[], size_t n) {
    gsl_vector *x_gsl = gsl_vector_alloc(n);
    gsl_vector_set(x_gsl, 0, x_init[0]);
    gsl_vector_set(x_gsl, 1, x_init[1]);

    const gsl_multiroot_fsolver_type *fs_type = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *solver_gsl = gsl_multiroot_fsolver_alloc(fs_type, 2);
    gsl_multiroot_fsolver_set(solver_gsl, &f, x_gsl);

    int status = GSL_CONTINUE;
    size_t iter = 0;
    size_t iter_max = 10000;
    while (status == GSL_CONTINUE && iter < iter_max) {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver_gsl);
        if (status)
            return std::nan("");

        status = gsl_multiroot_test_residual(solver_gsl->f, 1e-7);
    }

    if (iter == iter_max)
        return std::nan("");

    // We will only care about the first root (x or y)
    // We will never care about the second root s
    double root = gsl_vector_get(solver_gsl->x, 0);
    gsl_multiroot_fsolver_free(solver_gsl);
    gsl_vector_free(x_gsl);
    return root;
};

// Given x and initial y and s, find the best y
double root_finder_x(Shape *shape, double x, double y, double s) {
    struct arguments args = {x, shape};
    const size_t n = 2;
    gsl_multiroot_function f = {&_fsolverfunction_x, n, &args};
    double x_init[2] = {y, s};
    return _root_finder(f, x_init, n);
};


// Given y and initial x and s, find the best x
double root_finder_y(Shape *shape, double x, double y, double s) {
    struct arguments args = {y, shape};
    const size_t n = 2;
    gsl_multiroot_function f = {&_fsolverfunction_y, n, &args};
    double x_init[2] = {x, s};
    return _root_finder(f, x_init, n);
};


int _fsolverfunction_x(const gsl_vector *parameters, void *args, gsl_vector *residuals) {
    double x0 = ((struct arguments*) args)->value_fixed;
    Shape *shape = ((struct arguments*) args)->shape;

    double y = gsl_vector_get(parameters, 0);
    double s = gsl_vector_get(parameters, 1);

    const double f1 = y - shape->f_y(s);
    const double f2 = x0 - shape->f_x(s);

    gsl_vector_set(residuals, 0, f1);
    gsl_vector_set(residuals, 1, f2);

    return GSL_SUCCESS;
};


int _fsolverfunction_y(const gsl_vector *parameters, void *args, gsl_vector *residuals) {
    double y0 = ((struct arguments*) args)->value_fixed;
    Shape *shape = ((struct arguments*) args)->shape;

    double x = gsl_vector_get(parameters, 0);
    double s = gsl_vector_get(parameters, 1);

    const double f1 = y0 - shape->f_y(s);
    const double f2 = x - shape->f_x(s);

    gsl_vector_set(residuals, 0, f1);
    gsl_vector_set(residuals, 1, f2);

    return GSL_SUCCESS;
};
