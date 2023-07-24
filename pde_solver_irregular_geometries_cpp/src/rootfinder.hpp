#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "shapes.hpp"


/// @brief Given x and initial y and s, find the best y on the shape's surface
/// @param shape Shape pointer passed in for its parametric functions
/// @param x Fixed value of the x-dimension
/// @param y Initial value of the y-dimension
/// @param s Initial value of the s parameter
/// @return Root value y on the shape surface at x
double root_finder_x(Shape *shape, double x, double y, double s);


/// @brief Given y and initial x and s, find the best x on the shape's surface
/// @param shape Shape pointer passed in for its parametric functions
/// @param x Initial value of the x-dimension
/// @param y Fixed value of the y-dimension
/// @param s Initial value of the s parameter
/// @return Root value x on the shape surface at y
double root_finder_y(Shape *shape, double x, double y, double s);


#endif
