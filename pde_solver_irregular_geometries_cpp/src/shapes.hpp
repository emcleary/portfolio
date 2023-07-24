#ifndef SHAPE_H
#define SHAPE_H

#include <boost/math/constants/constants.hpp>
#include <iostream>
#include "utilities.hpp"

/// @brief Shape object defining shapes through parametric functions.
///        See children classes in the examples directory for
///        example implementations of parametric functions and boundary
///        conditions
class Shape {
public:
    /// @brief Constructor for Shape with a default ds value
    /// @param smin Minimum value of parametric s
    /// @param smax Maximum value of parametric s
    Shape(double smin, double smax) : m_smin(smin), m_smax(smax) {
        m_ds = (m_smax - m_smin) / 2;
    };

    /// @brief Constructor for Shape
    /// @param smin Minimum value of parametric s
    /// @param smax Maximum value of parametric s
    /// @param ds Initial step for parametric s when traversing the shape
    Shape(double smin, double smax, double ds) : m_smin(smin), m_smax(smax), m_ds(ds) {};

    virtual ~Shape(){};

    /// @brief Parametric function for x-coordinates
    /// @param s Parametric value
    /// @return Value of x at s
    virtual double f_x(double s) const = 0;

    /// @brief Parametric function for y-coordinates
    /// @param s Parametric value
    /// @return Value of y at s
    virtual double f_y(double s) const = 0;

    /// @brief Function used as boundary conditions on the shape surface in the PDE solver
    /// @param x x-coordinate
    /// @param y y-coordinate
    /// @return Value at the point
    virtual double boundary_conditions(double x, double y) const = 0;

    /// @brief Get minimum value of parameter s
    /// @return Value of smin
    double get_smin() { return m_smin; };

    /// @brief Get maximum value of parameter s
    /// @return Value of smax
    double get_smax() { return m_smax; };

    /// @brief Get step for parameter s
    /// @return Value of ds
    double get_ds() { return m_ds; };

private:
    /// @brief Minimum value of parameter s
    double m_smin;
    /// @brief Maximum value of parameter s
    double m_smax;
    /// @brief Step used for parameter s
    double m_ds;
};


/// @brief A point (x, y, s) for tracking input and outputs of parametric functions
struct point {
    double x;
    double y;
    double s;
    point(double xi, double yi, double si) : x(xi), y(yi), s(si) {};
};

#endif
