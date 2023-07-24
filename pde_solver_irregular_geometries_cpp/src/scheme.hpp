#ifndef SOLVER_H
#define SOLVER_H

#include "grid.hpp"
#include "utilities.hpp"

/// @brief A first order accurate finite difference scheme for solving 2-d 
///        Poissons's equation around a shape with an irregular geometry using 
///        blocking out techniques adapted from "Numerical Partial 
///        Differential Equations, Part 2", Chapter 11.1 by J. B. Thomas.
class Scheme {
public:
    /// @brief A numerical schemes used for the PDE solver
    /// @param dt Timestep
    Scheme(double dt) : m_dt(dt) {}; 

    /// @brief Initializes solutions and parameters based off the shape
    void initialize();

    /// @brief Set the grid used for the solver
    /// @param grid 2-d regular grid with an irregular shape blocked out
    void set_grid(Grid *grid);

    /// @brief Solve the PDE for n iterations
    /// @param n Number of iterations
    void iteration(int n);

    /// @brief Get solutions of the PDE
    /// @return Solution at iteration n
    array_double_2d get_solution() { return m_un; };

private:

    /// @brief Calculate distances from grid points to the irregular shape's surace.
    void calc_dist_from_surf();

    /// @brief Calculate parameters for the numerical scheme
    void calc_scheme_param();

    /// @brief Solution at iteration n
    array_double_2d m_un;

    /// @brief Solution at iteration n+1
    array_double_2d m_unp1;

    /// @brief Distance from the x-dimension to the shape surface from the left
    array_double_2d m_rho_l;
    /// @brief Distance from the x-dimension to the shape surface from the right
    array_double_2d m_rho_r;
    /// @brief Distance from the y-dimension to the shape surface from above
    array_double_2d m_rho_u;
    /// @brief Distance from the y-dimension to the shape surface from below 
    array_double_2d m_rho_d;

    /// @brief Boundary condition for grid points to the left of the shape surface
    array_double_2d m_un_rho_l;
    /// @brief Boundary condition for grid points to the right of the shape surface
    array_double_2d m_un_rho_r;
    /// @brief Boundary condition for grid points from above the shape surface
    array_double_2d m_un_rho_u;
    /// @brief Boundary condition for grid points from below the shape surface
    array_double_2d m_un_rho_d;

    /// @brief Weight for the numerical scheme
    array_double_2d m_s0, m_s1, m_s2, m_s3, m_s4;

    /// @brief Regular grid with irregular shape used for the solver
    Grid *m_grid;

    /// @brief Timestep
    double m_dt;

    /// @brief Number of iterations
    int m_iterations{};

};

#endif
