#ifndef MANAGER_H
#define MANAGER_H

#include "grid.hpp"
#include "scheme.hpp"
#include "shapes.hpp"

/// @brief Class for managing the numerical solver following the builder design patter
class Solver {
public:

    /// @brief A manager from numerically solving the PDE
    /// @param dt Timestep
    Solver(double dt) {
        m_scheme = new Scheme(dt);
    }

    ~Solver() {
        delete m_scheme;
    }

    /// @brief Sets a new shape
    /// @param shape Irregular shape blocked out in the grid
    void set_shape(Shape *shape) {
        if (!m_new_shape)
            m_new_shape = m_shape != shape;
        m_shape = shape;
    }


    /// @brief Sets a new grid
    /// @param grid Regular grid containing dimensions and resolution of the domain
    void set_grid(Grid *grid) {
        m_new_grid = m_grid != grid;
        m_grid = grid;
        if (m_grid->get_shape())
            set_shape(m_grid->get_shape());
    }


    /// @brief Solver the PDE, initializing if a new shape and/or grid is set
    /// @param num_iter Number of iterations to solve
    void solve(int num_iter) {
        if (!m_shape) {
            std::cout << "Must set the shape before solving!\n";
        }

        if (!m_grid) {
            std::cout << "Must set the grid before solving!\n";
        }

        if (m_new_shape) {
            m_grid->set_shape(m_shape);
        }

        if (m_new_grid) {
            m_scheme->set_grid(m_grid);
        }

        // must initialize/update all coefficients for the scheme
        if (m_new_shape || m_new_grid) {
            m_scheme->initialize();
        }

        m_scheme->iteration(num_iter);
        reset_bools();
    }

    /// @brief Write solution to a file
    /// @param filename Output filename
    void dump_solution(std::string filename) {
        FILE *myfile;
        const array_double_2d soln = m_scheme->get_solution();
        myfile = fopen(filename.c_str(), "w");
        size_t px = soln.shape()[0]; // number of points in x
        size_t py = soln.shape()[1];
        fprintf(myfile, "grid %lu %lu\n", px, py);
        for (size_t j=0; j < py; j++)
            for (size_t i=0; i < px; i++) {
                if (m_grid->is_on_shape(i, j) || m_grid->is_in_shape(i, j))
                    fprintf(myfile, "%.6e %.6e NAN\n", m_grid->get_x(i), m_grid->get_y(j));
                else
                    fprintf(myfile, "%.6e %.6e %.6e\n", m_grid->get_x(i), m_grid->get_y(j), soln[i][j]);
            }
        fclose(myfile);
    }

private:

    /// @brief Resets members that track setters for initialization purposes
    void reset_bools() {
        m_new_shape = false;
        m_new_grid = false;
    }

    /// @brief Tracks if a new grid was set
    bool m_new_grid = false;

    /// @brief Tracks if a new shape was set
    bool m_new_shape = false;

    /// @brief Shape blocked out in the grid
    Shape *m_shape = 0;

    /// @brief Grid used for the PDE solver
    Grid *m_grid = 0;

    /// @brief Numerical scheme and PDE used
    Scheme *m_scheme = 0;

};

#endif
