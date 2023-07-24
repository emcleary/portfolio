#include "solver.hpp"
#include <iostream>


void Scheme::initialize() {
    m_iterations = 0;

    // initialize solutions
    for (size_t i = 0; i < m_un.shape()[0]; ++i) {
        for (size_t j = 0; j < m_un.shape()[1]; ++j) {
            m_un[i][j] = 0.0d;
            m_unp1[i][j] = 0.0d;
        }
    }

    // Update whatever needs to be modified when setting a new grid
    calc_dist_from_surf();
    calc_scheme_param();
}


void Scheme::set_grid(Grid *grid) {
    m_grid = grid;

    // allocate/resize all arrays using the grid as needed
    m_grid->allocate_2d_array(&m_un);
    m_grid->allocate_2d_array(&m_unp1);

    m_grid->allocate_2d_array(&m_rho_l);
    m_grid->allocate_2d_array(&m_rho_r);
    m_grid->allocate_2d_array(&m_rho_u);
    m_grid->allocate_2d_array(&m_rho_d);

    m_grid->allocate_2d_array(&m_un_rho_l);
    m_grid->allocate_2d_array(&m_un_rho_r);
    m_grid->allocate_2d_array(&m_un_rho_u);
    m_grid->allocate_2d_array(&m_un_rho_d);

    m_grid->allocate_2d_array(&m_s0);
    m_grid->allocate_2d_array(&m_s1);
    m_grid->allocate_2d_array(&m_s2);
    m_grid->allocate_2d_array(&m_s3);
    m_grid->allocate_2d_array(&m_s4);
}


void Scheme::calc_dist_from_surf() {
    std::cout << "Calculating distances from shape surface\n";

    // initialize all arrays set here
    for (size_t i = 0; i < m_rho_d.shape()[0]; ++i) {
        for (size_t j = 0; j < m_rho_d.shape()[1]; ++j) {
            m_rho_l[i][j] = 1.0d;
            m_rho_r[i][j] = 1.0d;
            m_rho_u[i][j] = 1.0d;
            m_rho_d[i][j] = 1.0d;

            m_un_rho_l[i][j] = 1.0d;
            m_un_rho_r[i][j] = 1.0d;
            m_un_rho_u[i][j] = 1.0d;
            m_un_rho_d[i][j] = 1.0d;
        }
    }

    double distance;
    double xi, yj;
    double dx = m_grid->get_dx();
    double dy = m_grid->get_dy();
    bool in_shape_1, in_shape_2;

    // calculate values above and below surface points
    for (size_t i = 0; i < m_rho_d.shape()[0]; ++i) {
        xi = m_grid->get_x(i);
        for (const auto &yi : m_grid->get_surf_x(i)) {
            int j = m_grid->get_max_index_y(yi);
            if (m_grid->is_on_shape(i, j))
                continue;

            in_shape_1 = m_grid->is_in_shape(i, j);
            in_shape_2 = m_grid->is_in_shape(i, j + 1);
            if (!in_shape_1 && in_shape_2) {
                yj = m_grid->get_y(j);
                distance = yi - yj;
                m_rho_u[i][j] = distance / dy;
                m_un_rho_u[i][j] = m_grid->get_surface_value(xi, yi);
            } else if (in_shape_1 && !in_shape_2) {
                yj = m_grid->get_y(j + 1);
                distance = yj - yi;
                m_rho_d[i][j + 1] = distance / dy;
                m_un_rho_d[i][j + 1] = m_grid->get_surface_value(xi, yi);
            }
        }
    }

    // calculate values left and right of surface points
    for (size_t j = 0; j < m_rho_d.shape()[1]; ++j) {
        yj = m_grid->get_y(j);
        for (const auto &xj : m_grid->get_surf_y(j)) {
            int i = m_grid->get_max_index_x(xj);
            if (m_grid->is_on_shape(i, j))
                continue;

            in_shape_1 = m_grid->is_in_shape(i, j);
            in_shape_2 = m_grid->is_in_shape(i + 1, j);
            if (!in_shape_1 && in_shape_2) {
                xi = m_grid->get_x(i);
                distance = xj - xi;
                m_rho_r[i][j] = distance / dx;
                m_un_rho_r[i][j] = m_grid->get_surface_value(xj, yj);
            } else if (in_shape_1 && !in_shape_2) {
                xi = m_grid->get_x(i + 1);
                distance = xi - xj;
                m_rho_l[i + 1][j] = distance / dx;
                m_un_rho_l[i + 1][j] = m_grid->get_surface_value(xj, yj);
            }
        }
    }
};


void Scheme::calc_scheme_param() {
    std::cout << "Calculating scheme parameters\n";

    double dx = m_grid->get_dx();
    double dy = m_grid->get_dy();

    double cfl_dt = 0.002 / std::max(1./dx/dx, 1./dy/dy);
    if (m_dt > cfl_dt) {
        std::cout << "Timestep might be too large for a stable solution.\n";
        std::cout << "Make sure to double check this and manually reduce the timestep if necessary.\n";
        std::cout << cfl_dt << " is a reasonable estimate of a stable timestep.\n";
    }

    double rx = m_dt / dx / dx;
    double ry = m_dt / dy / dy;

    for (size_t i = 0; i < m_s0.shape()[0]; ++i) {
        for (size_t j = 0; j < m_s0.shape()[1]; ++j) {

            m_s0[i][j] = 1.0d - 2 * rx / m_rho_l[i][j] / m_rho_r[i][j] - 2 * ry / m_rho_u[i][j] / m_rho_d[i][j];
            m_s1[i][j] = 4 * rx / (1 + m_rho_l[i][j]) / m_rho_r[i][j] / (1 + m_rho_r[i][j]);
            m_s2[i][j] = 4 * rx / (1 + m_rho_r[i][j]) / m_rho_l[i][j] / (1 + m_rho_l[i][j]);
            m_s3[i][j] = 4 * ry / (1 + m_rho_d[i][j]) / m_rho_u[i][j] / (1 + m_rho_u[i][j]);
            m_s4[i][j] = 4 * ry / (1 + m_rho_u[i][j]) / m_rho_d[i][j] / (1 + m_rho_d[i][j]);

            // TODO: INITIAL CONDITION 0; ADD THIS SOMEWHERE ELSE TO ALLOW FOR OTHER ICs
            m_un[i][j] = 0.0d;
            if (m_grid->is_on_shape(i, j))
                m_un[i][j] = m_grid->get_surface_value(i, j);
            else if (m_grid->is_in_shape(i, j))
                m_un[i][j] = 1.0d;
        }
    }
};


void Scheme::iteration(int n) {
    int m;
    size_t i, j;
    size_t nx = m_grid->get_nx();
    size_t ny = m_grid->get_ny();

    for (m=0; m < n; ++m) {
#pragma omp parallel for private(i, j)
        for (i=0; i < nx + 1; ++i) {
            for (j=0; j < ny + 1; ++j) {
                // BCs on shape surface
                if (m_grid->is_on_shape(i, j) || m_grid->is_in_shape(i, j))
                    m_unp1[i][j] = m_un[i][j];

                // BCs on edges
                else if (i == 0)
                    m_unp1[i][j] = 0;
                else if (i == nx)
                    m_unp1[i][j] = 0;
                else if (j == 0)
                    m_unp1[i][j] = 0;
                else if (j == ny)
                    m_unp1[i][j] = 0;

                // Numerical scheme
                else {
                    m_unp1[i][j] = m_s0[i][j] * m_un[i][j]
                        + m_s1[i][j] * m_un[i+1][j] * m_un_rho_r[i][j]
                        + m_s2[i][j] * m_un[i-1][j] * m_un_rho_l[i][j]
                        + m_s3[i][j] * m_un[i][j+1] * m_un_rho_u[i][j]
                        + m_s4[i][j] * m_un[i][j-1] * m_un_rho_d[i][j];
                }
            }
        }

        // Swap pointers
        std::swap(m_un, m_unp1);

        ++m_iterations;
        if (m_iterations % 1000 == 0) {
            std::cout << "Iteration " << m_iterations << std::endl;
        }
    }
};
