#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <stdlib.h>

#include "grid.hpp"
#include "rootfinder.hpp"


Grid::Grid(double xmin, double xmax, int nx, double ymin, double ymax, int ny) {
    m_xmin = xmin;
    m_xmax = xmax;
    m_nx = nx;

    m_ymin = ymin;
    m_ymax = ymax;
    m_ny = ny;

    m_dx = (m_xmax - m_xmin) / m_nx;
    m_dy = (m_ymax - m_ymin) / m_ny;

    // Allocate and fill axes
    m_x.resize(boost::extents[m_nx + 1]);
    for (size_t i = 0; i < m_nx + 1; ++i) {
        m_x[i] = m_xmin + i * m_dx;
    }

    m_y.resize(boost::extents[m_ny + 1]);
    for (size_t i = 0; i < m_ny + 1; ++i) {
        m_y[i] = m_ymin + i * m_dy;
    }

    allocate_2d_array(&m_on_shape);
    allocate_2d_array(&m_in_shape);

    for (const auto &xi : m_x) {
        m_surf_x[xi] = std::vector<double>();
    }
    for (const auto &yi : m_y) {
        m_surf_y[yi] = std::vector<double>();
    }
};


void Grid::set_shape(Shape *shape) {
    m_shape = shape;

    // Recalculate surface
    calculate_surface_points();
    set_surface_points_x();
    set_surface_points_y();

    // Recalculate bools
    set_bools_in_shape();
    set_bools_on_shape();
};


void Grid::free_points() {
    for (auto p : m_points) {
        delete p;
    }
    m_points.clear();
}


void Grid::calculate_surface_points() {
    std::cout << "Calculating surface points\n";

    free_points();

    double smin = m_shape->get_smin();
    double smax = m_shape->get_smax();

    // Set the first point
    double sn = smin;
    double xn = m_shape->f_x(sn);
    double yn = m_shape->f_y(sn);
    point *pt = new point(xn, yn, sn);
    m_points.push_back(pt);

    // Set the remaining points
    double xerr = m_dx - m_atol;
    double yerr = m_dy - m_atol;
    while (sn < smax) {
        double xref = m_points.back()->x;
        double yref = m_points.back()->y;

        double t_ds = m_shape->get_ds();
        xn = m_shape->f_x(sn + t_ds);
        yn = m_shape->f_y(sn + t_ds);
        while (std::abs(xn - xref) > xerr || (std::abs(yn - yref) > yerr)) {
            t_ds = t_ds / 2;
            xn = m_shape->f_x(sn + t_ds);
            yn = m_shape->f_y(sn + t_ds);
        }
        sn = sn + t_ds;
        pt = new point(xn, yn, sn);
        m_points.push_back(pt);
    }
};


void Grid::set_surface_points_x() {
    std::cout << "Setting surface points x\n";

    for (const auto &xi : m_x) {
        m_surf_x[xi].clear();
    }

    for (size_t i = 0; i < m_points.size() - 1; ++i) {
        point *pa = m_points[i];
        point *pb = m_points[i + 1];
        double xa = pa->x;
        double ya = pa->y;
        double sa = pa->s;
        double xb = pb->x;
        double yb = pb->y;
        double sb = pb->s;

        if (std::abs(xa - xb) < m_atol)
            continue;

        int idx = -1;
        for (size_t j = 0; j < m_nx + 1; ++j) {
            bool b1 = ((xa + m_atol) > m_x[j]) & ((xb - m_atol) < m_x[j]);
            bool b2 = ((xa - m_atol) < m_x[j]) & ((xb + m_atol) > m_x[j]);
            if (b1 || b2)
                idx = j;
        }
        if (idx == -1)
            continue;

        double xj = m_x[idx];
        double yj = root_finder_x(m_shape, xj, ya, sa);
        if (!std::isnan(yj)) {
            m_surf_x[xj].push_back(yj);
            continue;
        }

        yj = root_finder_x(m_shape, xj, ya, sb);
        if (!std::isnan(yj)) {
            m_surf_x[xj].push_back(yj);
            continue;
        }

        yj = root_finder_x(m_shape, xj, yb, sa);
        if (!std::isnan(yj)) {
            m_surf_x[xj].push_back(yj);
            continue;
        }

        yj = root_finder_x(m_shape, xj, yb, sb);
        if (!std::isnan(yj)) {
            m_surf_x[xj].push_back(yj);
            continue;
        }

        std::cout << "No solution found!?\n";
        exit(EXIT_FAILURE);
    }
};


void Grid::set_surface_points_y() {
    std::cout << "Setting surface points y\n";

    for (const auto &yi : m_y) {
        m_surf_y[yi].clear();
    }

    for (size_t i = 0; i < m_points.size() - 1; ++i) {
        point *pa = m_points[i];
        point *pb = m_points[i + 1];
        double xa = pa->x;
        double ya = pa->y;
        double sa = pa->s;
        double xb = pb->x;
        double yb = pb->y;
        double sb = pb->s;

        if (std::abs(ya - yb) < m_atol)
            continue;

        int idx = -1;
        for (size_t j = 0; j < m_ny + 1; ++j) {
            bool b1 = ((ya + m_atol) > m_y[j]) & ((yb - m_atol) < m_y[j]);
            bool b2 = ((ya - m_atol) < m_y[j]) & ((yb + m_atol) > m_y[j]);
            if (b1 || b2)
                idx = j;
        }
        if (idx == -1)
            continue;

        double yj = m_y[idx];
        double xj = root_finder_y(m_shape, xa, yj, sa);
        if (!std::isnan(xj)) {
            m_surf_y[yj].push_back(xj);
            continue;
        }

        xj = root_finder_y(m_shape, xa, yj, sb);
        if (!std::isnan(xj)) {
            m_surf_y[yj].push_back(xj);
            continue;
        }

        xj = root_finder_y(m_shape, xb, yj, sa);
        if (!std::isnan(xj)) {
            m_surf_y[yj].push_back(xj);
            continue;
        }

        xj = root_finder_y(m_shape, xb, yj, sb);
        if (!std::isnan(xj)) {
            m_surf_y[yj].push_back(xj);
            continue;
        }

        std::cout << "No solution found!?\n";
        exit(EXIT_FAILURE);
    }
};


void Grid::set_bools_in_shape() {
    std::cout << "Settings bools in shape\n";

    // Initialize
    for (size_t j = 0; j < m_ny + 1; ++j) {
        for (size_t i = 0; i < m_nx + 1; ++i) {
            m_in_shape[i][j] = false;
        }
    }

    double min_elem;
    double max_elem;
    for (size_t j = 0; j < m_ny + 1; ++j) {
        double yj = m_y[j];
        for (size_t i = 0; i < m_nx + 1; ++i) {
            double xi = m_x[i];
            bool out_x = true;
            bool out_y = true;
            if (m_surf_x[xi].size() > 0) {
                min_elem = std::min_element(m_surf_x[xi].begin(), m_surf_x[xi].end())[0];
                max_elem = std::max_element(m_surf_x[xi].begin(), m_surf_x[xi].end())[0];
                out_x = (yj < min_elem + m_atol) || (yj > max_elem - m_atol);
            }

            if (m_surf_y[yj].size() > 0) {
                min_elem = std::min_element(m_surf_y[yj].begin(), m_surf_y[yj].end())[0];
                max_elem = std::max_element(m_surf_y[yj].begin(), m_surf_y[yj].end())[0];
                out_y = (xi < min_elem + m_atol) || (xi > max_elem - m_atol);
            }

            m_in_shape[i][j] = !(out_x || out_y);
        }
    }
};


void Grid::set_bools_on_shape() {
    std::cout << "Settings bools on shape\n";

    // Initialize
    for (size_t i = 0; i < m_nx + 1; ++i) {
        for (size_t j = 0; j < m_ny + 1; ++j) {
            m_on_shape[i][j] = false;
        }
    }

    bool b;
    int idx;
    for (size_t i = 0; i < m_nx + 1; ++i) {
        double xi = m_x[i];
        for (const auto &yi: m_surf_x[xi]) {
            idx = -1;
            for (size_t j = 0; j < m_ny + 1; ++j) {
                b = (yi - m_atol < m_y[j]) & (yi + m_atol > m_y[j]);
                if (b)
                    idx = j;
            }
            if (idx == -1)
                continue;

            m_on_shape[i][idx] = !m_in_shape[i][idx];
        }
    }

    for (size_t i = 0; i < m_ny + 1; ++i) {
        double yi = m_y[i];
        for (const auto &xi: m_surf_y[yi]) {
            idx = -1;
            for (size_t j = 0; j < m_nx + 1; ++j) {
                b = (xi - m_atol < m_x[j]) & (xi + m_atol > m_x[j]);
                if (b)
                    idx = j;
            }
            if (idx == -1)
                continue;

            m_on_shape[idx][i] = !m_in_shape[idx][i];
        }
    }
};


void Grid::check_in_bounds_x(size_t i) const {
    if (i > m_nx) {
        std::cout << "x-index i=" << i << " exceeds max nx=" << m_nx << "\n";
        exit(EXIT_FAILURE);
    }
}


void Grid::check_in_bounds_y(size_t j) const {
    if (j > m_ny) {
        std::cout << "y-index j=" << j << " exceeds max ny=" << m_ny << "\n";
        exit(EXIT_FAILURE);
    }
}


bool Grid::is_in_shape(size_t i, size_t j) const {
    check_in_bounds_x(i);
    check_in_bounds_y(j);
    return m_in_shape[i][j];
};


bool Grid::is_on_shape(size_t i, size_t j) const {
    check_in_bounds_x(i);
    check_in_bounds_y(j);
    return m_on_shape[i][j];
}


bool Grid::is_on_shape(double x, double y) const {
    bool on_shape = false;

    // Check y values at x
    if (m_surf_x.count(x)) {
        std::vector<double> arr = m_surf_x.find(x)->second;
        for (double yi: arr) {
            if (std::abs(y - yi) < m_atol)
                on_shape = true;
        }
    }

    // Check x values at y
    if (m_surf_y.count(y)) {
        std::vector<double> arr = m_surf_y.find(y)->second;
        for (double xi: arr) {
            if (std::abs(x - xi) < m_atol)
                on_shape = true;
        }
    }

    return on_shape;
}


double Grid::get_x(int i) const {
    check_in_bounds_x(i);
    return m_x[i];
}


double Grid::get_y(int j) const {
    check_in_bounds_y(j);
    return m_y[j];
}


int Grid::get_max_index_x(double x) const {
    int j = -1;
    for (size_t k = 0; k < m_nx + 1; ++k) {
        if (m_x[k] > x + m_atol) {
            j = k - 1;
            break;
        }
    }
    if (j == -1) {
        std::cout << "max index x cannot be found when x=" << x << "\n";
        exit(EXIT_FAILURE);
    }
    return j;
};


int Grid::get_max_index_y(double y) const {
    int j = -1;
    for (size_t k = 0; k < m_ny + 1; ++k) {
        if (m_y[k] > y + m_atol) {
            j = k - 1;
            break;
        }
    }
    if (j == -1) {
        std::cout << "max index y cannot be found when y=" << y << "\n";
        exit(EXIT_FAILURE);
    }
    return j;
};


double Grid::get_surface_value(double x, double y) const {
    if (!is_on_shape(x, y)) {
        std::cout << "get_surface_value(dou, dou)::point x=" << x << ", y=" << y << " is not on the shape surface!\n";
        exit(EXIT_FAILURE);
    }
    return m_shape->boundary_conditions(x, y);
}


double Grid::get_surface_value(size_t i, size_t j) const {
    double x = get_x(i);
    double y = get_y(j);
    if (!is_on_shape(i, j)) {
        std::cout << "get_surface_value(int, int)::point x=" << x << ", y=" << y << " is not on the shape surface!\n";
        exit(EXIT_FAILURE);
    }
    return m_shape->boundary_conditions(x, y);
}


std::vector<double> Grid::get_surf_x(size_t i) const {
    double x = get_x(i);
    return m_surf_x.at(x);
}


std::vector<double> Grid::get_surf_y(size_t j) const {
    double y = get_y(j);
    return m_surf_y.at(y);
}
