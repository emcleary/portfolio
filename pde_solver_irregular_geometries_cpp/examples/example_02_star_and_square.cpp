#include <boost/math/constants/constants.hpp>
#include <iostream>

#include "grid.hpp"
#include "scheme.hpp"
#include "shapes.hpp"
#include "solver.hpp"


class Star : public Shape {
public:
    Star() : Shape(0, 2*PI) {};

    double f_x(double s) const override {
        return (-9.0 * std::sin(2 * s) - 5.0 * std::sin(3 * s)) / 8.0;
    }

    double f_y(double s) const override {
        return ( 9.0 * std::cos(2 * s) - 5.0 * std::cos(3 * s)) / 8.0;
    }

    double boundary_conditions(double x, double y) const override {
        return std::log10(x * x + y * y);
    }
};


class Square : public Shape {
public:
    Square() : Shape(0, 4) {};

    double f_x(double s) const override {
        double dx = m_xmax - m_xmin;
        if (s < 1)
            return m_xmin;
        else if (s < 2)
            return m_xmin + (s - 1) * dx;
        else if (s < 3)
            return m_xmax;
        else
            return m_xmax - (s - 3) * dx;
    }

    double f_y(double s) const override {
        double dy = m_ymax - m_ymin;
        if (s < 1)
            return m_ymin + s * dy;
        else if (s < 2)
            return m_ymax;
        else if (s < 3)
            return m_ymax - (s - 2) * dy;
        else
            return m_ymin;
    }

    double boundary_conditions(double x, double y) const override {
        double xmid = 0.5 * (m_xmin + m_xmax);
        double ymid = 0.5 * (m_ymin + m_ymax);
        double mx = (x - xmid) / (m_xmax - m_xmin);
        double my = (y - ymid) / (m_ymax - m_ymin);
        return std::sin(PI * mx) * std::sin(PI * my);
    }

private:
    double m_xmin = -2.0;
    double m_xmax = 2.0;
    double m_ymin = -2.0;
    double m_ymax = 2.0;
};


int main() {

    double xmin = -4.0;
    double xmax = 4.0;
    int nx = 200;

    double ymin = -4.0;
    double ymax = 4.0;
    int ny = 200;

    double dt = 0.000001d;
    Solver solver(dt);

    Grid grid(xmin, xmax, nx, ymin, ymax, ny);
    solver.set_grid(&grid);

    // Example of using the builder design patter with shapes.
    // Simply set a new shape to the solver and run.
    // This can similarly be done with different grid objects
    // if the domain must change between runs.
    Star star;
    solver.set_shape(&star);
    solver.solve(50000);
    solver.dump_solution("solution_star.txt");

    Square square;
    solver.set_shape(&square);
    solver.solve(50000);
    solver.dump_solution("solution_square.txt");
    
    return 0;
}
