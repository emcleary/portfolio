#include <boost/math/constants/constants.hpp>
#include <iostream>

#include "grid.hpp"
#include "scheme.hpp"
#include "shapes.hpp"
#include "solver.hpp"


class Circle : public Shape {
public:
    Circle() : Shape(0, 2*PI) {};

    double f_x(double s) const override {
        return 0.5d + std::cos(s) / 5.0;
    }

    double f_y(double s) const override {
        return 0.5d + std::sin(s) / 5.0;
    }

    double boundary_conditions(double x, double y) const override {
        return std::sin(4.0d * PI * (x - 0.5d));
    }
};


int main() {

    double xmin = 0.0;
    double xmax = 1.0;
    int nx = 100;

    double ymin = 0.0;
    double ymax = 1.0;
    int ny = 100;

    double dt = 0.000001d;
    Solver solver(dt);

    Grid grid(xmin, xmax, nx, ymin, ymax, ny);
    solver.set_grid(&grid);

    Circle circle;
    solver.set_shape(&circle);
    solver.solve(10000);
    solver.dump_solution("solution_circle.txt");
    
    return 0;
}
