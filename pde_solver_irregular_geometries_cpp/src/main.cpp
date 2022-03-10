#pragma once
#include "../include/solver.h"
#include <iostream>

int main(int argc, char *argv[])
{
    std::cout << "\nCircle shape" << std::endl;
    Circle s1;
    Solver<Circle> s_circle = Solver(s1, 100, 0.000001d, 1, 1e-6d);
    s_circle.iteration(10000);
    s_circle.dump_solution("circle.txt");

    std::cout << "\nSquare shape" << std::endl;
    Square s2;
    Solver<Square> s_square = Solver(s2, 100, 0.000001d, 1, 1e-4d);
    s_square.iteration(25000);
    s_square.dump_solution("square.txt");

    std::cout << "\nStar shape" << std::endl;
    Star s3;
    Solver<Star> s_star = Solver(s3, 200, 0.000001d, 1, 1e-6d);
    s_star.iteration(75000);
    s_star.dump_solution("star.txt");

};
