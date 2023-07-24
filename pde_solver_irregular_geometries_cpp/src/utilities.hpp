#ifndef UTILITIES_H
#define UTILITIES_H

#include <boost/multi_array.hpp>

typedef boost::multi_array<double, 2> array_double_2d;
typedef boost::multi_array<double, 1> array_double_1d;
typedef boost::multi_array<bool, 2> array_bool_2d;
const double PI = boost::math::constants::pi<double>();

#endif
