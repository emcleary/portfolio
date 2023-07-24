#ifndef GRID_H
#define GRID_H

#include "shapes.hpp"
#include "utilities.hpp"
#include <map>



/// @brief A 2-dimensional uniform grid used for the PDE solver
class Grid {

public:
    /// @brief A 2-dimensional uniform grid used for the PDE solver
    /// @param xmin Minimum value of x
    /// @param xmax Maximum value of x
    /// @param nx Number of grid cells in x
    /// @param ymin Miminum value of y
    /// @param ymax Maximum value of y
    /// @param ny Number of grid cells in y
    Grid(double xmin, double xmax, int nx, double ymin, double ymax, int ny);

    ~Grid() {
        free_points();
    }

    /// @brief Set shape to the grid.
    /// @param shape An instance of a Shape child class
    void set_shape(Shape *shape);

    /// @brief Get the shape from the grid.
    /// @return A pointer to a Shape child class used by the grid
    Shape *get_shape() { return m_shape; };

    /// @brief Check if a grid point is inside the shape.
    /// @param i Index in the x-coordinate
    /// @param j Index in the y-coordinate
    /// @return true if the grid point is within the shape, false otherwise
    bool is_in_shape(size_t i, size_t j) const;

    /// @brief Check if a grid point is on the shape.
    /// @param i Index in the x-coordinate
    /// @param j Index in the y-coordinate
    /// @return true if the grid point is on the shape, false otherwise
    bool is_on_shape(size_t i, size_t j) const;

    /// @brief Check if a point is inside the shape.
    /// @param x x-coordinate within the x bounds
    /// @param y y-coordinate within the y bounds
    /// @return true if the grid point is on the shape, false otherwise
    bool is_on_shape(double x, double y) const;

    /// @brief Get x-coordinate at a given x index.
    /// @param i Index in the x-coordinate
    /// @return value of x at index i
    double get_x(int i) const;

    /// @brief Get y-coordinate at a given y index
    /// @param j Index in the y-coordinate
    /// @return value of y at index j
    double get_y(int j) const;

    /// @brief Get first x-coordinate index greater than the given value.
    /// @param x x-coordinate within the x bounds
    /// @return index from the x-coordinate
    int get_max_index_x(double x) const;

    /// @brief Get first y-coordinate index greater than the given value.
    /// @param y y-coordinate within the y bounds
    /// @return index from the y-coordinate
    int get_max_index_y(double y) const;

    /// @brief Get the number of grid cells in the x dimension
    /// @return Number of grid cells in the x dimension
    size_t get_nx() { return m_nx; };

    /// @brief Get the number of grid cells in the y dimension
    /// @return Number of grid cells in the y dimension
    size_t get_ny() { return m_ny; };

    /// @brief Get the grid cell length in the x dimension
    /// @return Length of grid cells in the x dimension
    double get_dx() { return m_dx; };

    /// @brief Get the grid cell length in the y dimension
    /// @return Length of grid cells in the y dimension
    double get_dy() { return m_dy; };

    /// @brief Get a vector of y-values that intersect the shape surface at value x at index i
    /// @param i Index in the x-coordinate
    /// @return vector of y-values at x corresponding to index i
    std::vector<double> get_surf_x(size_t i) const;

    /// @brief Get a vector of x-values that intersect the shape surface at value y at index j
    /// @param j Index in the y-coordinate
    /// @return vector of x-values at y corresponding to index j
    std::vector<double> get_surf_y(size_t j) const;

    /// @brief Calculate the boundary condition at a grid point on the shape surface
    /// @param x x-coordinate of a grid point
    /// @param y y-coordinate of a grid point
    /// @return value at the grid point on the surface
    double get_surface_value(double x, double y) const;

    /// @brief Calculate the boundary condition at a grid point on the shape surface
    /// @param i Index of the x-coordinate of a grid point
    /// @param j Index of the y-coordinate of a grid point
    /// @return value at the grid point on the surface
    double get_surface_value(size_t i, size_t j) const;

    /// @brief Resize/reshape a 2d array to the size of the grid
    /// @tparam T 2-d boost::multi_array of any type
    /// @param array Pointer to the array to be resized/reshaped
    template<typename T>
    void allocate_2d_array(T *array);


private:
    /// @brief Calculate a set of shape surface points
    void calculate_surface_points();
    /// @brief Calculate and store a set of grid points on the shape surface
    void set_surface_points_x();
    /// @brief Calculate and store a set of grid points on the shape surface
    void set_surface_points_y();
    /// @brief Store bools representing whether each grid point is inside or outside the shape
    void set_bools_in_shape();
    /// @brief Store bools representing whether each grid point is on the shape surface or not
    void set_bools_on_shape();
    /// @brief Delete surface points
    void free_points();

    /// @brief Check that the x-coordinate index is in bounds
    /// @param i Index of x-coordinate
    /// @throw EXIT_FAILURE if the index is out of bounds
    void check_in_bounds_x(size_t i) const;

    /// @brief Check that the y-coordinate index is in bounds
    /// @param j Index of y-coordinate
    /// @throw EXIT_FAILURE if the index is out of bounds
    void check_in_bounds_y(size_t j) const;


    /// @brief Object pointer representing the shape in the domain
    Shape *m_shape;
    /// @brief Minimum value of x-coordinate
    double m_xmin;
    /// @brief Maximum value of x-coordinate
    double m_xmax;
    /// @brief Minimum value of y-coordinate
    double m_ymin;
    /// @brief Maximum value of y-coordinate
    double m_ymax;
    /// @brief Number of grid cells in the x-coordinate
    size_t m_nx;
    /// @brief Number of grid cells in the y-coordinate
    size_t m_ny;
    /// @brief Length of grid cell in the x-coordinate
    double m_dx;
    /// @brief Length of grid cell in the y-coordinate
    double m_dy;
    /// @brief Array of grid points in the x-coordinate (size m_nx + 1)
    array_double_1d m_x;
    /// @brief Array of grid points in the x-coordinate (size m_ny + 1)
    array_double_1d m_y;
    /// @brief Array of bools representing grid points on the shape surface
    array_bool_2d m_on_shape;
    /// @brief Array of bools representing grid points inside the shape
    array_bool_2d m_in_shape;
    /// @brief Stores points in y that intersect the shape surface at each value of the x-coordinate
    std::map<double, std::vector<double>> m_surf_x;
    /// @brief Stores points in x that intersect the shape surface at each value of the y-coordinate
    std::map<double, std::vector<double>> m_surf_y;
    /// @brief Points on the shape surface found when traversing a parameterized function
    std::vector<point*> m_points;
    /// @brief Tolerance for comparing values
    double m_atol = 1e-6;

};


template<typename T>
void Grid::allocate_2d_array(T *array) {
    array->resize(boost::extents[m_nx + 1][m_ny + 1]);
}

#endif
