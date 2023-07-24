
add_executable(ex01_circle example_01_circle.cpp)
target_link_libraries(ex01_circle PRIVATE LibSrc gsl cblas)

add_executable(ex02_star_square example_02_star_and_square.cpp)
target_link_libraries(ex02_star_square PRIVATE LibSrc gsl cblas)
