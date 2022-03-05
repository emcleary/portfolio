import numpy as np
from solver import Solver

class Star(Solver):

    def boundary_conditions(self, x, y):
        return np.log(x*x + y*y)

    def initial_conditions(self):
        return np.zeros(self.X.shape)

    def f_x(self, s):
        return (-9*np.sin(2*s) - 5*np.sin(3*s)) / 8

    def f_y(self, s):
        return (9*np.cos(2*s) - 5*np.cos(3*s)) / 8

def main():
    q = Star(201, 0.000001, 0, 2*np.pi, 100, -4, 4, -4, 4)
    q.iteration(100000)
    q.plot_solution('star.png')

if __name__=='__main__':
    main()
