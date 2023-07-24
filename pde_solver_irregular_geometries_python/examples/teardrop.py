import numpy as np
import sys
sys.path.append('..')

from solver import Solver

class Star(Solver):

    def boundary_conditions(self, x, y):
        return np.asarray(np.sin(np.pi*y))

    def initial_conditions(self):
        return np.zeros(self.X.shape)

    def f_x(self, s):
        return np.cos(s)

    def f_y(self, s):
        return np.sin(s) * np.sin(s/2)**3

def main():
    q = Star(101, 0.000001, 0, 2*np.pi, 100, -3, 3, -3, 3)
    q.iteration(100000)
    q.plot_solution('teardrop.png')

if __name__=='__main__':
    main()
