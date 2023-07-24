import numpy as np
import sys
sys.path.append('..')

from solver import Solver

class Complex(Solver):

    def boundary_conditions(self, x, y):
        return x*x + y*y

    def initial_conditions(self):
        return np.zeros(self.X.shape)

    def f_x(self, s):
        return (2*np.sin(s) + (2./3) * np.sin(2*s) + (5./3) * np.sin(3*s)) / 1.5

    def f_y(self, s):
        return (np.sin(s) + 2*np.sin(5*s) + 8*np.sin(3*s)) / 6

def main():
    q = Complex(201, 0.0000001, 0, 2*np.pi, 1000, -3, 3, -3, 3)
    q.iteration(10000)
    q.plot_solution('complex.png')

if __name__=='__main__':
    main()
