import numpy as np
from solver import Solver

class Square(Solver):

    def boundary_conditions(self, x, y):
        return np.sin(np.pi * (x-0.5) / 0.4) \
            * np.sin(np.pi * (y-0.5) / 0.4)

    def initial_conditions(self):
        return np.zeros(self.X.shape)

    def f_x(self, s):
        xmin = 0.3
        xmax = 0.7
        if s < 1:
            return xmin
        elif s < 2:
            return xmin + (s-1)*(xmax-xmin)
        elif s < 3:
            return xmax
        else:
            return xmax + (s-3)*(xmin-xmax)

    def f_y(self, s):
        ymin = 0.3
        ymax = 0.7
        if s < 1:
            return ymin + s*(ymax-ymin)
        elif s < 2:
            return ymax
        elif s < 3:
            return ymax + (s-2)*(ymin-ymax)
        else:
            return ymin

def main():
    q = Square(101, 0.000001, 0, 4, 100, atol=1e-4)
    q.iteration(100000)
    q.plot_solution('square.png')

if __name__=='__main__':
    main()
