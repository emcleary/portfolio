import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from itertools import product
from time import time

class Solver:
    """A 2-dimensional, first order accurate finite difference
    scheme for elliptic PDE with irregular geometries

    Reference: "Numerical Partial Differential Equations: Conservation
    Laws and Elliptic Equations", Ch. 11, by J. W. Thomas

    Attributes
    ----------
    un : float array, 2-dimensional
        state vector representing the solution at timestep "n"
    unp1 : float array, 2-dimensional
        state vector representing the solution at timestep "n+1"
    t : float
        time, initialized to zero
    dt : float
        timestep
    iterations : int
        number of iterations done
    xmin : float
        minimum value of the x-axis
    xmax : float
        maximum value of the x-axis
    ymin : float
        minimum value of the y-axis
    ymax : float
        maximum value of the y-axis
    dx : float
        stepsize of the x-axis
    dy : float
        stepsize of the y-axis
    x : float array, 1-dimensional
        gridpoints on the x-axis
    y : float array, 1-dimensional
        gridpoints on the y-axis
    X : float array, 2-dimensional
        the x domain
    Y : float array, 2-dimensional
        the y domain
    si : float
        starting parameter
    sf : float
        final parameter
    ds : float
        max stepsize for the parameter
    atol : float
        tolerance
    rho_l : float array, 2-dimensional
        distance to the surface in the -x direction, 1 by default
    rho_r : float array, 2-dimensional
        distance to the surface in the +x direction, 1 by default
    rho_d : float array, 2-dimensional
        distance to the surface in the -y direction, 1 by default
    rho_u : float array, 2-dimensional
        distance to the surface in the +y direction, 1 by default
    un_rho_l : float array, 2-dimensional
        boundary condition of a neighboring surface in the -x direction
    un_rho_r : float array, 2-dimensional
        boundary condition of a neighboring surface in the +x direction
    un_rho_d : float array, 2-dimensional
        boundary condition of a neighboring surface in the -y direction
    un_rho_u : float array, 2-dimensional
        boundary condition of a neighboring surface in the +y direction
    s0 : float array, 2-dimensional
        array of scheme parameters
    s1 : float array, 2-dimensional
        array of scheme parameters
    s2 : float array, 2-dimensional
        array of scheme parameters
    s3 : float array, 2-dimensional
        array of scheme parameters
    s4 : float array, 2-dimensional
        array of scheme parameters
    b_on_shape : boolean array, 2-dimensional
        tracks gridpoints on the shape's surface
    b_in_shape : boolean array, 2-dimensional
        tracks gridpoints within the shape
    surf_x : dictionary
        tracks points on x-gridlines on the shape's surface
    surf_y : dictionary
        tracks points on y-gridlines on the shape's surface

    Methods
    -------
    f_x(s)
        parameteric function calculating x at parameter "s"
    f_y(s)
        parameteric function calculating y at parameter "s"
    initial_conditions()
        Sets initial conditions
    boundary_conditions()
        Calculate boundary conditions on the surface (edges of the
        domain are fixed to 0).
    iteration(n, tol)
        Runs "n" iterations of the scheme
    plot_solution(filename)
        Plots the solution and saves the figure to "filename".

    """

    def __init__(self, n, dt, si, sf, ns, xmin=0, xmax=1, ymin=0, ymax=1, atol=1e-6):
        """Constructs attributes for the solver object.

        Arguments
        ---------
        n : int
            number of gridpoints on both the x- and y-axis
        dt : float
            timestep
        si : float
            starting parameter
        sf : float
            final parameter
        ns : int
            number of gridcells
        xmin : float, optional
            minimum value of the x-axis
        xmax : float, optional
            maximum value of the x-axis
        ymin : float, optional
            minimum value of the y-axis
        ymax : float, optional
            maximum value of the y-axis
        atol : float, optional
            tolerance
        """

        self.atol = atol
        self.t = 0
        self.iterations = 0
        self.dt = dt

        # Parametric variable
        self.si = si
        self.sf = sf
        self.ds = (self.sf - self.si) / (ns+1)

        # Grid
        self.n = n
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.x = np.linspace(self.xmin, self.xmax, self.n)
        self.y = np.linspace(self.ymin, self.ymax, self.n)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        assert self.dx == self.dy, 'Only tested for cases with a regular grid'
        self.r = self.dt / self.dx / self.dx
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

        # Distance from gridpoint to shape surface, if adjacent; 1 by default
        self.rho_l = np.ones(self.X.shape)
        self.rho_r = np.ones(self.X.shape)
        self.rho_u = np.ones(self.X.shape)
        self.rho_d = np.ones(self.X.shape)

        # Overrides adjecent gridpoint with boundary condition of the surface, if needed;
        # 1 by default
        self.un_rho_l = np.ones(self.X.shape)
        self.un_rho_r = np.ones(self.X.shape)
        self.un_rho_u = np.ones(self.X.shape)
        self.un_rho_d = np.ones(self.X.shape)

        # Scheme parameters
        self.s0 = None
        self.s1 = None
        self.s2 = None
        self.s3 = None
        self.s4 = None

        # Booleans for gridpoints
        self.b_on_shape = np.zeros(self.X.shape, dtype=bool)
        self.b_in_shape = np.zeros(self.X.shape, dtype=bool)

        # State variables
        self.un = np.zeros(self.X.shape)
        self.unp1 = np.zeros(self.X.shape)

        # Surface points
        self.surf_x = None
        self.surf_y = None

        # Fill/update arrays with values needed
        self._calc_surf_points()
        self._set_bools_on_shape()
        self._set_bools_in_shape()
        self._calc_dist_from_surf()
        self._calc_scheme_parameters()


    def _set_bools_on_shape(self):
        """Determines which gridpoints are on the surfaces of the shape,
        filling the b_on_shape attribute.
        """
        def is_on_shape(a, s):
            b1 = a - self.atol < s
            b2 = a + self.atol > s
            return any(b1 * b2)

        for j in range(self.n):
            y = self.y[j]
            for i in range(self.n):
                x = self.x[i]
                self.b_on_shape[i,j] = is_on_shape(x, self.surf_y[y]) \
                    or is_on_shape(y, self.surf_x[x])


    def _set_bools_in_shape(self):
        """Determines which gridpoints are within the shape, filling the
        b_in_shape attribute. (Breath First Search)
        """

        queue = [(0, 0)]
        b_out_shape = np.zeros(self.X.shape, dtype=bool)

        def is_in_shape(c, d, surf):
            if c < d:
                b1 = c + self.atol < surf
                b2 = d - self.atol > surf
            else:
                b1 = c - self.atol > surf
                b2 = d + self.atol < surf
            return any(b1 * b2)

        while queue:
            i,j = queue.pop(0)  ## FIFO

            if self.b_in_shape[i,j]:
                continue
            if self.b_on_shape[i,j]:
                continue
            if b_out_shape[i,j]:
                continue
            b_out_shape[i,j] = True

            x = self.x[i]
            y = self.y[j]

            if i > 0:
                queue.append((i-1,j))
                n = self.x[i-1]
                surf = self.surf_y[y]
                if surf:
                    self.b_in_shape[i-1,j] = is_in_shape(x, n, surf)

            if i < self.n-1:
                queue.append((i+1,j))
                n = self.x[i+1]
                surf = self.surf_y[y]
                if surf:
                    self.b_in_shape[i+1,j] = is_in_shape(x, n, surf)

            if j > 0:
                queue.append((i,j-1))
                n = self.y[j-1]
                surf = self.surf_x[x]
                if surf:
                    self.b_in_shape[i,j-1] = is_in_shape(y, n, surf)

            if j < self.n-1:
                queue.append((i,j+1))
                n = self.y[j+1]
                surf = self.surf_x[x]
                if surf:
                    self.b_in_shape[i,j+1] = is_in_shape(y, n, surf)

        # Update b_in_shape
        self.b_in_shape = ~b_out_shape * ~self.b_on_shape


    def _calc_surf_points(self):
        """ Calculates points on the surface of the shape along each gridline.
        """

        # Setup a bunch of points by following the trajectory;
        # each x- or y-line must have a pair of points surrounding it
        sn = self.si
        xn = self.f_x(sn)
        yn = self.f_y(sn)
        points = [(xn, yn, sn)]
        while sn < self.sf:
            ds = self.ds
            xn = self.f_x(sn + self.ds)
            yn = self.f_y(sn + self.ds)
            while np.abs(xn - points[-1][0]) > (self.dx - self.atol) or np.abs(yn - points[-1][1]) > (self.dy - self.atol):
                ds /= 2
                xn = self.f_x(sn + ds)
                yn = self.f_y(sn + ds)
            sn += ds
            points.append((xn, yn, sn))

        def equations_x(p, y):
            x, s = p
            return (y-self.f_y(s), x-self.f_x(s))

        def equations_y(p, x):
            y, s = p
            return (y-self.f_y(s), x-self.f_x(s))

        # Use these points to find intersections with the grid on each axis
        self.surf_x = {xi:[] for xi in self.x}
        for pa, pb in zip(points[1:], points[:-1]):
            x_pa, y_pa, s_pa = pa
            x_pb, y_pb, s_pb = pb
            if np.isclose(x_pa, x_pb, self.atol):
                continue

            # Find a gridline of x
            if x_pa > x_pb:
                b = (x_pa + self.atol > self.x) * (x_pb - self.atol < self.x)
            else:
                b = (x_pa - self.atol < self.x) * (x_pb + self.atol > self.x)

            if any(b):
                assert sum(b) == 1
                x = self.x[b][0]
            else:
                continue

            # Calculate the y not on a gridline
            y = None
            for yi, si in product([y_pa, y_pb], [s_pa, s_pb]):
                r = root(equations_y, [yi, si], x)
                if r.success:
                    y = r.x[0]
                    self.surf_x[x].append(y)
                    break
            assert y is not None
        
        self.surf_y = {yi:[] for yi in self.y}
        for pa, pb in zip(points[1:], points[:-1]):
            x_pa, y_pa, s_pa = pa
            x_pb, y_pb, s_pb = pb
            if np.isclose(y_pa, y_pb, self.atol):
                continue

            # Find a gridline of y
            if y_pa > y_pb:
                b = (y_pa + self.atol > self.y) * (y_pb - self.atol < self.y)
            else:
                b = (y_pa - self.atol < self.y) * (y_pb + self.atol > self.y)

            if any(b):
                assert sum(b) == 1
                y = self.y[b][0]
            else:
                continue

            # Calculate the x not on a gridline
            x = None
            for xi, si in product([x_pa, x_pb], [s_pa, s_pb]):
                r = root(equations_x, [xi, si], y)
                if r.success:
                    x = r.x[0]
                    self.surf_y[y].append(x)
                    break
            assert x is not None


    def _calc_dist_from_surf(self):
        """Calculated distances from the shape surface to neighboring
        gridpoints, filling the rho_* and un_rho_* attributes.
        """
        
        # Find distances from grid points to surface points,
        # b_in_shape will be used to identify directions to set the measurement.

        # Loop over surface points
        for idx_x, (x, ys) in enumerate(self.surf_x.items()):
            for y in ys:
                gt = self.y > y + self.atol
                idx_y = np.where(gt)[0][0] - 1
                if self.b_on_shape[idx_x,idx_y]:
                    continue
                if ~self.b_in_shape[idx_x,idx_y] and self.b_in_shape[idx_x,idx_y+1]:
                    d = y - self.y[idx_y]
                    self.rho_u[idx_x,idx_y] = d / self.dy
                    self.un_rho_u[idx_x,idx_y] = self.boundary_conditions(self.x[idx_x], self.y[idx_y]+d)
                elif ~self.b_in_shape[idx_x,idx_y+1] and self.b_in_shape[idx_x,idx_y]:
                    d = self.y[idx_y+1] - y
                    self.rho_d[idx_x,idx_y+1] = d / self.dy
                    self.un_rho_d[idx_x,idx_y+1] = self.boundary_conditions(self.x[idx_x], self.y[idx_y+1]-d)

        for idx_y, (y, xs) in enumerate(self.surf_y.items()):
            for x in xs:
                gt = self.x > x + self.atol
                idx_x = np.where(gt)[0][0] - 1
                if self.b_on_shape[idx_x,idx_y]:
                    continue
                if ~self.b_in_shape[idx_x,idx_y] and self.b_in_shape[idx_x+1,idx_y]:
                    d = x - self.x[idx_x]
                    self.rho_r[idx_x,idx_y] = d / self.dx
                    self.un_rho_r[idx_x,idx_y] = self.boundary_conditions(self.x[idx_x]+d, self.y[idx_y])
                elif self.b_in_shape[idx_x,idx_y] and ~self.b_in_shape[idx_x+1,idx_y]:
                    d = self.x[idx_x+1] - x
                    self.rho_l[idx_x+1,idx_y] = d / self.dx
                    self.un_rho_l[idx_x+1,idx_y] = self.boundary_conditions(self.x[idx_x+1]-d, self.y[idx_y])


    def _calc_scheme_parameters(self):
        """Calculated constant parameters used in the scheme and the initial
        state.
        """

        self.s0 = 1 - 2*self.r/self.rho_l/self.rho_r - 2*self.r/self.rho_u/self.rho_d
        self.s1 = 4 * self.r / (1+self.rho_l) / self.rho_r / (1+self.rho_r)
        self.s2 = 4 * self.r / (1+self.rho_r) / self.rho_l / (1+self.rho_l)
        self.s3 = 4 * self.r / (1+self.rho_d) / self.rho_u / (1+self.rho_u)
        self.s4 = 4 * self.r / (1+self.rho_u) / self.rho_d / (1+self.rho_d)

        self.un = self.initial_conditions()
        self.un[self.b_on_shape] = self.boundary_conditions(self.X[self.b_on_shape], self.Y[self.b_on_shape])
        self.un[self.b_in_shape] = 1.0

    def iteration(self, n=100, tol=1e-4):
        """Performs n iterations for the first-order accurate scheme.
        
        Arguments
        ---------
        n : int, optional
            number of iterations
        tol : float, optional
            terminates iterations when this error tolerance between iterations
            is met

        """
        for _ in range(n):
            self.iterations += 1
            self.t += self.dt

            # Calculate state at the next timestep
            self.unp1 = self.s0 * self.un
            self.unp1[:-1,:] += self.s1[:-1,:] * self.un[1:,:] * self.un_rho_r[:-1,:]
            self.unp1[1:,:] += self.s2[1:,:] * self.un[:-1,:] * self.un_rho_l[1:,:]
            self.unp1[:,:-1] += self.s3[:,:-1] * self.un[:,1:] * self.un_rho_u[:,:-1]
            self.unp1[:,1:] += self.s4[:,1:] * self.un[:,:-1] * self.un_rho_d[:,1:]

            # Preserve boundary conditions on the shape and parameters inside the shape
            self.unp1[self.b_on_shape] = self.un[self.b_on_shape]
            self.unp1[self.b_in_shape] = self.un[self.b_in_shape]

            # BCs on the edges fixed at 0
            self.unp1[0,:] = 0
            self.unp1[-1,:] = 0
            self.unp1[:,0] = 0
            self.unp1[:,-1] = 0

            # Swap pointers
            self.un, self.unp1 = self.unp1, self.un

            # Calculate error between iterations
            error = np.linalg.norm(self.un - self.unp1)
            if self.iterations % 1000 == 0:
                print(f'iteration {self.iterations}, error {error}')
            if error < tol:
                break

        print(f'finished on iteration {self.iterations} with error {error}')
        return self.un

    def boundary_conditions(self, x, y):
        """Set boundary conditions of the surface of the shape."""
        return np.sin(4*np.pi*(x - 0.5))

    def initial_conditions(self):
        """Set initial condition of the whole state array."""
        return np.zeros(self.X.shape)

    def f_x(self, s):
        """Parametric function for the x-coordinate of the shape."""
        return 0.5 + np.cos(s)/5

    def f_y(self, s):
        """Parametric function for the y-coordinate of the shape."""
        return 0.5 + np.sin(s)/5

    def plot_solution(self, filename=None):
        """Plots the solution"""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        c = np.copy(self.un)
        c[self.b_in_shape] = np.nan
        cs = ax.contourf(self.X, self.Y, c)
        cbar = fig.colorbar(cs)
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename, bbox_inches='tight')

    def dump_solution(self, filename):
        """Writes solution to a file"""
        with open(filename, 'w') as file:
            file.write(f'grid {self.n} {self.n}\n')
            for j in range(self.n):
                for i in range(self.n):
                    xi = self.x[i]
                    yj = self.y[j]
                    if self.b_in_shape[i][j]:
                        file.write(f'{xi:.6e} {yj:.6e} NAN\n')
                    else:
                        u = self.un[i, j]
                        file.write(f'{xi:.6e} {yj:.6e} {u:.6e}\n')


def main():
    s = Solver(101, 0.000001, 0, 2*np.pi, 1)
    s.iteration(10000)
    s.plot_solution('solver.png')


if __name__=='__main__':
    main()
