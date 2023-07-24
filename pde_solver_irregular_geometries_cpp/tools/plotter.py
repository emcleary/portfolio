import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

def get_cmd_line_arguments():
    parser = ArgumentParser()
    parser.add_argument('filename', type=str, help='File containing data to plot')
    parser.add_argument('output', nargs='?', type=str, help='')
    return parser.parse_args()


def main():
    args = get_cmd_line_arguments()
    x = []
    y = []
    u = []
    with open(args.filename, 'r') as file:
        _, nx, ny = file.readline()[:-1].split(' ')
        nx = int(nx)
        ny = int(ny)
        for line in file.readlines():
            xi, yi, ui = line[:-1].split(' ')
            if ui == 'NAN':
                ui = np.nan
            x.append(float(xi))
            y.append(float(yi))
            u.append(float(ui))

    x = np.asarray(x).reshape(ny, nx)
    y = np.asarray(y).reshape(ny, nx)
    u = np.asarray(u).reshape(ny, nx)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    cs = ax.contourf(x, y, u)
    cbar = fig.colorbar(cs)
    if args.output:
        plt.savefig(args.output, bbox_inches='tight')
    else:
        plt.show()


if __name__=='__main__':
    main()
