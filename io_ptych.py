import numpy as np
from pathlib import Path

def read_vz(filepath, filename, nx, ny): #does not need nx and ny, can be inferred from data but can be useful to check the input
    vz_data = []
    with open(Path(filepath) / filename, 'r') as f:
        for line in f:
            vz_data.append(line.strip().split(","))
    vz_data = np.array(vz_data, dtype=float)
    x_coords = np.unique(vz_data[:, 0])
    y_coords = np.unique(vz_data[:, 1]) 

    dx = x_coords[1] - x_coords[0]
    dy = y_coords[1] - y_coords[0]

    if nx != len(x_coords) or ny != len(y_coords):
        raise ValueError(f"Expected nx={nx} and ny={ny}, but got {len(x_coords)} and {len(y_coords)} from the data, check the input file.")
   
    vz_values = vz_data[:, 2].reshape((nx, ny))
    return x_coords, y_coords, vz_values, dx, dy

def write_4d_data(filename, data):

    frame = np.asarray(data)

    with open(Path(filename), "a", encoding="ascii", newline="\n") as f:
        np.savetxt(f, frame.ravel(order="F"), fmt="%.15e")



def write_4d_data_hdf5(filepath, data):
    pass


# Debug output writers corresponding to timstemv1.6_debug.f90.
def _debug_file(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="ascii", newline="\n")


def write_reciprocal_grid(path, kx, ky, kr):
    with _debug_file(path) as f:
        f.write("# i j kx ky kr\n")
        for i in range(len(kx)):
            for j in range(len(ky)):
                f.write(f"{i+1:8d} {j+1:8d} {kx[i]:24.16E} {ky[j]:24.16E} {kr[i, j]:24.16E}\n")


def write_radial_grid(path, grid, name):
    with _debug_file(path) as f:
        f.write(f"# index {name}\n")
        for i, value in enumerate(grid, 1):
            f.write(f"{i:12d} {value:24.16E}\n")


def write_probe_positions(path, probe_x, probe_y):
    with _debug_file(path) as f:
        f.write("# probe_number x y\n")
        n = 1
        for x in probe_x:
            for y in probe_y:
                f.write(f"{n:12d} {x:24.16E} {y:24.16E}\n")
                n += 1


def write_complex_grid(path, array):
    with _debug_file(path) as f:
        f.write("# i j real imag intensity\n")
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                value = array[i, j]
                f.write(f"{i+1:8d} {j+1:8d} {value.real:24.16E} {value.imag:24.16E} {abs(value)**2:24.16E}\n")


def write_normalization(path, number, x, y, bounds, integral, factor):
    ix_start, ix_end, iy_start, iy_end = bounds
    with _debug_file(path) as f:
        f.write("# probe_number x y ix_start ix_end iy_start iy_end integral norm_factor\n")
        f.write(f"{number:12d} {x:24.16E} {y:24.16E} "
                f"{ix_start+1:8d} {ix_end+1:8d} {iy_start+1:8d} {iy_end+1:8d} "
                f"{integral:24.16E} {factor:24.16E}\n")


def write_detector(path, intensity):
    with _debug_file(path) as f:
        f.write("# i j intensity\n")
        for i in range(intensity.shape[0]):
            for j in range(intensity.shape[1]):
                f.write(f"{i+1:8d} {j+1:8d} {intensity[i, j]:24.16E}\n")


def write_xyz(path, x, y, array, fmt=".10f"):
    with _debug_file(path) as f:
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                f.write(f"{x[i]:{fmt}} , {y[j]:{fmt}} , {array[i, j]:{fmt}}\n")


def write_aperture(path, x, y, aperture):
    # Fortran fmt3 = (F10.6,A2,F10.6,A2,F15.11)
    with _debug_file(path) as f:
        for i in range(aperture.shape[0]):
            for j in range(aperture.shape[1]):
                f.write(f"{x[i]:10.6f}, {y[j]:10.6f}, {aperture[i, j]:15.11f}\n")
