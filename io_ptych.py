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
        np.savetxt(f, frame, fmt="%.15e")



def write_4d_data_hdf5(filepath, data):
    pass
