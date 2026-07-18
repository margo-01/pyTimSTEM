import numpy as np

def read_vz(filename, nx, ny): #does not need nx and ny, can be inferred from data but can be useful to check the input
    vz_data = []
    with open(filename, 'r') as f:
        for line in f:
            vz_data.append(line.strip().split(","))
    vz_data = np.array(vz_data, dtype=float)
    x_coords = np.unique(vz_data[:, 0])
    y_coords = np.unique(vz_data[:, 1])

    if nx != len(x_coords) or ny != len(y_coords):
        raise ValueError(f"Expected nx={nx} and ny={ny}, but got {len(x_coords)} and {len(y_coords)} from the data, check the input file.")
   
    vz_values = vz_data[:, 2].reshape((ny, nx))
    return x_coords, y_coords, vz_values



