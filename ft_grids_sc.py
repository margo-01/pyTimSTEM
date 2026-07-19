import numpy as np

def adjust_vz(vz, x_coords, y_coords, ref_x, ref_y): #set reference point to zero, make potential positive l.281

    ix = np.argmin(np.abs(x_coords - ref_x))
    iy = np.argmin(np.abs(y_coords - ref_y))

    reference_value = vz[ix, iy]

    return -(vz - reference_value)
    

def build_supercell(vz, cells_x, cells_y): # here the adjusted vz , l.308 onwards
    periodic_part = vz[:-1, :-1]

    repeated = np.tile(periodic_part, (cells_x, cells_y))

    # Restore one periodic boundary point on each side.
    vz_sc = np.empty(
        (repeated.shape[0] + 1, repeated.shape[1] + 1),
        dtype=repeated.dtype,
    )

    vz_sc[:-1, :-1] = repeated
    vz_sc[-1, :-1] = repeated[0, :]
    vz_sc[:-1, -1] = repeated[:, 0]
    vz_sc[-1, -1] = repeated[0, 0]

    return vz_sc

def make_real_space_grid():
    pass

def make_reciprocal_space_grid():
    pass

def fft_forward():
    pass

def fft_backward():
    pass






