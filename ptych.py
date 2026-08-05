import numpy as np
from scipy.integrate import simpson

import ft_grids_sc as ftg
# CHECK THIS PART!!!

def transmission_function(sigma, vz): #l 473
    t = np.exp(1j * sigma * vz)
    return t

def aberration_function(wvln, c3, c5, defocus, kr): #l. 597 #for loop, allocating the array with np.zeros() ?
    expchai_f = np.exp(-1j* np.pi* wvln* kr**2 * 
                ( 0.5 * c3 * wvln**2 * kr**2   - defocus
        + (1.0 / 3.0) * c5 * wvln**4 * kr**4 )
)
    return expchai_f
    

def aperature_function(nx, ny, kr, k_max): #l. 608
    aperature = np.zeros((nx, ny))
    aperature[kr <= k_max] = 1.0
    return aperature

#   
# WRITES A DIFFERENT bfdiscsiye in the timstem output   write(outfile,*) "bfdiscsize: ", bfdiscsize+1
# probe positions l.530 - l.558

def calculate_probe_positions(scan_start_x, scan_end_x, scan_start_y, scan_end_y, probe_step_x, probe_step_y):
    np_x = int((scan_end_x - scan_start_x) / probe_step_x) + 1
    np_y = int((scan_end_y - scan_start_y) / probe_step_y) + 1

    if abs((scan_end_x-scan_start_x)-np.real(np_x-1)*probe_step_x)> abs((scan_end_x-scan_start_x)-np.real(np_x)*probe_step_x):
        np_x=np_x+1
    if abs((scan_end_y-scan_start_y)-np.real(np_y-1)*probe_step_y)> abs((scan_end_y-scan_start_y)-np.real(np_y)*probe_step_y):
        np_y=np_y+1


    num_probe_positions = np_x * np_y

    rp_x = scan_start_x + np.arange(np_x) * probe_step_x #does this have the fi
    rp_y = scan_start_y + np.arange(np_y) * probe_step_y

           
    return rp_x, rp_y, num_probe_positions


    


#build probe function for a specific probe position 

def build_rec_probe_function(kxgrid, kygrid, aperture_f, expchai_f, r_px, r_py):  #l.689
    translation = np.exp(-1j * 2.0 * np.pi * (kxgrid[:, None] * r_px + kygrid[None, :] * r_py))
    rec_probe = aperture_f * expchai_f * translation
    return rec_probe

def find_norm_bounds(scan_start_x, scan_end_x, scan_start_y, scan_end_y, xgrid, ygrid): 
    ix_start, iy_start = ftg.nearest_grid_index(
        xgrid,
        ygrid,
        scan_start_x,
        scan_start_y,
    )

    ix_end, iy_end = ftg.nearest_grid_index(
        xgrid,
        ygrid,
        scan_end_x,
        scan_end_y,
    )

    return ix_start, ix_end, iy_start, iy_end

#probe normalisation only within the scan grid (?) l.697 - 713
def normalize_probe(probe, xgrid, ygrid, ix_start, ix_end, iy_start, iy_end):

    probe_sqr = np.abs(
        probe[ix_start : ix_end + 1, iy_start : iy_end + 1]
    ) ** 2

    integrated_intensity = simpson(
        simpson(
            probe_sqr,
            x=ygrid[iy_start : iy_end + 1],
            axis=1,
        ),
        x=xgrid[ix_start : ix_end + 1],
    )
        
    

    if integrated_intensity <= 0:
        raise ValueError(
            "Integrated probe intensity is not positive"
        )

    probe = probe / np.sqrt(integrated_intensity)

    return probe

def centered_crop(psi_f, detector_size): #l 728
   # if detector_size > min(psi_f.shape):
     #   raise ValueError("Detector crop is larger than the diffraction grid")
    start_x = (psi_f.shape[0] - detector_size) // 2
    start_y = (psi_f.shape[1] - detector_size) // 2
    crop = psi_f[start_x:start_x + detector_size, start_y:start_y + detector_size]
    crop = np.abs(crop) ** 2
    return crop


















