import numpy as np

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

#   bfdiscsize = 2*int(k_max/dkx)+8 l.423
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

    rp_x = scan_start_x + np.arange(np_x) * probe_step_x
    rp_y = scan_start_y + np.arange(np_y) * probe_step_y

           
    return rp_x, rp_y, num_probe_positions


    


#build probe function for a specific probe position 

def build_rec_probe_function(kxgrid, kygrid, dkx, dky, aperture_f, expchai_f, r_px, r_py, nrp):  #l.689
    translation = np.exp(-1j * 2.0 * np.pi * (kxgrid[:, None] * r_px[nrp] + kygrid[None, :] * r_py[nrp]))
    rec_probe = aperture_f * expchai_f * translation
    return rec_probe

