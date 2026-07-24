import numpy as np
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

#build probe function

#calcuate probe positions 

