import numpy as np
# CHECK THIS PART!!!

def transmission_function(sigma, vz): #l 473
    t = np.exp(1j * sigma * vz)
    return t

def aberration_function(sigma, c3, c5, defocus, kxgrid, kygrid): #l. 597
    pass

def aperature_function(kxgrid, kygrid, k_max): #l. 608
    pass
#   bfdiscsize = 2*int(k_max/dkx)+8 l.423
# WRITES A DIFFERENT bfdiscsiye in the timstem output   write(outfile,*) "bfdiscsize: ", bfdiscsize+1



# probe positions l.530 - l.558

#build probe function

#calcuate probe positions 

