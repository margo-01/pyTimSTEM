import numpy as np
# CHECK THIS PART!!!

ELEMENTARY_CHARGE_SI = 1.602176565e-19
ELECTRON_MASS_SI = 9.10938291e-31
SPEED_OF_LIGHT_SI = 299792458.0


def electron_parameters(energy_kev: float):
    energy_atomic = energy_kev / 0.027211385 # convert energy from keV to atomic units l.234

    wavelength = 860.7963871 / np.sqrt(
        energy_atomic * (energy_atomic + 37538.0) 
    ) #calculate relativistic wavelength l. 237

    gamma = (
        1.0
        + energy_kev
        * 1000.0
        * ELEMENTARY_CHARGE_SI
        / (ELECTRON_MASS_SI * SPEED_OF_LIGHT_SI**2)
    )

    sigma = gamma * wavelength / 6.283185307

    return wavelength, sigma


def transmission_function(sigma, vz): #l 473
    t = np.exp(1j * sigma * vz)
    return t


# probe positions l.530 - l.558

#build probe function

#calcuate probe positions 
