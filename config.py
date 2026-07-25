from dataclasses import dataclass

import omegaconf
import yaml
import numpy as np

@dataclass
class SimConfig:
    data_path: str
    potential_filename: str

    unit_cell_nx: int #number of grid points in the unit cell along x
    unit_cell_ny: int #number of grid points in the unit cell along y

    cells_x: int # super cell size in x direction (number of unit cells)
    cells_y: int # super cell size in y direction (number of unit cells) 
    # supercell needs to be approximately square to avoid artifacts in the diffraction pattern

    energy_kev: float

    scan_start_x: float # starting x position for the scan in bohr 
    scan_end_x: float # ending x position for the scan in bohr
    scan_start_y: float # starting y position for the scan in bohr
    scan_end_y: float # ending y position for the scan in bohr

    probe_step_x: float # step size for the probe in x direction in bohr
    probe_step_y: float # step size for the probe in y direction in bohr

    c3: float
    c5: float
    defocus: float
    alpha_max_mrad: float # maximum convergence angle in mrad

    reference_x: float = 0.0 #reference poiint for the projected potential in bohr
    reference_y: float = 0.0 


    @classmethod
    def from_yaml(cls, yaml_file):
        cfg_file = omegaconf.OmegaConf.load(yaml_file)
        default_config = omegaconf.OmegaConf.structured(cls)
        merged_config = omegaconf.OmegaConf.merge(default_config, cfg_file)
        return omegaconf.OmegaConf.to_object(merged_config)

    def convert_units(self):
        ELEMENTARY_CHARGE_SI = 1.602176565e-19
        ELECTRON_MASS_SI = 9.10938291e-31
        SPEED_OF_LIGHT_SI = 299792458.0

        energy_atomic = self.energy_kev / 0.027211385 # convert energy from keV to atomic units l.234

        wavelength = 860.7963871 / np.sqrt(
            energy_atomic * (energy_atomic + 37538.0) 
        ) #calculate relativistic wavelength l. 237

        gamma = (
            1.0
            + self.energy_kev
            * 1000.0
            * ELEMENTARY_CHARGE_SI
            / (ELECTRON_MASS_SI * SPEED_OF_LIGHT_SI**2)
        )

        self.wavelength = wavelength
        self.gamma = gamma
        self.sigma = gamma * wavelength / 6.283185307

        self.c3 *= 18897259.890 # convert c3 from mm to atomic units
        self.c5 *= 18897259.890 # convert c5 from mm to atomic units
        self.defocus *= 1.889725989 # convert defocus from nm to atomic units
        self.alpha_max_mrad *= 0.001 # convert alpha_max from mrad to rad

        self.k_max=self.alpha_max_mrad/self.wavelength


