from dataclasses import dataclass

import omegaconf
import yaml

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