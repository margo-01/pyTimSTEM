import numpy as np

import config as config
import ft_grids_sc
import io

# Load configuration from YAML file
cfg = config.SimConfig.from_yaml("test_input.yaml")

# Read the projected potential data 
x_coords, y_coords, vz, dx, dy = io.read_vz(cfg.data_path, cfg.unit_cell_nx, cfg.unit_cell_ny)

# Make supercell Vz from the unit cell Vz data 
vz = ft_grids_sc.adjust_vz(vz, x_coords, y_coords, cfg.reference_x, cfg.reference_y)
vz_sc = ft_grids_sc.build_supercell(vz, cfg.cells_x, cfg.cells_y)

# Set up the grids


#