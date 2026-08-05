import numpy as np
from pathlib import Path

import config as config
import ft_grids_sc
import io_ptych
import ptych 

# Load configuration from YAML file
cfg = config.SimConfig.from_yaml("test_input.yaml")
cfg.convert_units()

# Read the projected potential data 
x_coords, y_coords, vz, dx, dy = io_ptych.read_vz(
    cfg.data_path,
    cfg.potential_filename,
    cfg.unit_cell_nx,
    cfg.unit_cell_ny,
)
# Make supercell Vz from the unit cell Vz data 
vz = ft_grids_sc.adjust_vz(vz, x_coords, y_coords, cfg.reference_x, cfg.reference_y)
vz_sc = ft_grids_sc.build_supercell(vz, cfg.cells_x, cfg.cells_y)

# Set up the grids

nx, ny = vz_sc.shape
x, y, Lx, Ly = ft_grids_sc.make_real_space_grids(nx, ny, dx, dy, cfg.cells_x, cfg.cells_y)

kxgrid, kygrid, dkx, dky = ft_grids_sc.make_reciprocal_space_grid(
        nx,
        ny,
        Lx,
        Ly,
    )

kr = ft_grids_sc.make_recip_radial_grid(
    kxgrid,
    kygrid,
)

# Ptych setup
t = ptych.transmission_function(cfg.sigma, vz_sc)
chai = ptych.aberration_function(cfg.wavelength, cfg.c3_au, cfg.c5_au, cfg.defocus_au, kr)
aperture = ptych.aperature_function(nx, ny, kr, cfg.k_max)
rp_x, rp_y, num_probe_positions = ptych.calculate_probe_positions(cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, cfg.probe_step_x, cfg.probe_step_y)

bfdiscsize = 2*int(cfg.k_max/dkx)+8 #l.423
detector_size = bfdiscsize + 1

#checking if an output file already exists in the directory

filename = Path(cfg.data_path) / "4D_000.raw"
# check if the file exists in the directory
if filename.exists():
    raise FileExistsError(f"The file {filename} already exists. Please delete the existing file.")

# loop over probe positions

for probe_x in rp_x:
    for probe_y in rp_y:
        probe_f = ptych.build_rec_probe_function(
            kxgrid, kygrid, aperture, chai, probe_x, probe_y
        )
        probe = ft_grids_sc.fft_backward(probe_f, dx, dy)
        ix_start, ix_end, iy_start, iy_end = ptych.find_norm_bounds(
            cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, x, y
        )
        probe = ptych.normalize_probe(probe, x, y, ix_start, ix_end, iy_start, iy_end)

        exit_wave = probe * t
        psi_t = ft_grids_sc.fft_forward(exit_wave, dx, dy) #in fortran, this part is calculated on the full grid 
        pattern = ptych.centered_crop(psi_t, detector_size)
        io_ptych.write_4d_data(filename, pattern)







