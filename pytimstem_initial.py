import numpy as np

import config as config
import ft_grids_sc
import io_ptych
import ptych 

# Load configuration from YAML file
cfg = config.SimConfig.from_yaml("test_input.yaml")

# Read the projected potential data 
x_coords, y_coords, vz, dx, dy = io_ptych.read_vz(cfg.data_path, cfg.unit_cell_nx, cfg.unit_cell_ny)

# Make supercell Vz from the unit cell Vz data 
vz = ft_grids_sc.adjust_vz(vz, x_coords, y_coords, cfg.reference_x, cfg.reference_y)
vz_sc = ft_grids_sc.build_supercell(vz, cfg.cells_x, cfg.cells_y)

# Set up the grids

x, y, Lx, Ly = ft_grids_sc.make_real_space_grids(cfg.unit_cell_nx, cfg.unit_cell_ny, dx, dy, cfg.cells_x, cfg.cells_y)
kxgrid, kygrid, kr, dkx, dky = ft_grids_sc.make_reciprocal_space_grid(cfg.unit_cell_nx, cfg.unit_cell_ny, dx, dy)


# Ptych setup
t = ptych.transmission_function(cfg.sigma, vz_sc)
chai = ptych.aberration_function(cfg.wavelength, cfg.c3, cfg.c5, cfg.defocus, cfg.kr)
aperture = ptych.aperature_function(cfg.unit_cell_nx, cfg.unit_cell_ny, cfg.kr, cfg.k_max)
rp_x, rp_y = ptych.make_probe_positions(cfg.probe_positions, cfg.unit_cell_nx, cfg.unit_cell_ny, cfg.cells_x, cfg.cells_y, dx, dy)

#loop over probe positions


for ix, probe_x in enumerate(rp_x):
    for iy, probe_y in enumerate(rp_y):
        probe_f = ptych.build_probe_fourier(...)
        probe = ptych.fft_backward(probe_f, dx, dy)
        probe = ptych.normalize_probe(probe, ...)

        exit_wave = probe * t
        diffraction_wave = ft_grids_sc.fft_forward(exit_wave, dx, dy)
        diffraction = np.abs(diffraction_wave) ** 2

        patterns[ix, iy] = centered_crop(diffraction, detector_size)






