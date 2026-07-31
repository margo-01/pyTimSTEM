import numpy as np

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
chai = ptych.aberration_function(cfg.wavelength, cfg.c3_au, cfg.c5_au, cfg.defocus_au, cfg.kr)
aperture = ptych.aperature_function(cfg.unit_cell_nx, cfg.unit_cell_ny, cfg.kr, cfg.k_max)
rp_x, rp_y = ptych.calculate_probe_positions(cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, cfg.probe_step_x, cfg.probe_step_y)

bfdiscsize = 2*int(cfg.k_max/dkx)+8 #l.423
detector_size = bfdiscsize + 1


#loop over probe positions


for ix, probe_x in enumerate(rp_x):
    for iy, probe_y in enumerate(rp_y):
        pattern = np.zeros((bfdiscsize, bfdiscsize))
        probe_f = ptych.build_rec_probe_function(kxgrid, kygrid, aperture, chai, rp_x, rp_y, iy)
        probe = ft_grids_sc.fft_backward(probe_f, dx, dy)
        probe = ptych.normalize_probe(probe, x, y, cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, dx, dy)

        exit_wave = probe * t
        diffraction_wave = ft_grids_sc.fft_forward(exit_wave, dx, dy)
        diffraction = np.abs(diffraction_wave) ** 2

        pattern[ix, iy] = ptych.centered_crop(diffraction, detector_size)
        io_ptych.write_4d_data()







