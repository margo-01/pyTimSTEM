from pathlib import Path
import time
import numpy as np
from scipy.integrate import simpson

import config as config
import ft_grids_sc
import io_ptych
import ptych 
from utils import ProgressLogger

# Load configuration from YAML file
cfg = config.SimConfig.from_yaml("test_input.yaml")

output_dir = Path(cfg.data_path) / "pyTimSTEM_out"
output_dir.mkdir(parents=True, exist_ok=True)

cfg.output_dir = output_dir

progress_logger = ProgressLogger(cfg)
progress_logger.write("Welcome, this is pyTimSTEM! Getting ready to run your simulation.")
cfg.convert_units()
progress_logger.write(f"I converted the microscope parameter units to atomic units, c3: {cfg.c3_au}, c5: {cfg.c5_au}, defocus: {cfg.defocus_au}, alpha_max: {cfg.alpha_max_rad}, k_max: {cfg.k_max}, sigma: {cfg.sigma}")
# Read the projected potential data 
progress_logger.write(f"I read the projected potential data from the {cfg.potential_filename} file in the {cfg.data_path} directory")

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
progress_logger.write(f"I am setting up the grids now")

nx, ny = vz_sc.shape
x, y, Lx, Ly = ft_grids_sc.make_real_space_grids(nx, ny, dx, dy)

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
progress_logger.write(f"The real space grid has {nx} x {ny} points, and the reciprocal space grid has {kxgrid.shape[0]} x {kygrid.shape[0]} y points, with dkx: {dkx}, dky: {dky}")

# Grid diagnostics, matching timstemv1.6_debug.f90.
io_ptych.write_xyz(output_dir / "post_uc_vz.out", x_coords, y_coords, vz)
io_ptych.write_xyz(output_dir / "supercell_vz.timstem", x, y, vz_sc)
io_ptych.write_reciprocal_grid(output_dir / "debug_reciprocal_grid.dat", kxgrid, kygrid, kr)

## ALL CORRECT SO FAR CHECKED AGAINST FORTRAN ##

dr = 0.25 * min(dx, dy)
dkr = 0.25 * min(dkx, dky)
rmax = np.sqrt(Lx**2 + Ly**2)
nrmax = int(rmax / dr)
nkmax = int(cfg.k_max / dkr)
rad_grid = np.arange(nrmax) * (rmax / (nrmax - 1))
rad_gridk = np.arange(nkmax) * (cfg.k_max / (nkmax - 1))

# Ptych setup
t = ptych.transmission_function(cfg.sigma, vz_sc)
chai = ptych.aberration_function(cfg.wavelength, cfg.c3_au, cfg.c5_au, cfg.defocus_au, kr)
aperture = ptych.aperature_function(nx, ny, kr, cfg.k_max)
io_ptych.write_aperture(
    output_dir / "aperturestem.timstem",
    x_coords,
    y_coords,
    aperture[:len(x_coords), :len(y_coords)],
)

## APERTURE CHECKED AGAINST FORTRAN, ALL CORRECT ##
rp_x, rp_y, num_probe_positions = ptych.calculate_probe_positions(cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, cfg.probe_step_x, cfg.probe_step_y)
io_ptych.write_probe_positions(output_dir / "debug_probe_positions.dat", rp_x, rp_y)

## PROBE POSITIONS FINE ##

bfdiscsize = 2*int(cfg.k_max/dkx)+8 #l.423
detector_size = bfdiscsize + 1

progress_logger.write(f"I am using a bfdiscsize of {bfdiscsize+1}.")


## BFDISCSIZE CHECKED AGAINST FORTRAN, CORRECT ##

#checking if an output file already exists in the directory



filename = output_dir / "4D_000.timstem"# check if the file exists in the directory

if filename.exists():
    raise FileExistsError(f"The file {filename} already exists. Please delete the existing file.")

# loop over probe positions
progress_logger.write(f"We're all set! Starting the loop over probe positions. Total number of probe positions: {num_probe_positions}")
loop_start = time.perf_counter()

for ix_probe, probe_x in enumerate(rp_x):
    for iy_probe, probe_y in enumerate(rp_y):
        probe_number = ix_probe * len(rp_y) + iy_probe + 1
        probe_f = ptych.build_rec_probe_function(
            kxgrid, kygrid, aperture, chai, probe_x, probe_y
        )
        probe = ft_grids_sc.fft_backward(probe_f, dx, dy)
        if probe_number == 1:
            io_ptych.write_complex_grid(output_dir / "debug_probe_fourier.dat", probe_f)
            io_ptych.write_complex_grid(output_dir / "debug_probe_real.dat", probe)
## FT AND REAL PROBES SEEM CORRECT ##

        ix_start, ix_end, iy_start, iy_end = ptych.find_norm_bounds(
            cfg.scan_start_x, cfg.scan_end_x, cfg.scan_start_y, cfg.scan_end_y, x, y
        )
        probe_sqr = abs(probe[ix_start:ix_end + 1, iy_start:iy_end + 1]) ** 2
        integral = simpson(simpson(probe_sqr, x=y[iy_start:iy_end + 1], axis=1),
                           x=x[ix_start:ix_end + 1])
        norm_factor = 1.0 / np.sqrt(integral)
        probe = ptych.normalize_probe(probe, x, y, ix_start, ix_end, iy_start, iy_end)
        if probe_number == 1:
            io_ptych.write_normalization(
                output_dir / "debug_probe_normalization.dat", probe_number,
                probe_x, probe_y, (ix_start, ix_end, iy_start, iy_end),
                integral, norm_factor,
            )
## NORMALIZATION WORKS BUT NUMERICAL DIFFERENCES WITH FORTRAN ##

        exit_wave = probe * t
        psi_t = ft_grids_sc.fft_forward(exit_wave, dx, dy) #in fortran, this part is calculated on the full grid 
        pattern = ptych.centered_crop(psi_t, detector_size)
        if probe_number == 1:
            io_ptych.write_complex_grid(output_dir / "debug_exit_fourier.dat", psi_t)
            ## FT OF EXIT WAVE SEEMS CORRECT BUT NUMERICAL DIFFERENCES ##
            io_ptych.write_detector(output_dir / "debug_detector.dat", pattern)
## PATTERN IS CORRECT!! ##
        io_ptych.write_4d_data(filename, pattern)
loop_elapsed = time.perf_counter() - loop_start
progress_logger.write(f"Loop over probe positions completed in {loop_elapsed:.2f} seconds.")

progress_logger.write(f"All done! I saved the 4D dataset to {filename}. Please check the output and see if it is reasonable")    






