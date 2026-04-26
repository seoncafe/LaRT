#!/bin/bash

python make_amr_sphere_radial.py --level_min 2 --level_max 7 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 7 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 8 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 9 --refine_boundary --boundary_level_max 5 --plot --compare \
				-o 'sphere_amr1.fits.gz'
python make_amr_sphere_radial.py --level_min 2 --level_max 10 --refine_boundary --boundary_level_max 5 --plot --compare \
				-o 'sphere_amr2.fits.gz'
python make_amr_sphere_radial.py --level_min 3 --level_max 10 --refine_boundary --boundary_level_max 5 --plot --compare \
				-o 'sphere_amr3.fits.gz'
python make_amr_sphere_radial.py --level_min 5 --level_max 10 --refine_boundary --boundary_level_max 6 \
				--velocity rotating_galaxy_halo --vrot 100 --rinner 0.1 \
				-o 'sphere_Vrot100_amr.fits.gz' --plot
python make_amr_sphere_radial.py --level_min 5 --level_max 10 --refine_boundary --boundary_level_max 6 \
				--velocity rotating_galaxy_halo --vrot 300 --rinner 0.1 \
				-o 'sphere_Vrot300_amr.fits.gz' --plot
