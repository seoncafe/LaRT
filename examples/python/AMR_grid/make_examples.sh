#!/bin/bash

python make_amr_sphere_radial.py --level_min 2 --level_max 7 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 7 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 8 --refine_boundary --boundary_level_max 5 --plot --compare
python make_amr_sphere_radial.py --level_min 2 --level_max 9 --refine_boundary --boundary_level_max 5 --plot --compare -o 'sphere_amr1.fits.gz'
python make_amr_sphere_radial.py --level_min 2 --level_max 10 --refine_boundary --boundary_level_max 5 --plot --compare -o 'sphere_amr2.fits.gz'
python make_amr_sphere_radial.py --level_min 3 --level_max 10 --refine_boundary --boundary_level_max 5 --plot --compare -o 'sphere_amr3.fits.gz'
python make_amr_sphere_radial.py --level_min 3 --level_max 10 --refine_boundary --boundary_level_max 6 --plot --compare -o 'sphere_amr4.fits.gz'
