################################################################################
# This file is part of CosTuuM
# Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CosTuuM is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
###############################################################################

import numpy as np
import CosTuuM
import scipy.interpolate as interpol

## required input parameters ##################################################
dust_type = CosTuuM.SILICON
dust_size = 1.0e-6  # m
radiation_wavelength = 1.0e-5  # m
theta_in = 0.5 * np.pi  # radians
theta_out = 0.5 * np.pi  # radians
phi = np.pi  # radians - relative angle between incoming and outgoing ray
dust_shape = CosTuuM.SingleShapeShapeDistribution(2.0)  # oblate grain
# alignment distribution for the dust grains
# currently, CosTuuM only has one type of distribution that assumes no
# alignment below some size and a named type of alignment above that size
# we choose Mishchenko alignment above 0.1 um. Since we only have one grain
# size, the grains will have this alignment.
dust_alignment = CosTuuM.SizeBasedAlignmentDistribution(
    1.0e-7, CosTuuM.MISHCHENKO_ALIGNMENT
)
# Use builtin Draine refractive indices for an assumed dust temperature of 20K
dust_properties = CosTuuM.DraineDustProperties(dust_temperature=20.0)

## optional input parameters ##################################################
# All parameters are set to their default values
# Minimum and maximum allowed order of spherical basis functions used in the
# expansion
minimum_order = 10
maximum_order = 100
# Gauss-Legendre factor. For spherical basis functions of a given order N, the
# number of Gauss-Legendre quadrature points used for numerical integration is
# set to this factor times N
gauss_legendre_factor = 2
# Relative tolerance value used to decide when the T-matrix calculation is
# converged
tolerance = 1.0e-4
# Maximum allowed memory usage of the algorithm, in bytes. Set this to a value
# lower than the available memory on the machine if you don't want to run out
# of memory. The algorithm will try to reduce its memory footprint to fit,
# although there is no guarantee that this will work for large expansion orders
maximum_memory_size = 10000000000  # ~10GB
# Number of shared-memory parallel threads to use in the calculation.
# The recommended value is a number close to the number of available cores on
# your machine.
number_of_threads = 4
# Optional log files that can be used to assess the performance of the
# algorithm
quicksched_graph_log = None
quicksched_task_log = None
quicksched_task_type_log = None
memory_log = None
# Make the algorithm more verbose than it already is?
verbose = False

result = CosTuuM.get_scattering_matrix_table(
    types=dust_type,
    sizes=dust_size,
    wavelengths=radiation_wavelength,
    theta_in=theta_in,
    theta_out=theta_out,
    phi=phi,
    shape_distribution=dust_shape,
    alignment_distribution=dust_alignment,
    dust_properties=dust_properties,
    minimum_order=minimum_order,
    maximum_order=maximum_order,
    gauss_legendre_factor=gauss_legendre_factor,
    tolerance=tolerance,
    maximum_memory_size=maximum_memory_size,
    number_of_threads=number_of_threads,
    quicksched_graph_log=quicksched_graph_log,
    quicksched_task_log=quicksched_task_log,
    quicksched_task_type_log=quicksched_task_type_log,
    memory_log=memory_log,
    verbose=verbose,
)

print(result)
