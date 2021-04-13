from astroplan import download_IERS_A
import classes_methods.Helper_fun as fun
from Transit_List import (etc_calculator, full_transit_calculation,
                          get_default_constraints, single_transit_calculation)
from classes_methods.obsob import estimate_snr
# download_IERS_A()

# TODO: Optimization
# ETC calls take the majority of the time
# Either reduce number of calls significantly (interpolate inbetween?)
# Or somehow send a batch request that should be faster in theory.


date = "" # today
max_delta_days = 100
name = "WASP-107 b"
constraints = get_default_constraints()

# For a single system
eclipse = single_transit_calculation(date, max_delta_days, name, constraints)

# For all systems
eclipses_list = full_transit_calculation(date, max_delta_days, constraints)
eclipses_list = estimate_snr(eclipses_list, minimum_snr=100)
# eclipses_list = etc_calculator(eclipses_list)
fun.save_pickled("eclipses_list.pkl", eclipses_list)

pass
