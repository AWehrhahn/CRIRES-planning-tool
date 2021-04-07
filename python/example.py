from astropy import units as u
import astroplan

import classes_methods.Helper_fun as fun
from Transit_List import single_transit_calculation, get_default_constraints, full_transit_calculation, etc_calculator

# TODO: Optimization
# ETC calls take the majority of the time
# Either reduce number of calls significantly (interpolate inbetween?)
# Or somehow send a batch request that should be faster in theory.

date = "" # today
max_delta_days = 100
name = "WASP-107 b"
constraints = get_default_constraints()

eclipse = single_transit_calculation(date, max_delta_days, name, constraints)

eclipses_list = full_transit_calculation(date, max_delta_days, constraints)
eclipses_list = etc_calculator(eclipses_list)
fun.pickle_dumper_objects("eclipses_list.pkl", eclipses_list)

pass
