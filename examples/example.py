from astroplan import download_IERS_A
import crires_planning_tool.classes_methods.Helper_fun as fun
from crires_planning_tool.Transit_List import (etc_calculator, full_transit_calculation,
                          get_default_constraints, single_transit_calculation)
from crires_planning_tool.classes_methods.obsob import estimate_snr
# download_IERS_A()

# TODO: Optimization
# ETC calls take the majority of the time
# Either reduce number of calls significantly (interpolate inbetween?)
# Or somehow send a batch request that should be faster in theory.


date = "" # today
max_delta_days = 1000
name = "WASP-107 b"
constraints = get_default_constraints()

# For a single system
# With the SNR estimate
# Only a single exposure, so snr is just snr_median
eclipse = single_transit_calculation(date, max_delta_days, name, constraints)
eclipse = estimate_snr([eclipse], minimum_snr=100)[0]
snr = [obs["snr_median"] for obs in eclipse.eclipse_observable]
eclipses_list = [eclipse]

# # With the ETC calculation
# # Note that the snr for the whole transit is snr_median * n_exposures_possible
# eclipse_etc = single_transit_calculation(date, max_delta_days, name, constraints)
# eclipse_etc = etc_calculator([eclipse_etc])[0]
# snr_etc = [obs["snr_median"] * obs["n_exposures_possible"] for obs in eclipse_etc.eclipse_observable]

# For all systems
# eclipses_list = full_transit_calculation(date, max_delta_days, constraints)
# eclipses_list = estimate_snr(eclipses_list, minimum_snr=100)
# eclipses_list = etc_calculator(eclipses_list)
# fun.save_pickled("eclipses_list.pkl", eclipses_list)

data = fun.to_pandas(eclipses_list)
fun.save_csv("{path}/test.csv", data)
fun.save_excel("{path}/test.xlsx", data)

