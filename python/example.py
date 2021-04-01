from astropy import units as u
import astroplan

from Transit_List import single_transit_calculation, get_default_constraints, full_transit_calculation

# TODO: Optimization
# Planet Observability takes a lot of time
# Already parralleized it, but still not great

date = "" # today
max_delta_days = 100
name = "WASP-107 b"
constraints = get_default_constraints()

# planet = single_transit_calculation(date, max_delta_days, name, constraints)

eclipses_list = full_transit_calculation(date, max_delta_days, constraints)

pass
