from astropy.time import Time
from astropy import units as u

from crires_planning_tool.transit_planner import transit_calculation

planet = ["WASP-107 b", "TRAPPIST-1 b"]
date_start = Time.now()
date_end = date_start + 1000 * u.day

df = transit_calculation(planet, date_start, date_end)
df.to_csv("test.csv")

pass

