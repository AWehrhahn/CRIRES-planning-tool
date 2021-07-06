from astropy.time import Time
from astropy import units as u

from crires_planning_tool.transit_planner import transit_calculation
from crires_planning_tool.interactive_graph import create_interactive_graph

planet = "TRAPPIST-1"
date_start = Time.now()
date_end = date_start + 1000 * u.day

df = transit_calculation(planet, date_start, date_end)
df.to_csv("test.csv")

create_interactive_graph(df, "plot1.html", "plot2.html")

pass

