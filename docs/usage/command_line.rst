Command Line Interface
======================

The CRIRES Planning Tool features a command line Interface to access the functionality.
The recommended way to call it is `python -m crires-planning-tool`, 
but it can also be accessed by calling `python transit_planner.py` directly.

usage: crires-planning-tool [-h] [-O OBSERVER] [-C CATALOG] [-m MODE]
                            [-o OUTPUT] [-f] [-p] [--plot-file-1 PLOT_FILE_1]
                            [--plot-file-2 PLOT_FILE_2] [-v] [-s]
                            begin end [planets [planets ...]]

CRIRES+ planning tool to determine the observable transits of one or more
planets

positional arguments:
  begin                 first date to check for transits
  end                   last date to check for transits
  planets               Names of the planets to calculate in the format 'star
                        letter'. If just the name of the star is given it will
                        get all the planets of that system. If mode is set to
                        'criteria' this defines the criteria instead of the
                        planets, look at the NASA Exoplanet Archive
                        documentation of the TAP interface for available
                        criteria.

optional arguments:
  -h, --help            show this help message and exit
  -O OBSERVER, --observer OBSERVER
                        Location of the observer, by default 'paranal'. This
                        can be any name supported by astropy.
  -C CATALOG, --catalog CATALOG
                        Name of the data catalog to use, by default 'nexa'
                        (Nasa Exoplanet Archive)
  -m MODE, --mode MODE  Which mode to use 'planets' or 'criteria'. In mode
                        'planets' we define the specific planets to
                        investigate, in mode 'criteria' we specify conditions
                        that the stars and planets need to fullfill.
  -o OUTPUT, --output OUTPUT
                        Save the data to this file
  -f, --file            Load planet names from a file instead, one planet name
                        per line. Lines starting with # are considered
                        comments
  -p, --plot            If set will create an interactive plot
  --plot-file-1 PLOT_FILE_1
                        filename for the first plot
  --plot-file-2 PLOT_FILE_2
                        filename for the second plot
  -v, --verbose         verbosity level up to -vv
  -s, --silent          removes all console output, except for errors
