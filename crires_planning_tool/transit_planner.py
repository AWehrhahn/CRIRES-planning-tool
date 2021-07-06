#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Does all the transit planning

@author: Ansgar Wehrhahn
"""

import argparse
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Value
from os.path import dirname, join

import astroplan as ap
import dateparser
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from tqdm import tqdm

try:
    from .interactive_graph import create_interactive_graph
    from .nasa_exoplanets_archive import NasaExoplanetsArchive
except ImportError:
    from interactive_graph import create_interactive_graph
    from nasa_exoplanets_archive import NasaExoplanetsArchive

logger = logging.getLogger(__name__)


def get_default_constraints():
    """
    Defines the default constraints if none are explicitly set.
    Those are Altitude above 30 deg, airmass less than 1.7, 
    astronomical twilight at the observing location, and at 
    least 45 deg seperation from the moon

    Returns
    -------
    constraints : list
        list of constraints
    """    
    # Altitude constraints definition
    altcons = ap.AltitudeConstraint(min=+30 * u.deg, max=None)

    # Airmass constraints definition
    airmasscons = ap.AirmassConstraint(min=None, max=1.7)

    # Astronomical Nighttime constraints definition: 
    # begin and end of each night at paranal as 
    # AtNightConstraint.twilight_astronomical
    night_cons_per_night = ap.AtNightConstraint.twilight_astronomical()

    # Moon Constraint
    mooncons = ap.MoonSeparationConstraint(min=+45 * u.deg, max=None)

    constraints = [night_cons_per_night, altcons, airmasscons, mooncons]
    return constraints


def calculate_airmass(obstime, observer, target):
    """
    Get the airmass at obstime

    Parameters
    ----------
    obstime : Time
        observation times
    observer : astroplan.Observer
        telescope location
    target : astroplan.FixedTarget
        Star location

    Returns
    -------
    airmass : Quantity
        airmass at obstime
    """
    altaz = observer.altaz(obstime, target)
    airmass = altaz.secz
    return airmass


@u.quantity_input
def estimate_snr(
    magnitude: u.mag, exptime: u.second, airmass: u.one, method="crires", mband="K"
):
    """
    Computes SNR of stellar spectrum from magnitude, exposure time (in seconds) and airmass
    @author fabian lesjak

    Parameters
    ----------
    magnitude : Quantity
        stellar magnitude in the band specified by mband
    exptime : Quantity
        exposure time
    airmass : Quantity
        airmass
    method : {"crires", "carmenes"}, optional
        which formula to use, by default "crires"
    mband : str, optional
        which spectral band to use, by default "K"

    Returns
    -------
    snr : Quantity
        Signal to Noise Ratio
    """
    # SNR with airmass = 1
    mag = magnitude.to_value(u.mag)
    et = exptime.to_value(u.second)
    am = airmass.to_value(u.one)

    # This is from old Carmenes documentation, factor 1.1774 so that it agrees better with Ansgars result
    if method == "carmenes":
        if mband == "J":
            snr = 1.1774 * 100 / np.sqrt(40 * 10 ** ((mag - 4.2) / 2.5)) * np.sqrt(et)
        else:
            raise ValueError("Use Jmag for calculation of Carmenes SNR.")
    elif method == "crires":
        if mband == "K":
            snr = 449.4241 * np.sqrt(10 ** (-mag / 2.5)) * np.sqrt(et) - 6.3144
        else:
            raise ValueError("Use Kmag for calculation of Crires SNR.")
    else:
        raise ValueError("Method not recognized. Use Crires or Carmenes.")

    # Scale to airmass = airm
    extcof = 0.05  # extinction coefficient, see Lombardi et al., 2018
    snr *= 10 ** (extcof / 5 * (1 - airmass))

    return snr


def single_planet(data, date_start, date_end, constraints, observer, verbose=0):
    """
    Performs the observability and SNR estimation for a single planet

    Parameters
    ----------
    data : dict
        data for this planet
    date_start : Time
        start date
    date_end : Time
        end date
    constraints : list
        list of observability constraints
    observer : astroplan.Observer
        telescope location
    verbose : int, optional
        how much information to print, by default 0

    Returns
    -------
    result : dict
        contains everything about this planet, one entry per observable transit
    """
    primary_eclipse_time = Time(data["pl_tranmid"], format="jd")
    orbital_period = data["pl_orbper"] * u.day
    eclipse_duration = data["pl_trandur"] * u.hour
    name = data["pl_name"]
    # Update the progress bar
    # Define the target coordinates
    coord = SkyCoord(ra=data["ra"], dec=data["dec"], unit=u.deg)
    target = ap.FixedTarget(coord, name=name)
    # Define the star-planet system
    system = ap.EclipsingSystem(
        primary_eclipse_time=primary_eclipse_time,
        orbital_period=orbital_period,
        duration=eclipse_duration,
        name=name,
    )
    # Find all eclipses of this system
    n_max_eclipses = (date_end - date_start) / orbital_period
    n_max_eclipses = max(n_max_eclipses, 1)
    eclipses = system.next_primary_ingress_egress_time(
        date_start, n_eclipses=n_max_eclipses
    )
    # If the last eclipse is past the final date, just cut it off
    if eclipses[-1][0] > date_end:
        eclipses = eclipses[:-1]
    if len(eclipses) == 0:
        if verbose >= 1:
            logger.warning(f"No observable transits found for planet {name}")
        return None

    # Check if they are observable
    is_observable = ap.is_event_observable(
        constraints, observer, target, times_ingress_egress=eclipses
    )[0]
    observable_eclipses = eclipses[is_observable]
    n_eclipses = len(observable_eclipses)
    if n_eclipses == 0:
        if verbose >= 0:
            logger.warning(f"No observable transits found for planet {name}")
        return None

    # Calculate the SNR for those that are left
    magnitude = data["sy_kmag"] * u.mag
    exptime = eclipse_duration
    # exptime = np.diff(observable_eclipses.jd, axis=1).ravel() * u.day
    obstime = Time(np.mean(observable_eclipses.jd, axis=1), format="jd")
    airmass = calculate_airmass(obstime, observer, target)
    snr = estimate_snr(magnitude, exptime, airmass)

    # Save results for later
    result = {
        "name": [name] * n_eclipses,
        "snr": snr.to_value(1),
        "exptime": [eclipse_duration.to_value(u.second)] * n_eclipses,
        "time": obstime.mjd,
        "time_begin": observable_eclipses[:, 0].datetime,
        "time_end": observable_eclipses[:, 1].datetime,
        "stellar_effective_temperature": [data["st_teff"]] * n_eclipses,
    }
    return result


def transit_calculation(
    planets_or_criteria,
    date_start,
    date_end,
    constraints=None,
    observer="paranal",
    catalog="nexa",
    verbose=0,
    mode="planets",
    parallel=True,
):
    """
    Determine the SNR of all observable transits for one or more planets between two times

    Parameters
    ----------
    planets : list, str
        planet names in the format 'star letter'
    date_start : Time
        earliest time for a possible transit
    date_end : Time
        last time for a possible transit
    constraints : list, optional
        list of astroplan constraints. If not set uses a set of reasonable constraints.
    observer : str, ap.Observer, optional
        the location of the observatory/telescope, by default "paranal"
    catalog : str, optional
        the name of the data catalog, by default "nexa", which stands for the Nasa Exoplanet Archive
    verbose : int, optional
        amount of information displayed during runtime, 0 (default), 1 (some info), 2 (technical details), -1 (none)
    parallel : bool, optional
        if True will perform the observability and SNR calculation for all planets in parallel, by default True.

    Returns
    -------
    df : pd.DataFrame
        the results with one entry per transit
        entries are
          - name: name of the star and planet letter
          - snr: estimated snr for this transit
          - exptime: exposure time in seconds
          - time: mid transit time in MJD
          - time_begin: ingress time in MJD
          - time_end: egress time in MJD
    """

    assert mode in [
        "planets",
        "criteria",
    ], f"Expected one of ['planets', 'criteria'], but got {mode} for parameter 'mode'"
    assert catalog in [
        "nexa"
    ], f"Expected one of ['nexa',], but got {catalog} for parameter 'catalog'"

    # Fix the inputs to match our expectations
    if constraints is None:
        constraints = get_default_constraints()

    if not isinstance(observer, ap.Observer):
        observer = ap.Observer.at_site(observer)
    if isinstance(planets_or_criteria, str):
        planets_or_criteria = [planets_or_criteria]

    if planets_or_criteria is None:
        if mode == "criteria":
            if verbose >= 0:
                logger.warning(
                    "No criteria were set for the observation, using default values instead. Set 'planets_or_criteria' to an empty list if you want to use all planets instead []."
                )
            planets_or_criteria = [
                "pl_bmassj < 1",
                "dec < 10",
                "sy_jmag < 10",
                "sy_hmag < 10",
                "st_teff < 6000",
            ]
        else:
            raise TypeError(
                "Expected type str or list for parameter 'planets_or_criteria', but got None"
            )

    date_start = Time(date_start)
    date_end = Time(date_end)
    assert date_start < date_end, "Starting date must be earlier than the end date"

    if verbose >= 0:
        if mode == "planets":
            logger.info(
                f"Calculating observable transits for planets: {planets_or_criteria}"
            )
        elif mode == "criteria":
            logger.info(
                f"Calculating observable transits for planets that fullfill: {planets_or_criteria}"
            )
        else:
            raise ValueError
        logger.info(f"between {date_start} and {date_end}")
    if observer != "paranal" or verbose >= 1:
        if observer.name is not None:
            observer_name = observer.name
        else:
            observer_name = observer.location
        logger.info(f"at observing location {observer_name}")
    if catalog != "nexa" or verbose >= 1:
        logger.info(f"using data catalog {catalog}")

    if catalog == "nexa":
        # TODO: this is slow, add a cache?
        nexa = NasaExoplanetsArchive(verbose=verbose)
        if mode == "planets":
            catalog = nexa.query_objects(planets_or_criteria)
        elif mode == "criteria":
            catalog = nexa.query_criteria(planets_or_criteria)
        else:
            raise ValueError
    else:
        raise ValueError

    # Drop masked values
    mask = np.ma.getmaskarray(catalog["pl_tranmid"])
    mask |= np.ma.getmask(catalog["pl_orbper"])
    mask |= np.ma.getmask(catalog["pl_trandur"])
    mask |= np.ma.getmask(catalog["sy_kmag"])
    mask |= np.ma.getmask(catalog["st_teff"])
    catalog = catalog[~mask]

    if verbose >= 0:
        logger.info(
            f"Successfully loaded planet data from catalog and found {len(catalog)} potential planets"
        )

    results = {
        "name": [],
        "snr": [],
        "exptime": [],
        "time": [],
        "time_begin": [],
        "time_end": [],
        "stellar_effective_temperature": [],
    }
    if parallel:
        with ProcessPoolExecutor() as executor:
            futures = []
            for data in catalog:
                args = (data,
                        date_start,
                        date_end,
                        constraints,
                        observer,
                        verbose,)
                futures += [
                    executor.submit(
                        single_planet,
                        data,
                        date_start,
                        date_end,
                        constraints,
                        observer,
                        verbose,
                    )
                ]
            for future in tqdm(
                as_completed(futures),
                desc="Planets",
                total=len(futures),
                disable=verbose < 0,
            ):
                results_single = future.result()
                if results_single is not None:
                    for key in results.keys():
                        results[key] += [results_single[key]]
    else:
        for data in tqdm(catalog, desc="Planets"):
            results_single = single_planet(
                data, date_start, date_end, constraints, observer, verbose
            )
            if results_single is not None:
                for key in results.keys():
                    results[key] += [results_single[key]]

    # Store everything in a dataframe
    if len(results["name"]) == 0:
        raise ValueError("No observable transits found")
    for column in results.keys():
        results[column] = np.concatenate(results[column])
    df = pd.DataFrame(data=results)
    if verbose >= 1:
        logger.debug(df)
    return df


def main():
    """
    Command line interface for the crires planning tool
    all parameters are set via sys.argv

    Returns
    -------
    df : pd.DataFrame
        The dataframe containing all observable transits
    """

    parser = argparse.ArgumentParser(description="CRIRES+ planning tool")
    # Add all arguments
    parser.add_argument(
        "begin", type=dateparser.parse, help="first date to check for transits"
    )
    parser.add_argument(
        "end", type=dateparser.parse, help="last date to check for transits"
    )
    parser.add_argument(
        "planets",
        type=str,
        nargs="*",
        help="Names of the planets to calculate in the format 'star letter'",
    )
    parser.add_argument(
        "-O",
        "--observer",
        type=str,
        help="Location of the observer, by default 'paranal'",
        default="paranal",
    )
    parser.add_argument(
        "-C",
        "--catalog",
        type=str,
        help="Name of the data catalog to use, by default 'nexa' (Nasa Exoplanet Archive)",
        default="nexa",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        help="Which mode to use 'planets' or 'criteria'",
        default="planets",
    )
    parser.add_argument("-o", "--output", type=str, help="Save the data to this file")
    parser.add_argument(
        "-f",
        "--file",
        action="store_true",
        help="Load planet names from a file instead, one planet name per line. Lines starting with # are considered comments",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="If set will create an interactive plot",
    )
    parser.add_argument("--plot-file-1", help="filename for the first plot")
    parser.add_argument("--plot-file-2", help="filename for the first plot")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="verbosity level up to -vv")
    parser.add_argument("-s", "--silent", action="store_true", help="removes all console output, except for errors")
    args = parser.parse_args()

    # Process the arguments
    verbose = args.verbose if args.silent is False else -1
    planets = args.planets
    date_start = args.begin
    date_end = args.end
    observer = args.observer
    catalog = args.catalog
    mode = args.mode
    output = args.output
    plot = args.plot
    plot_file_1 = args.plot_file_1
    plot_file_2 = args.plot_file_2
    if plot_file_1 is not None or plot_file_2 is not None:
        plot = True


    if args.file:
        filename = planets[0]
        if verbose >= 1:
            print(f"Reading input planets/criteria from file {filename}")
        with open(filename, "r") as f:
            planets = f.readlines()
            planets = [p.strip() for p in planets]
            planets = [p for p in planets if not p.startswith("#")]

    if len(planets) == 0:
        planets = None

    if verbose < 0:
        warnings.filterwarnings("ignore", category=UserWarning, append=True)

    # Run the main method
    df = transit_calculation(
        planets,
        date_start,
        date_end,
        observer=observer,
        catalog=catalog,
        verbose=verbose,
        mode=mode,
    )

    # Output as desired
    if output is not None:
        df.to_csv(output)
        if verbose >= 1:
            print(f"Stored results in file {output}")
    if verbose == 0 and output is None:
        # for verbose >= 1, we already print it in the method
        print(df)

    if plot:
        create_interactive_graph(df, plot_file_1, plot_file_2)
    return df

if __name__ == "__main__":
    main()
