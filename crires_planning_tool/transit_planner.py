from os.path import join, dirname

import astroplan as ap
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm


try:
    from .nasa_exoplanets_archive import NasaExoplanetsArchive
except ImportError:
    from nasa_exoplanets_archive import NasaExoplanetsArchive


def get_default_constraints():
    """ Altitude constraints definition """
    altcons = ap.AltitudeConstraint(min=+30 * u.deg, max=None)

    """ Airmass constraints definition """
    airmasscons = ap.AirmassConstraint(min=None, max=1.7)

    """ Astronomical Nighttime constraints definition: begin and end of each night at paranal as AtNightConstraint.twilight_astronomical """
    night_cons_per_night = ap.AtNightConstraint.twilight_astronomical()

    """ Moon Constraint """
    mooncons = ap.MoonSeparationConstraint(min=+45 * u.deg, max=None)

    constraints = [night_cons_per_night, altcons, airmasscons, mooncons]
    return constraints


def load_catalog(catalog="nexa"):
    filenames = {
        "nexa": "PlanetList.csv",
        "custom": "PlanetList_edit.csv",
    }
    # Use a standard catalog if available
    # Otherwise interpret it as a filename
    try:
        filename = join(dirname(__file__), "csv_files", filenames[catalog])
    except KeyError:
        filename = catalog

    candidate_list = pd.read_csv(filename)
    return candidate_list


def calculate_airmass(obstime, observer, target):
    altaz = observer.altaz(obstime, target)
    airmass = altaz.secz
    return airmass


@u.quantity_input
def estimate_snr(
    magnitude: u.mag, exptime: u.second, airmass: u.one, method="crires", mband="K"
):
    """
    Computes SNR of stellar spectrum from magnitude, exposure time (in seconds) and airmass
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
            print("Use Jmag for calculation of Carmenes SNR.")
            snr = np.nan
    elif method == "crires":
        if mband == "K":
            snr = 449.4241 * np.sqrt(10 ** (-mag / 2.5)) * np.sqrt(et) - 6.3144
        else:
            print("Use Kmag for calculation of Crires SNR.")
            snr = np.nan

    else:
        print("Method not recognized. Use Crires or Carmenes.")
        snr = np.nan

    # Scale to airmass = airm
    extcof = 0.05  # extinction coefficient, see Lombardi et al., 2018
    snr *= 10 ** (extcof / 5 * (1 - airmass))

    return snr


def transit_calculation(
    planets,
    date_start,
    date_end,
    constraints=None,
    observer="paranal",
    catalog="nexa",
    verbose=0,
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

    assert date_start < date_end, "Starting date must be earlier than the end date"

    # Fix the inputs to match our expectations
    if constraints is None:
        constraints = get_default_constraints()

    if isinstance(observer, str):
        observer = ap.Observer.at_site(observer)
    if isinstance(planets, str):
        planets = [planets]

    date_start = Time(date_start)
    date_end = Time(date_end)

    if verbose >= 0:
        print(f"Calculating observable transits for planets: {planets}")
        print(f"between {date_start} and {date_end}")
    if observer != "paranal" or verbose >= 1:
        print(f"at observing location {observer}")
    if catalog != "nexa" or verbose >= 1:
        print(f"using data catalog {catalog}")

    if catalog == "nexa":
        # TODO: this is slow, add a cache?
        nexa = NasaExoplanetsArchive(verbose=verbose)
        catalog = nexa.query_objects(planets)
    else:
        raise ValueError(
            f"Could not parse selected catalog, expected 'nexa', but got {catalog} instead"
        )

    if verbose >= 1:
        print("Successfully loaded planet data from catalog")

    results = {
        "name": [],
        "snr": [],
        "exptime": [],
        "time": [],
        "time_begin": [],
        "time_end": [],
    }
    with tqdm(catalog, desc="Planets", disable=verbose < 0) as progress:
        for data in progress:
            primary_eclipse_time = Time(data["pl_tranmid"], format="jd")
            orbital_period = data["pl_orbper"] * u.day
            eclipse_duration = data["pl_trandur"] * u.hour
            name = data["pl_name"]
            # Update the progress bar
            progress.desc = f"Planet {name}"
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
            eclipses = system.next_primary_ingress_egress_time(
                date_start, n_eclipses=n_max_eclipses
            )
            if eclipses[-1][0] > date_end:
                eclipses = eclipses[:-1]

            # Check if they are observable
            is_observable = ap.is_event_observable(
                constraints, observer, target, times_ingress_egress=eclipses
            )[0]
            observable_eclipses = eclipses[is_observable]

            # Calculate the SNR for those that are left
            magnitude = data["sy_kmag"] * u.mag
            exptime = eclipse_duration
            # exptime = np.diff(observable_eclipses.jd, axis=1).ravel() * u.day
            obstime = Time(np.mean(observable_eclipses.jd, axis=1), format="jd")
            airmass = calculate_airmass(obstime, observer, target)
            snr = estimate_snr(magnitude, exptime, airmass)

            # Save results for later
            n_eclipses = len(snr)
            if n_eclipses == 0:
                if verbose >= 0:
                    print(f"WARNING: No observable transits found for planet {name}")
                continue
            results["name"] += [[name] * n_eclipses]
            results["snr"] += [snr.to_value(1)]
            results["exptime"] += [[eclipse_duration.to_value(u.second)] * n_eclipses]
            results["time"] += [obstime.mjd]
            results["time_begin"] += [observable_eclipses[:, 0].mjd]
            results["time_end"] += [observable_eclipses[:, 1].mjd]
        progress.desc = "Planets"
    # Store everything in a dataframe
    if len(results["name"]) == 0:
        raise ValueError("No observable transits found")
    for column in results.keys():
        results[column] = np.concatenate(results[column])
    df = pd.DataFrame(data=results)
    if verbose >= 1:
        print(df)
    return df


if __name__ == "__main__":
    # TODO: add a menu
    # one planet all transits in a time frame
    # multiple planets in a time frame
    import argparse
    import dateparser

    parser = argparse.ArgumentParser(description="CRIRES+ planning tool")
    parser.add_argument(
        "start", type=dateparser.parse, help="first date to check for transits"
    )
    parser.add_argument(
        "end", type=dateparser.parse, help="last date to check for transits"
    )
    parser.add_argument(
        "planets",
        type=str,
        nargs="+",
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
    parser.add_argument("-o", "--output", type=str, help="Save the data to this file")
    parser.add_argument(
        "-f",
        "--file",
        action="store_true",
        help="Load planet names from a file instead, one planet name per line. Lines starting with # are considered comments",
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("-s", "--silent", action="store_true")
    args = parser.parse_args()

    verbose = args.verbose if args.silent is False else -1
    planets = args.planets
    date_start = args.start
    date_end = args.end
    observer = args.observer
    catalog = args.catalog

    if args.file:
        filename = planets[0]
        if verbose >= 1:
            print(f"Reading input planets from file {filename}")
        with open(filename, "r") as f:
            planets = f.readlines()
            planets = [p.strip() for p in planets]
            planets = [p for p in planets if not p.startswith("#")]

    # Run the main method
    df = transit_calculation(
        planets,
        date_start,
        date_end,
        observer=observer,
        catalog=catalog,
        verbose=verbose,
    )

    if args.output is not None:
        df.to_csv(args.output)
        if verbose >= 1:
            print(f"Stored results in file {args.output}")
    if verbose == 0 and not args.output:
        # for verbose >= 1, we already print it in the method
        print(df)

    pass
