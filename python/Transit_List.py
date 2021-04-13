#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:55:10 2020

Original File Transit List. This file contains all the routines to use the functions, classes and class methods to compute observability of targets,
transits of exoplanets and to call the Exposure Time Calculator (ETC) from ESO to calculate S/N signal to noise ratio for observations with CRIRES+.

The different functionalities can be accessed via a menu popping up when running this file.

More functionalities can of course be added and the ETC part can be extended to other instruments used at the VLT
in Paranal, Chile.

Part of this tool could also be extracted and used for other observatories like the first part about observability
under certain constraints.

This tool was created as part of my master thesis: "Planning observations of terrestrial Exoplanets around M type Stars
with CRIRES+" during the peculiar times of Covid-19 in 2020.

Documentation can be found in the README file of this bundle or on GitHub or in my master thesis which can be found here: "COMING SOON".

Examples on how to use:
-----------------------

call Transit_List.py - The tool will welcome you and check if you have internet connection. If your internet connection is working,
you should see now a menu like:

    Choose one of the following options:
 1: Run full transit calculation
 2: run call ETC part for a list of transits
 3: run single transit planning
 4: run single target planning
 5: Plotting data of some result file
Enter number:


The different obtions will start the following procedures
1: Run full transit calculation -

2: run call ETC part for a list of transits -

3: run single transit planning -

4: run single target planning -

5: Plotting data of some result file -




@author: jonaszbinden
GitHub: jonaszubindu
"""


import copy
import datetime
import json
import logging
import os
import sys

import astroplan
import astropy
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astroplan import Observer, download_IERS_A, get_IERS_A_or_workaround
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon, get_sun
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.iers import IERS_Auto
from astropy.visualization import astropy_mpl_style, quantity_support
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

import classes_methods.Helper_fun as fun

from tqdm import tqdm
from concurrent.futures import (
    ThreadPoolExecutor,
    ProcessPoolExecutor,
    as_completed,
    wait,
)
from classes_methods.classes import (
    Eclipses,
    Exoplanets,
    Nights,
    load_eclipses_from_file,
)
from classes_methods.misc import misc

# """ Update most recent IERS data """
# iers.Conf.iers_auto_url.set(
#     "https://datacenter.iers.org/data/9/finals2000A.all"
# )  # temporary fix for iers data
# download_IERS_A(show_progress=True)
IERS_Auto()

plt.style.use(astropy_mpl_style)
quantity_support()

# TODO: make this copy prints (see maybe PySME for an implementation)
logging.basicConfig(
    filename="Transit_List.log",
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s-%(levelname)s-%(message)s",
)
logger = logging.getLogger(__name__)


def check_connection(
    host="http://exoplanetarchive.ipac.caltech.edu/",
):  # Nasa Exoplanet Archive
    """Check Internet Connection to Nasa Exoplanet Archive"""
    req = requests.get(host)  # Python 3.x
    if req.ok:
        print("Connected to {}".format(host))
    else:
        raise Warning(
            "Check connection to {}, response code:{}".format(host, req.status_code)
        )
    return req


def get_default_constraints():
    """ Altitude constraints definition """
    altcons = astroplan.AltitudeConstraint(min=+30 * u.deg, max=None)

    """ Airmass constraints definition """
    airmasscons = astroplan.AirmassConstraint(min=None, max=1.7)

    """ Astronomical Nighttime constraints definition: begin and end of each night at paranal as AtNightConstraint.twilight_astronomical """
    night_cons_per_night = astroplan.AtNightConstraint.twilight_astronomical()

    """ Moon Constraint """
    mooncons = astroplan.MoonSeparationConstraint(min=+45 * u.deg, max=None)

    constraints = [night_cons_per_night, altcons, airmasscons, mooncons]
    return constraints


def parse_date(date, max_delta_days):
    max_delta_days = int(max_delta_days)
    if date == "" or date is None:
        date = datetime.date.today()
    else:
        date = datetime.date.fromisoformat(date)

    dt = datetime.timedelta(days=1)
    date_end = max_delta_days * dt + date

    return date, max_delta_days, date_end


def parallel_func(planet, max_delta_days, obs_time, nights_paranal, constraints):
    """ This is used in  full_transit_calculation for parrallel processing """
    eclipse = Eclipses(max_delta_days, planet)
    eclipse.observability(
        obs_time, nights_paranal, constraints=constraints, check_eclipse=True,
    )
    return eclipse


def full_transit_calculation(date, max_delta_days, constraints, catalog="nexa_new"):

    date, max_delta_days, d_end = parse_date(date, max_delta_days)

    print(
        f"*** Running full transit analysis for transits between {date} and {d_end} ***"
    )

    # initialize list to sort exoplanet candidates and check if data are available to calculate observability.
    exoplanets = Exoplanets()
    exoplanets.load_catalog(catalog)
    """ Check all Planets if they have the necessary properties in the database to process for transit observations """
    exoplanets.filter_data()
    # work with Parse_list_transits in the same way as in the Transit_Example to retrieve all necessary transit information
    n_planets = len(exoplanets.data_complete)
    print(f"Found {n_planets} planets in Nasa Archive with Transit data")

    """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    nights_paranal = Nights(date, max_delta_days, load_from_pickle=0)

    """ Running methods to compute observability of single transits during the given timespan. """
    midnight = datetime.time(0, 0, 0)
    obs_time = Time(datetime.datetime.combine(nights_paranal.date[0], midnight))

    eclipses_list = []
    with ProcessPoolExecutor() as executor:
        futures = []
        for _, planet in tqdm(
            exoplanets.data_complete.iterrows(),
            total=n_planets,
            desc="Launched Processes",
        ):
            futures += [
                executor.submit(
                    parallel_func,
                    planet,
                    max_delta_days,
                    obs_time,
                    nights_paranal,
                    constraints,
                )
            ]

        for future in tqdm(
            as_completed(futures), total=n_planets, desc="Collected Processes"
        ):
            result = future.result()
            if len(result.eclipse_observable) != 0:
                eclipses_list += [result]

    return eclipses_list


def call_etc_for_list_of_transits(date, max_delta_days, filename):
    date, max_delta_days, d_end = parse_date(date, max_delta_days)

    date = datetime.date.fromisoformat(date)
    eclipses_list = load_eclipses_from_file(filename, max_delta_days)
    return eclipses_list


def etc_calculator(eclipses_list, minimum_snr=100):
    """
    For each observable planet:

    Do stuff with the input file 'etc-form.json' here:
    use: ETC.update_etc_form(**kwargs) from Etc_form_class

    Then write the whole file again as a json file: etc-form.json
    with ETC.write_etc_format_file()
    and run it with Etc_form_class.run_etc_calculator

    """

    """
    Calculates the median of the signal to noise ratio achievable in transits that allow around or more than 20 single exposures
    during the transit. The Number of exposures possible is added to the list eclipse_observable and eclipse_mid_observable
    for comparison. Each exposure is optimised to have NDIT between 16 and 32 with a minimum S/N = 100. The resulting S/N ratios
    of each exposure are used to compute the overall median. More values like DIT, NDIT, SN of each exposure for each transit
    could be stored as well, but this has not been implemented yet. If one gets stuck in this loop due to internet connection or
    unexpected errors, one may rerun the code from here, instead of rerunning everything again.
    """

    # DEBUG: Non Parrallel version
    for i, planet in tqdm(enumerate(eclipses_list), total=len(eclipses_list)):
        for j, eclipse in tqdm(enumerate(planet.eclipse_observable), total=len(planet.eclipse_observable)):
            eclipses_list[i].eclipse_observable[j] = fun.snr_estimate_nexposures(eclipse, planet, snr=minimum_snr)

    # with ProcessPoolExecutor(max_workers=None) as executor:
    #     futures = {}
    #     for i, planet in enumerate(eclipses_list):
    #         for j, eclipse in enumerate(planet.eclipse_observable):
    #             futures[executor.submit(fun.snr_estimate_nexposures, eclipse, planet, snr=minimum_snr)] = (i, j)

    #     for future in tqdm(as_completed(futures), total=len(futures), desc="Eclipses"):
    #         i, j = futures[future]
    #         eclipses_list[i].eclipse_observable[j] = future.result()

    return eclipses_list


def single_transit_calculation(
    date, max_delta_days, name, constraints, catalog="nexa_new", minimum_SN=100
):
    date, max_delta_days, d_end = parse_date(date, max_delta_days)

    print(
        f"*** Running single transit analysis for transits between {date} and {d_end} ***"
    )

    exoplanets = Exoplanets()
    exoplanets.load_catalog(catalog)
    exoplanets.filter_data()

    select = exoplanets.data_complete["pl_name"] == name
    planet = exoplanets.data_complete[select]
    if len(planet) == 0:
        raise Warning(f"Planet with name {name} not found in {catalog}")
    planet = planet.iloc[0]

    """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
    nights_paranal = Nights(date, max_delta_days, load_from_pickle=0)

    """ ETC part for single candidate exposure time calculation and number of possible exposures. """
    midnight = datetime.time(0, 0, 0)
    obs_time = Time(datetime.datetime.combine(nights_paranal.date[0], midnight))
    eclipse = Eclipses(max_delta_days, planet)
    eclipse.observability(
        obs_time, nights_paranal, constraints=constraints, check_eclipse=1
    )

    if len(eclipse.eclipse_observable) == 0:
        raise Warning(f"No observable eclipse found for {name}")

    with ThreadPoolExecutor(max_workers=None) as executor:
        futures = []
        for ecl in eclipse.eclipse_observable:
            futures += [
                executor.submit(
                    fun.snr_transit_observation_optimization,
                    ecl,
                    planet,
                    snr=minimum_SN,
                )
            ]
        wait(futures)

    # Store final Data
    name = eclipse.name.split(" ")
    name = name[0] + name[1]
    filename = f"{name}_events_processed_{date}_{max_delta_days}d.pkl"
    fun.save_pickled(filename, eclipse)
    return eclipse


if __name__ == "__main__":
    check_connection()

    """ Ask for menu input """
    k = misc.user_menu(
        menu=(
            "Run full transit calculation",
            "Run call ETC part for a list of transits",
            "Run single transit planning",
            "Run single target planning",
            "Plotting data of some result file",
        )
    )

    """ Location and UTC offset Paranal """
    paranal = Observer.at_site("paranal", timezone="Chile/Continental")

    constraints = get_default_constraints()

    """ Catalog to get planetary data from, nexa_old -> provided catalog by astroquery.NasaExoplanetArchive or nexa_new -> alpha version of new catalog: Planetary Systems Composite Data"""
    catalog = "nexa_new"

    minimum_snr = 100
    use_etc_calculator = "n"

    ##########################################################################################################

    if k == 1:
        """ Computes the observability of a list of Candidates """
        date = misc.ask_for_value(
            msg="Enter start date like 2020-06-01 or press enter to use todays date "
        )
        max_delta_days = misc.ask_for_value(
            msg="Enter number of days to compute transits for "
        )
        use_etc_calculator = misc.ask_for_value(
            msg="Do you want to call the ETC calculator to process the results S/N ratio? (WARNING : Only works with stable internet connection!) y/n "
        )

        full_transit_calculation(date, max_delta_days, constraints, catalog=catalog)

    ##########################################################################################################

    if k == 2:
        use_etc_calculator = "y"
        date = misc.ask_for_value(
            msg="Enter date like 2020-05-28 of the file of transit data you want to use "
        )
        max_delta_days = misc.ask_for_value(
            msg="Enter timespan in days to the appropriate file "
        )
        filename = f"Eclipse_events_processed_{date}_{max_delta_days}d.pkl"
        ans = misc.ask_for_value(
            msg=f"Do you want to load file : {filename} to feed to ETC? [y,n]"
        )
        if ans == "n":
            filename = misc.ask_for_value(msg="Enter filename with data to plot:  ")
        elif ans == "y":
            pass
        else:
            print("Gosh, what do you want then???")
            sys.exit()

        eclipses_list = call_etc_for_list_of_transits(date, max_delta_days, filename)

    ##########################################################################################################
    """ ETC part to process list of observable candidates """
    if use_etc_calculator == "y":
        eclipses_list = etc_calculator(eclipses_list, minimum_snr=minimum_snr)
        filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(d, max_delta_days)
        fun.save_pickled(filename, eclipses_list)

    ##########################################################################################################

    if k == 3:

        """ Only run single Planet candidate for observability and compute exact number of possible exposures thorugh the ETC """

        date = misc.ask_for_value(
            msg="Enter start date like 2020-06-01 or press enter to use todays date "
        )
        max_delta_days = misc.ask_for_value(
            msg="Enter number of days to compute transits for "
        )
        name = misc.ask_for_value(
            msg="Enter name like KELT-10 b of planet you want check observability "
        )

        Planet = single_transit_calculation(
            date, max_delta_days, name, constraints, eclipses_list
        )

    ##########################################################################################################

    if k == 4:

        """ Run analysis of single target for its observability, if this should also be available for target lists, I need some more information about where these lists would come from. """

        print("Coming soon")
        sys.exit()

        # print('Run through target properties manually or load from name...')
        # name = misc.ask_for_value(msg='name of target ')
        # target = NasaExoplanetArchive.query_star(name)

        # d = misc.ask_for_value(msg='name of target ')
        # Max_Delta_days = misc.ask_for_value(msg='name of target ')

        # st_Teff = misc.ask_for_value(msg='name of target ')
        # st_jmag = misc.ask_for_value(msg='name of target ')
        # number_of_time_steps_per_night = misc.ask_for_value(msg='how many time steps per night, press enter for default value=1000')
        # DIT = misc.ask_for_value(msg='name of target ')
        # NDIT = misc.ask_for_value(msg='name of target ')

        # """ Generates the class object Nights and calculates the nights for paranal between d and d_end """
        # Nights_paranal = Nights(d, Max_Delta_days, LoadFromPickle=0)

        # if delta_midnight == '':
        #     delta_midnight = np.linspace(-12, 12, 1000) * u.hour # defines number of timesteps per 24 hours

        # """
        # Some docs here
        # """

        # target = Target(name, star_Teff, star_jmag)
        # target.target_observability(Nights_paranal, constraints=constraints, delta_midnight=delta_midnight)

        # obs_time = misc.ask_for_value(msg='Choose a observation time for which you want the signal to noisu ratio S/N')
        # try :
        #     fun.SN_Ratio_Target(obs_time, target)
        # except Warning as w:
        #     print(w)
        #     print('Something went wrong in:{}:{}, taking next observation...'.format(Planet.name, eclipse['obs time']))
        #     logging.exception(w)
        #     logging.error('Something went wrong in:{}:{}, taking next observation...'.format(Planet.name, eclipse['obs time']))
        #     break
        # except Exception as e:
        #     #go back to menu
        #     print(e)
        #     logging.exception(e)
        #     print('\a')
        #     print('Catched random exception')

        #     print('Shall we save what has been computed so far to a picklefile? You may load that pickle file anytime later and continue from there. Just use the function pickled_items to load manually and investigate the output with next(output) or use load_planets_from_pickle to generate a list of Eclipses instances again like Eclipses_List')
        #     save = input('Do you want to save? y/n ')
        #     if save == 'y':
        #         d = d.isoformat() # start date from which the nights in paranal are calculated
        #         filename = 'Eclipse_events_processed_{}_{}d.pkl'.format(d, Max_Delta_days)
        #         fun.pickle_dumper_objects(filename, Eclipses_List)
        #         sys.exit()
        #     elif save == 'n':
        #         sys.exit()
        #     else:
        #         sys.exit()

    ##########################################################################################################

    if k == 1 and use_etc_calculator == "n":
        sys.exit()

    if k in [1, 2, 3, 4]:
        """ Storing data and plotting data """
        if k in [1, 2, 4]:
            data = eclipses_list
        elif k == 3:
            data = Planet
        else:
            raise Exception("What happened?")

        date, max_delta_days, d_end = parse_date(date, max_delta_days)
        filename = f"Eclipse_events_processed_{date}_{max_delta_days}d.pkl"

        ranking, df_gen, df_frame, num_trans = fun.data_sorting_and_storing(
            data, filename, write_to_csv=1
        )
        ranked_events, Obs_events = fun.postprocessing_events(
            date, max_delta_days, Nights, data
        )
        fun.xlsx_writer(filename, df_gen, df_frame, Obs_events)

    k2 = misc.user_menu(
        menu=(
            "Plot candidates over full period",
            "Plot single night of (mutual) target(s)",
            "Get target finder image ",
        )
    )

    if k == 5:
        """ Plotting data of some result file """

        if k2 in [1, 2]:
            filename = misc.ask_for_value(
                msg="Enter filename with data to plot:  [only works with non customized file names]"
            )
            if filename.split(".")[-1] == "pkl":
                pass
            else:
                filename = filename.split(".")[0] + ".pkl"

            d = datetime.date.fromisoformat(filename.split("_")[-2])
            Max_Delta_days = int((filename.split("_")[-1].split(".")[0]).split("d")[0])
            eclipses_list = load_eclipses_from_file(filename, Max_Delta_days)

    if k2 == 1 and k == 5:
        """ Plotting candidates over full period """
        ranking, df_gen, df_frame, _ = fun.data_sorting_and_storing(
            eclipses_list, write_to_csv=0
        )
        ranked_events, Obs_events = fun.postprocessing_events(
            d, Max_Delta_days, Nights, eclipses_list
        )
        fun.xlsx_writer(filename, df_gen, df_frame, Obs_events)
        ranking = fun.plotting_transit_data(
            d, Max_Delta_days, ranking, eclipses_list, Nights, ranked_events
        )
    elif k2 == 1:
        ranking = fun.plotting_transit_data(
            d, Max_Delta_days, ranking, eclipses_list, Nights, ranked_events
        )

    if k2 == 2:
        """ Plot single night of (mutual) target(s) """
        d = misc.ask_for_value(
            msg="Enter date like 2020-05-28 of the night you want to investigate, CAUTION: dates are regarded with respect to UTC "
        )
        d = datetime.date.fromisoformat(d)
        midn = datetime.time(0, 0, 0)
        d = datetime.datetime.combine(d, midn)
        name = misc.ask_for_value(
            msg=f"Do you want to plot all mutually observable targets for the night of the {d}? Press enter, otherwise write the name of the target you want to plot "
        )
        found = 0
        if name != "":
            # if name == target.name:
            #     fun.plot_night(d, location = paranal.location, obs_obj = target)
            # else:
            for planet in eclipses_list:
                if planet.name == name:
                    fun.plot_night(d, location=paranal.location, obs_obj=planet)
                    found = 1
        if found == 0:
            print("Did not get valid name, plotting all candidates...")
            fun.plot_night(d, location=paranal.location, obs_obj=eclipses_list)

    if k2 == 3:
        """ Get target finder image """
        name = misc.ask_for_value(msg="Write the name of the target you want to plot ")
        try:
            fun.find_target_image(name)
        except NameError:
            name = misc.ask_for_value(msg=f"{name} does not exist, try again ")
            fun.find_target_image(name)
