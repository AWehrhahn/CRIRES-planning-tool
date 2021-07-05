#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:53:09 2020

This file contains helper Functions used in Transit_List.py

@author: jonaszbinden
GitHub: jonaszubindu
"""
import copy
import datetime
import logging
import pickle
from os.path import dirname, join, realpath

import numpy as np
import pandas as pd
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_finder_image
from astropy import units as u
from astropy.coordinates import AltAz, get_moon, get_sun
from astropy.time import Time
from astropy.visualization import astropy_mpl_style, quantity_support
from matplotlib import pyplot as plt

try:
    from .classes import Nights
    from .etc_form import EtcForm
except ImportError:
    from classes_methods.classes import Nights
    from classes_methods.etc_form import EtcForm

plt.style.use(astropy_mpl_style)
quantity_support()
""" Location and UTC offset Paranal """
paranal = Observer.at_site("paranal", timezone="Chile/Continental")
logger = logging.getLogger(__name__)


def parse_date(date, max_delta_days):
    max_delta_days = int(max_delta_days)
    if date == "" or date is None:
        date = datetime.date.today()
    elif isinstance(date, datetime.date):
        pass
    elif isinstance(date, datetime.datetime):
        date = date.date()
    else:
        date = datetime.date.fromisoformat(str(date))

    dt = datetime.timedelta(days=1)
    date_end = max_delta_days * dt + date

    return date, max_delta_days, date_end

def etc_calculator_texp(obs_obj, obs_time, snr=100):
    """
        Optimizes NDIT for the S/N minimum defined by ''snr'' for a given DIT for a certain
        observation target ''obs_obj'' at a certain observation time ''obs_time''.

        Parameters
        ----------
        obs_obj : class object
            class object of a observation target.

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        snr : float
            Minimum S/N ratio that should be reached in a complete exposure

        Returns
        -------
        Exposure_time : float
            Time to compute a full set of exposures with NDIT single exposures of duration DIT.

        DIT : float
            Exposure time for single exposure.

        NDIT : float
            number of single exposures to create one complete set of exposure.

        output : namespace object
            Object containing the full output data from the ETC.

        ETC : etc_form object
            etc_form class instance with input data for the ETC.

    """
    ndit_opt = 24  # NDIT should optimally be between 16-32
    etc = EtcForm(inputtype="snr-Templ", instrument="crires")
    gsmag = obs_obj.star_jmag
    if gsmag < 9.3 * u.mag:
        gsmag = (
            9.3 * u.mag
        )  # old crires, brighter targets were dimmed by a filter down to 9.3
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(
        obs_obj, obs_time
    )  # add moon_target_sep
    input_data = etc.update_form(
        snr=snr,
        temperature=obs_obj.star_teff,
        brightness=obs_obj.star_jmag,
        airmass=airmass,
        moon_target_sep=moon_target_sep,
        moon_phase=moon_phase,
        gsmag=gsmag,
    )

    output = etc.call(input_data)
    ndit = etc.get_ndit(output)

    # TODO: can we predict this better, so we don't have to send several requests?
    # Routine to change ndit to 16-32 and change dit accordingly:
    cycles = 0
    while ndit < 16 or ndit > 32:
        exposure_time = ndit * input_data["timesnr"]["dit"]
        dit_new = exposure_time / ndit_opt  # determine DIT for NDIT=24
        input_data["timesnr"]["dit"] = dit_new
        logging.info("executed cycle:{}, DIT:{}, NDIT{}".format(cycles, dit_new, ndit))
        try:
            # recalculate the new NDIT
            output = etc.call(input_data)
            ndit = etc.get_ndit(output)
        except Warning:
            raise Warning("DIT seams not feasable input value: {}".format(dit_new))
        if ndit == 0:
            raise Warning("NDIT not available from etc-calculator")
        if cycles > 5:
            raise Warning("too many tries to bring NDIT between 16-32")
        cycles += 1

    dit = input_data["timesnr"]["dit"]
    exposure_time = ndit * dit  # seconds
    logging.info(
        f"Final values: Exposure time:{exposure_time}, DIT: {dit}, NDIT:{ndit}"
    )

    return exposure_time, dit, ndit, output, etc


def etc_calculator_snr(obs_obj, obs_time, ndit, dit):
    """
        Calculates solely the S/N ratio for a given ''dit'' and ''ndit'' for a certain observation
        target ''obs_obj'' at a certain observation time ''obs_time''. CAUTION: This function has not been used
        or tested much yet.

        Parameters
        ----------
        obs_obj : class object
            class object of an observation target.

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        ndit : int
            Number of frames taken during a full single exposure.
        dit : float
            Exposure time for each frame.

        Returns
        -------
        output : namespace object
            Object containing the full output data from the ETC.

        ETC : etc_form object
            etc_form class instance with input data for the ETC.
    """
    etc = EtcForm(inputtype="ndit-Templ", instrument="crires")
    gsmag = obs_obj.star_jmag
    if gsmag < 9.3 * u.mag:
        gsmag = 9.3 * u.mag  # Check again why that one is
    moon_target_sep, moon_phase, airmass, _ = airmass_moon_sep_obj_altaz(
        obs_obj, obs_time
    )  # add moon_target_sep
    input_data = etc.update_form(
        temperature=obs_obj.star_Teff,
        brightness=obs_obj.star_jmag,
        airmass=airmass,
        dit=dit,
        ndit=ndit,
        moon_target_sep=moon_target_sep,
        moon_phase=moon_phase,
        gsmag=gsmag,
    )

    output = etc.call(input_data)

    return output, etc


def calculate_snr_ratio(snr_data):
    """
        Calculates the median of the signal to noise S/N ratio data ''sn_data''.

        Parameters
        ----------
        sn_data : list
            Containing the S/N ratio data of which the median should be calculated.

        Returns
        -------
        median_SN : float
            Median of the S/N ratio data.
        min_SN : float
            minimum S/N.
        max_SN : float
            maximum S/N.
    """
    snr_data = np.asarray(snr_data)
    median_snr = np.median(snr_data)
    min_snr = np.min(snr_data)
    max_snr = np.max(snr_data)
    return median_snr, min_snr, max_snr


def extract_out_data(outputs):
    """
        Function to extract the S/N ratio data from the ''output'' file generated by the ETC.

        Parameters
        ----------
        outputs : namespace object or list
            Object or list of objects containing the full output data from the ETC.

        Returns
        -------
        SN_data : list
            Contains a list of all data from the output(s) of the ETC.
    """
    SN_data = []
    if isinstance(outputs, list):
        for output in outputs:
            for data in output["data"]["orders"]:
                for det in data["detectors"]:
                    SN_data.extend(det["data"]["snr"]["snr"]["data"])
    else:
        output = outputs
        for data in output["data"]["orders"]:
            for det in data["detectors"]:
                SN_data.extend(det["data"]["snr"]["snr"]["data"])
    return SN_data


def airmass_moon_sep_obj_altaz(obs_obj, obs_time, location=paranal.location):
    """
        This function calculates the moon target separation, moon phase (moon sun separation), airmass factor and local coordinates to observe
        the object ''obs_obj'' at ''obs_time'' at the location given in ''location'', which is normally paranal.

        Parameters
        ----------
        obs_obj : class object
            instance of class Eclipses, Can also be any other .

        obs_time : astropy.time.core.Time
            Time in UTC of observation.

        location : astropy.coordinates.EarthLocation, optional
            location of the observatory. The default is paranal.location.

        Returns
        -------
        moon_target_sep : float
            angular seperation between moon and target in degrees.

        moon_phase : float
            angular seperation between moon and sun in degrees.

        airmass : float
            Airmass factor at observation time and location.

        obs_altazs : astropy.coordinates.AltAz object
            Azimuth and Altitude in deg at which the object can be observed at the chosen time and location.

    """
    if isinstance(obs_obj, FixedTarget):
        obs_coor = obs_obj.coord
    else:
        obs_coor = obs_obj.coordinates.coord
    frame_obs = AltAz(obstime=obs_time, location=location)
    obs_altazs = obs_coor.transform_to(frame_obs)
    airmass = obs_altazs.secz
    moon = get_moon(obs_time).transform_to(frame_obs)
    sun = get_sun(obs_time).transform_to(frame_obs)
    # calculates the moon target separation
    moon_target_sep = moon.separation(obs_altazs)
    # moon_target_sep = moon_target_sep.deg * u.deg
    moon_phase = sun.separation(moon)
    # moon_phase = moon_phase_angle(time=obs_time)
    # moon_phase = moon_phase.to(u.deg)
    z = 90 * u.deg - obs_altazs.alt
    zmoon = 90 * u.deg - moon.alt
    sep_min = np.abs(z - zmoon)
    sep_max = np.abs(z + zmoon)
    moon_target_sep = np.clip(moon_target_sep, sep_min, sep_max)
    moon_target_sep = (moon_target_sep, moon.alt)

    return moon_target_sep, moon_phase, airmass, obs_altazs


def load_pickled(filename):
    """ Unpickle a file of pickled data. """
    path = join(dirname(__file__), "../picklefiles")

    with open(join(path, filename), "rb") as f:
        return pickle.load(f)[0]


def save_pickled(filename, *args):
    """
        Simple function to store class objects or list of class objects ''Objects'' as .pkl file under ''filename''.

        Parameters
        ----------
        filename : str
            filename under which the data should be stored.

        Objects : class object or list of class objects
            class object to store.

    """
    path = join(dirname(__file__), "../picklefiles")

    with open(join(path, filename), "wb") as out:
        pickle.dump(args, out)

    logger.debug(f"Successfully pickled file {filename}")


##########################################################################################################


def snr_transit_observation_optimization(eclipse, planet, snr=100):
    """
        Calculates exactly how many exposures are possible to take during a single transit and adds the data
        to the object Eclipses.eclipse_observable. This function gets only called for single targets cause
        it calls the ETC calculator several times. This is to make sure that changes in the exposure time
        during the transit can be accounted for.

        Parameters
        ----------
        eclipse : object from list Eclipses.eclipse_observable.
            contains all the data about a single transit observation.

        planet : instance of class Eclipses.
            Eclipses instance with all the data about the planet in question.
            
        snr : float
            Minimum S/N ratio that should be reached in a complete exposure    

    """

    """ Checking if eclipse has already been processed """
    try:
        test_empty = eclipse["n_exposures_possible"]
        if test_empty != None:
            print("{} has already been processed, skipping...".format(planet.name))
    except KeyError:
        print("{} gets fed to ETC calculator for best observations".format(planet.name))
        logging.info(
            "{} gets fed to ETC calculator for best observations".format(planet.name)
        )
        obs_time = eclipse["eclipse_mid"]["time"]

        Transit_dur = planet.transit_duration.to(u.second).value  # in seconds

        """ Initial Calculation of Exposure time for Transit: eclipse """
        # for different S/N ratio add argument 'snr'= , to every Etc_calculator_Texp function
        Exposure_time, _, _, output, _ = etc_calculator_texp(
            planet, obs_time
        )  # obtimising NDIT for each single exposure with S/N min = 100 in seconds

        """ Find how many single exposures are possible to take """

        Exposure_times = [Exposure_time]
        SN_data_overall = []
        num_exp = 1
        range_obs_times = Exposure_time
        time_between_exposures = 10  # buffer between two exposures in seconds
        n = 1
        median_SN_single_exp = []

        delta = datetime.timedelta(
            seconds=int(np.ceil(Exposure_time / 2 + time_between_exposures))
        )
        obs_time_up = obs_time + delta
        obs_time_down = obs_time - delta

        while Transit_dur > range_obs_times:
            print("number of exposures: {}".format(num_exp))
            Exposure_time_up, _, _, output, _ = etc_calculator_texp(
                planet, obs_time_up
            )  # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            Exposure_time_down, _, _, output, _ = etc_calculator_texp(
                planet, obs_time_down
            )  # obtimising NDIT for each single exposure with S/N min = 100 in seconds
            Exposure_times.append(Exposure_time_up)
            Exposure_times.append(Exposure_time_down)
            delta_up = datetime.timedelta(
                seconds=int(np.ceil(Exposure_time_up + time_between_exposures))
            )
            delta_down = datetime.timedelta(
                seconds=int(np.ceil(Exposure_time_down + time_between_exposures))
            )
            obs_time_up = obs_time_up + delta_up
            obs_time_down = obs_time_down - delta_down
            num_exp = 2 * n
            range_obs_times = (obs_time_up - obs_time_down).sec
            n += 1

            # Observations.extend(obs_times)
            SN_data = extract_out_data(output)
            SN_data_overall.extend(SN_data)
            median_SN, _, _ = calculate_snr_ratio(SN_data)
            median_SN_single_exp.append(median_SN)
        Exposure_times.sort()

        if num_exp < 20 and range_obs_times > Transit_dur:
            print(
                "Time to reach 20 exposure exceeded, number of possible exposures: {}".format(
                    num_exp
                )
            )
            eclipse["n_exposures_possible"] = num_exp
            eclipse[
                "Comment"
            ] = "Reaching 20 exposures with S/N = 100 exceeds Transit Duration"
        elif num_exp >= 20:  # and range_obs_times <= Transit_dur:
            print(
                "Reached 20 Exposures in: {} seconds during a {} seconds long Transit".format(
                    np.ceil(range_obs_times), Transit_dur
                )
            )
            Median_SN, _, _ = calculate_snr_ratio(SN_data_overall)
            median_SN_single_exp.sort()
            eclipse[
                "n_exposures_possible"
            ] = num_exp  # estimates the number of exposures possible according to the transit duration and the maximum exposure time calculated reaching 20 exposures
            eclipse["Time necessary to reach 20 exposures [s]"] = np.ceil(
                range_obs_times
            )
            eclipse["S/N overall median"] = Median_SN
            eclipse["Minimum S/N"] = median_SN_single_exp[0]
            eclipse["Maximum S/N"] = median_SN_single_exp[-1]
            eclipse["Minimum Exposure Time"] = Exposure_times[0]
            eclipse["Maximum Exposure Time"] = Exposure_times[-1]
            eclipse["List of Exposure Times"] = Exposure_times


##########################################################################################################


def snr_estimate_nexposures(eclipse, planet, snr=100):
    """
        Calculates the exposure time to reach 16 < NDIT < 32 for Transit mid, begin and end
        and from the maximum exposure time estimates the number of possible exposure during the whole transit.
        Also stores data about minimum, maximum and medium S/N ratio and adds the data
        to the object Eclipses.eclipse_observable.

        Parameters
        ----------
        eclipse : object from list Eclipses.eclipse_observable.
            contains all the data about a single transit observation.

        planet : instance of class Eclipses.
            Eclipses instance with all the data about the planet in question.
        
        snr : float
            Minimum S/N ratio that should be reached in a complete exposure

    """

    """ Checking if eclipse has already been processed """

    logging.info(f"{planet.name} gets fed to ETC calculator for best observations")

    obs_time = eclipse["eclipse_mid"]["time"]
    transit_dur = planet.transit_duration.to_value(u.second)
    # for different S/N ratio add argument 'snr'= , to every Etc_calculator_Texp function
    # obtimising NDIT for each single exposure with S/N min = 100 in seconds
    exposure_time, _, _, output, _ = etc_calculator_texp(planet, obs_time, snr=snr)

    SN_data = extract_out_data(output)
    median_SN, min_SN, max_SN = calculate_snr_ratio(SN_data)
    num_exp_possible = int(np.floor(transit_dur / exposure_time))
    # estimates the number of exposures possible according to the transit duration and the maximum exposure time
    eclipse["n_exposures_possible"] = num_exp_possible
    eclipse["snr_median"] = median_SN
    eclipse["snr_minimum"] = min_SN
    eclipse["snr_maximum"] = max_SN
    eclipse["average_exposure_time"] = exposure_time
    return eclipse


##########################################################################################################


def to_pandas(eclipses_List):
    """
        Store the data from the eclipses_list in a dataframe and/or a csv file

        Parameters
        ----------
        eclipses_list : list
            Contains all the objects from class Eclipses that have been loaded.

        filename : str, optional
            If this function is called independent of the observability runs, include filename
            from which Eclipse data should get loaded for sorting and processing.

        Returns
        -------
        frame : pd.Dataframe
            Data organized into a dataframe

    """
    frame = []
    for planet in eclipses_List:
        for eclipse in planet.eclipse_observable:
            data = {
                "name": eclipse["name"],
                "time": eclipse["obs_time"],
                "transit_duration": eclipse["transit_duration"],
                "stellar_effective_temperature": eclipse[
                    "stellar_effective_temperature"
                ].to_value("K"),
                "stellar_magnitude_j": eclipse["magnitude_j"].to_value("mag"),
                "n_exposures_possible": eclipse["n_exposures_possible"],
                "snr_median": eclipse["snr_median"],
                "snr_maximum": eclipse["snr_maximum"],
                "average_exposure_time": eclipse["average_exposure_time"],
                "is_observable": eclipse["is_primary_eclipse_observable"],
            }
            for key in ["begin", "mid", "end"]:
                ekey = f"eclipse_{key}"
                data[f"time_{key}"] = eclipse[ekey]["time"]
                data[f"airmass_{key}"] = eclipse[ekey]["airmass"].to_value(1)
                data[f"az_{key}"] = eclipse[ekey]["az"].to_value("deg")
                data[f"alt_{key}"] = eclipse[ekey]["alt"].to_value("deg")
                data[f"moon_sep_{key}"] = eclipse[ekey]["moon sep"].to_value("deg")
                data[f"moon_phase_{key}"] = eclipse[ekey]["moon phase"].to_value("deg")
            frame += [data]

    frame = pd.DataFrame(frame)
    return frame


##########################################################################################################


def plotting_transit_data(eclipses_list):
    """
        Plotting final data in ''Eclipses_List'' for the time span given in ''Nights'' or from date ''d'' for ''Max_Delta_days'' days.
        For now this only works for Eclipses, will include later functionality to plot general
        targets and maybe different functions to plot different kinds of data. Might contain more types of output
        than it has now.

        Parameters
        ----------
        Eclipses_List : list
            contains Eclipses class objects, which should be plotted.
    """

    if not isinstance(eclipses_list, list):
        eclipses_list = [eclipses_list]

    if not isinstance(eclipses_list, pd.DataFrame):
        eclipses_list = to_pandas(eclipses_list)

    x = eclipses_list["time"].mjd.values
    y = eclipses_list["n_exposures_possible"]
    plt.plot(x, y)
    plt.xlabel("Time")
    plt.ylabel("#Exposures")
    plt.show()


##########################################################################################################


def find_target_image(name):
    """
        Plotting function to get a quick find image for the target.

        Parameters
        ----------
        name : str
            Name of target.

        Returns
        -------
        Figure showing the target relative on the night sky.

    """
    messier1 = FixedTarget.from_name(name)
    ax, hdu = plot_finder_image(messier1)
    plt.show()


def plot_night(date, location, obs_obj, mix_types=1):
    """
        Plots the targets of a single night, depending where on the night sky they appear at which time.

        Parameters
        ----------
        date : datetime or astropy.time.Time
            date for which the night shall be plottet.

        location : astroplan.Observer.location or astropy.coordinates.EarthLocation
            Location of observatory.

        obs_obj : class object or list of class objects
            Class object containing information about coordinates, observation times.

        mix_types : int (obtional)
            set to zero if you want to only compare mutual transits in the same night.

    """
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = ax1.twiny()

    delta_midnight = np.linspace(-12, 12, 1000) * u.hour
    d = Time(date)
    night = d + delta_midnight
    frame_obs = AltAz(obstime=night, location=location)

    moon = get_moon(night)
    sun = get_sun(night)
    moon_altazs = moon.transform_to(frame_obs)
    sun_altazs = sun.transform_to(frame_obs)

    ax1.fill_between(
        delta_midnight,
        0 * u.deg,
        90 * u.deg,
        sun_altazs.alt < -0 * u.deg,
        color="0.5",
        zorder=0,
    )
    ax1.fill_between(
        delta_midnight,
        0 * u.deg,
        90 * u.deg,
        sun_altazs.alt < -18 * u.deg,
        color="k",
        zorder=0,
    )

    ax1.plot(delta_midnight, sun_altazs.alt, color=[0.75] * 3, label="Sun")
    ax1.plot(delta_midnight, moon_altazs.alt, color=[0.75] * 3, ls="--", label="Moon")

    no_list = 0
    warning = 0

    if type(obs_obj) == list and len(obs_obj) > 1:
        """ plotting for list of objects """
        no_ecl_observable = 0
        for obs_obj in obs_obj:
            if hasattr(
                obs_obj, "eclipse_observable"
            ):  # If object is related to eclipses
                for eclipse in obs_obj.eclipse_observable:
                    eclipse1 = copy.deepcopy(eclipse)
                    if eclipse1["eclipse_mid"]["time"].datetime.date() == date.date():
                        obs_time = eclipse1["eclipse_mid"]["time"]
                        t = obs_time.datetime.time()
                        h = (t.hour + t.minute / 60 + t.second / 3600) * u.hour
                        delta_eclipse = np.linspace(
                            h - obs_obj.transit_duration.to(u.hour) / 2,
                            h + obs_obj.transit_duration.to(u.hour) / 2,
                            100,
                        )
                        delta_eclipse_frame = np.linspace(
                            -obs_obj.transit_duration / 2,
                            +obs_obj.transit_duration / 2,
                            100,
                        )
                        transit = Time(obs_time) + delta_eclipse_frame
                        frame_ecl = AltAz(obstime=transit, location=location)
                        obs_ecl = obs_obj.Coordinates.coord.transform_to(frame_ecl)

                        obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                        if any(
                            obs_altazs[obs_altazs.alt > 85 * u.deg]
                        ):  # checking if target reaches zenith angle
                            warning = 1
                        im = ax1.scatter(
                            delta_midnight,
                            obs_altazs.alt,
                            label=obs_obj.name,
                            lw=0,
                            s=8,
                            cmap="viridis",
                            vmin=-10,
                            vmax=10,
                        )  # plot candidate
                        if eclipse["n_exposures_possible"] >= 20:
                            ax1.scatter(
                                delta_eclipse, obs_ecl.alt, color="red", lw=3, s=8
                            )
                        no_ecl_observable = 0
                    else:
                        no_ecl_observable = 1

            elif mix_types == 1:
                obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                im = ax1.scatter(
                    delta_midnight,
                    obs_altazs.alt,
                    label=obs_obj.name,
                    lw=0,
                    s=8,
                    cmap="viridis",
                    vmin=-10,
                    vmax=10,
                )  # Plot candidate

    elif type(obs_obj) == list and len(obs_obj) == 1:
        obs_obj = obs_obj[0]
        no_list = 1
    else:
        pass
        no_list = 1

    if no_list == 1:
        no_ecl_observable = 0
        """ Plotting for single objects """
        if hasattr(obs_obj, "eclipse_observable"):  # If object is related to eclipses
            for eclipse in obs_obj.eclipse_observable:
                eclipse1 = copy.deepcopy(eclipse)
                if eclipse1["eclipse_mid"]["time"].datetime.date() == date.date():
                    obs_time = eclipse1["eclipse_mid"]["time"]
                    t = obs_time.datetime.time()
                    h = (t.hour + t.minute / 60 + t.second / 3600) * u.hour
                    delta_eclipse = np.linspace(
                        h - obs_obj.transit_duration / 2,
                        h + obs_obj.transit_duration / 2,
                        100,
                    )
                    delta_eclipse_frame = np.linspace(
                        -obs_obj.transit_duration / 2,
                        +obs_obj.transit_duration / 2,
                        100,
                    )
                    transit = Time(obs_time) + delta_eclipse_frame
                    frame_ecl = AltAz(obstime=transit, location=location)
                    obs_ecl = obs_obj.Coordinates.coord.transform_to(frame_ecl)

                    obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
                    if any(
                        obs_altazs[obs_altazs.alt > 85 * u.deg]
                    ):  # checking if target reaches zenith angle
                        warning = 1
                    im = ax1.scatter(
                        delta_midnight,
                        obs_altazs.alt,
                        c=obs_altazs.secz.value,
                        label=obs_obj.name,
                        lw=0,
                        s=8,
                        cmap="viridis",
                        vmin=-10,
                        vmax=10,
                    )  # plot candidate
                    if eclipse["n_exposures_possible"] >= 20:
                        ax1.scatter(
                            delta_eclipse, obs_ecl.alt, color="red", lw=3, s=8
                        )  # plot transit
                else:
                    no_ecl_observable = 1

    if no_ecl_observable == 1:
        obs_altazs = obs_obj.Coordinates.coord.transform_to(frame_obs)
        im = ax1.scatter(
            delta_midnight,
            obs_altazs.alt,
            label=obs_obj.name,
            lw=0,
            s=8,
            cmap="viridis",
            vmin=-10,
            vmax=10,
        )  # Plot candidate
        # phi = np.linspace(0, np.pi, 20)
        # second_xticks = obs_altazs.az[np.int16(np.floor(500*(1+np.tanh(phi/2))))]
        # ax2.set_xlim(obs_altazs.az[0], obs_altazs.az[-1])
        # ax2.set_xticks(second_xticks)
        # ax2.set_xlabel('Azimuth Target [deg]')
    if no_list == 1:
        fig.colorbar(im).set_label("Airmass")

    fig.legend(loc="upper left")
    ax1.set_xlim(-12 * u.hour, 12 * u.hour)
    ax1.set_xticks((np.arange(13) * 2 - 12) * u.hour)
    if warning == 1:
        ax1.set_title(f"{date.date()}, WARNING: Target close to Zenith!")
    else:
        ax1.set_title(f"{date.date()}")
    ax1.set_ylim(0 * u.deg, 90 * u.deg)
    ax1.set_xlabel("Hours from EDT Midnight")
    ax1.set_ylabel("Altitude [deg]")
    plt.show()

    path = join(dirname(__file__), "../Plots")

    fig.savefig(f"{path}/{d}-single_night.eps")


##########################################################################################################

def save_csv(filename : str, data: pd.DataFrame):
    path = join(dirname(__file__), "../csv_files")
    filename = filename.split(".")[0] + ".csv"
    filename = filename.format(path=path)
    filename = realpath(filename)

    if not isinstance(data, pd.DataFrame):
        data = to_pandas(data)
    data.to_csv(filename)

    print(f"Data written to {filename}")

def save_excel(filename : str, data: pd.DataFrame):
    path = join(dirname(__file__), "../csv_files")
    filename = filename.split(".")[0] + ".xlsx"
    filename = filename.format(path=path)
    filename = realpath(filename)

    if not isinstance(data, pd.DataFrame):
        data = to_pandas(data)
    data.to_excel(filename)

    print(f"Data written to {filename}")
