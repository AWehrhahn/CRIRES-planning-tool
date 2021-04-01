#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:49:49 2020

Classes for Transit_List observability

@author: jonaszbinden
"""

import math as m
from os.path import dirname, join

import astroplan
import astropy
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed, wait
from tqdm import tqdm

# from astroplan import EclipsingSystem
# import astropy.coordinates
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon

# from astropy import table
from astropy.time import Time
from astropy.visualization import astropy_mpl_style, quantity_support

plt.style.use(astropy_mpl_style)
quantity_support()
# from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
# import astroquery.open_exoplanet_catalogue as oec
import datetime
import logging
import pickle
import time

import pandas as pd
from astroplan import FixedTarget, Observer, download_IERS_A, get_IERS_A_or_workaround
from astropy.coordinates import get_sun
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

import classes_methods.Helper_fun as fun
from classes_methods import csv_file_import
from classes_methods.Helper_fun import help_fun_logger

# from threading import Thread


""" Location and UTC offset Paranal """
paranal = Observer.at_site("paranal", timezone="Chile/Continental")

##########################################################################################################


class Exoplanets:
    """
        Stores list of available planets for processing. The list of candidates gets requested with the file
    """

    @help_fun_logger
    def __init__(self):
        """
            Initialize lists to sort the planets retrieved from Nasa Exoplanet Archive according
            to the available data to compute Transits.

            Exoplanets_List_Nasa : contains all planets found in the Nasa Exoplanet Archive with names from a name-list.

            Exoplanet_not_Nasa : contains all the names for which no data was found in Nasa Exoplanet Archive

            Transit_data_missing : contains all the planet names for which not all the necessary datas and infos about the planet were available.

            Parse_Planet_Nasa : contains all planet tables for which the necessary data about the planet have been found.
                                This list gets used for later processing of the planets to asses their observability.

        """
        self.Exoplanets_List_Nasa = []
        self.Exoplanet_not_Nasa = []

        self.Transit_data_missing = []
        self.Parse_planets_Nasa = []
        self.Fail = []

    @help_fun_logger
    def Planet_finder(self, catalog):
        """ 
        
            Checking if Planet can be found in Nasa Exoplanet Archive 
            
            Parameters
        	----------
        	name : str
                name of planet, loaded from the file PlanetList.csv or any other file 		
                containing all the accepted planets from the Nasa Exoplanet Archive
        
        
        """
        # Planet_try = NasaExoplanetArchive.query_planet(name, all_columns=True)
        path = join(dirname(__file__), "../csv_files")
        filenames = {
            "nexa_new": join(path, "PlanetList.csv"),
            "nexa_old": join(path, "PlanetList_old.csv"),
            "custom": join(path, "PlanetList_edit.csv"),
        }
        filename = filenames[catalog]
        candidate_list = pd.read_csv(filename)

        self.exoplanet_list = candidate_list
        self.exoplanets_list_nasa = [p for _, p in candidate_list.iterrows()]
        self.exoplanet_not_nasa = []

    ##########################################################################################################

    @help_fun_logger
    def hasproperties(self, catalog):
        """
            Checks if the planet tables have all the necessary data for later processing.

            Returns
            -------
            None.

        """
        if catalog == "nexa_old":

            for planet in self.Exoplanets_List_Nasa:
                flag = None
                print("Checking Planet " + planet["pl_name"] + " for Transit Data")
                if np.ma.is_masked(planet["pl_tranmid"][0]) is True:
                    print(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit mid\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["pl_orbper"][0]) is True:
                    print(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about orbit period\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["pl_trandur"][0]) is True:
                    print(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit duration\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["ra"][0]) is True:
                    print(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit duration\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["dec"][0]) is True:
                    print(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit duration\n"
                    )
                    flag = True
                try:
                    sky_coords = SkyCoord.from_name(planet["pl_name"][0])
                except Exception:
                    try:
                        sky_coords = planet["sky_coord"][0]
                        print("no Sky coordinates found for " + planet["pl_name"][0])
                    except:
                        flag = True
                if not sky_coords:
                    flag = True

                if not flag:
                    self.Parse_planets_Nasa.append(planet)
                    print(
                        "Planet " + planet["pl_name"] + " added to Transit_data_avail\n"
                    )
                else:
                    self.Transit_data_missing.append(planet["pl_name"])
                    print(
                        "Planet "
                        + planet["pl_name"]
                        + " added to Transit_data_missing\n"
                    )

        elif catalog == "nexa_new":

            for planet in self.Exoplanets_List_Nasa:
                # print("Checking Planet " + planet["pl_name"] + " for Transit Data")
                if (
                    m.isnan(planet["ra"])
                    or m.isnan(planet["pl_trandur"])
                    or m.isnan(planet["pl_orbper"])
                    or m.isnan(planet["pl_tranmid"])
                    or m.isnan(planet["dec"])
                ):
                    self.Transit_data_missing.append(planet["pl_name"])
                    # print(
                    #     f"Planet {planet['pl_name']} added to Transit_data_missing"
                    # )
                else:
                    self.Parse_planets_Nasa.append(planet)
                    # print(
                    #     f"Planet {planet['pl_name']} added to Transit_data_avail"
                    # )


##########################################################################################################


class Nights(object):
    @help_fun_logger
    def __init__(self, d, Max_Delta_days, LoadFromPickle=0):
        """
            Calculates the nights at paranal for a certain start date ''d'' and end date, reached after ''Max_Delta_days''.
            Retrieves the sun coordinates for each night from astroplan to determine the night times.
    
            Parameters
            ----------
            d : datetime.date
                Start date from which the nights at paranal are computed.
    
            Max_Delta_days : int
                Time span for which the nights at paranal are computed.
    
            LoadFromPickle : int
                If ''LoadFromPickle'' = 1, checks if Night data for given time span is available.
        """
        dt = datetime.timedelta(days=1)
        d_end = d + dt * Max_Delta_days
        self.Max_Delta_days = Max_Delta_days
        if LoadFromPickle == 1:
            """ Check if there exist pkl files with night data for the preferred range, if not: initializes empty lists for Nights instance """
            try:
                d_str = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(
                    d_str, self.Max_Delta_days
                )
                nights = fun.pickled_items(filename)
                self.start = nights.__next__()
                self.end = nights.__next__()
                self.Max_Delta_days = nights.__next__()
                self.date = nights.__next__()
                self.night = nights.__next__()
                self.loaded = 1
            except Exception as e:
                print(e)
                self.loaded = 0
                d_str = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(
                    d_str, self.Max_Delta_days
                )
                print(
                    "No Night data found for {}, computing nights...".format(filename)
                )

        self.start = d
        self.end = d_end
        self.Max_Delta_days = (d_end - d).days
        # self.date = []
        # self.night = []
        self.loaded = 0

        # list with datetime objects for the midnights in question
        self.date = [self.start + dt * k for k in range(self.Max_Delta_days)]

    ##########################################################################################################

    @help_fun_logger
    def Calculate_nights_paranal(
        self, delta_midnight, observatory=paranal, WriteToPickle=0
    ):
        """
            Calculates the nights at ''observatory'', default=paranal for a certain start date and end date. Retrieves the sun coordinates
            for each night from astroplan.

            Parameters
            ----------
            delta_midnight : numpy.linspace
                array containing a grid of timesteps for which the nights datetimes should get computed.

            observatory : astroplan.Observer (optional)
                contains EarthLocation and timezone info about the observer at, default is paranal.

            WriteToPickle : int
                Object Nights gets written into a pickle file.

            Returns
            -------
            self.dates : list
                Contains datetime.date objects for each night between d and d_end.

            self.coords : list
                Contains dict with sun coordinates for each time in Nights.night.

            self.night : list
                Contains lists with nighttimes for each night. The number timesteps for each nights is defined in
                delta_midnight.

        """

        if self.loaded == 1:
            """ If the nights at paranal with startdate d could be loaded from file, yield: """
            print(
                "Nights loaded from file, continueing with processing planets for observability"
            )
        else:

            """ All times are in UTC respectively """
            print(
                "Calculating the nights of paranal from the {} until the {}".format(
                    self.start, self.end
                )
            )
            self.night = []
            midnight = datetime.time(0, 0, 0)
            for date in self.date:

                # list with datetime objects for the midnights in question
                midnight_datetime = datetime.datetime.combine(date, midnight)
                # Time object for each midnight gets created in UTC, midnight in UTC.
                midnight_at_site_UTC = observatory.datetime_to_astropy_time(
                    midnight_datetime
                )

                Night_paranal = midnight_at_site_UTC + delta_midnight  # in UTC
                # compute frame AltAz for get_sun
                frame_24_h_paranal = AltAz(
                    obstime=Night_paranal, location=observatory.location
                )
                sunaltazs_24_h = get_sun(Night_paranal).transform_to(frame_24_h_paranal)

                night = []
                for n, _ in enumerate(delta_midnight):
                    """ Calculates the night times for each night """
                    if sunaltazs_24_h[n].alt < -18 * u.deg:
                        night.append(str(sunaltazs_24_h[n].obstime.value))
                self.night.append(Time(night))

            if WriteToPickle == 1:
                """Write Nights_paranal_table to file"""

                d = self.start
                d = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(d, self.Max_Delta_days)
                fun.pickle_dumper_objects(filename, self)


##########################################################################################################


class Eclipses:
    """
        Initialization of Eclipse class. For a planet the necessary data for eclipse observation get initialized here.

        Parameters:
        -------------------
        name : string
            'pl_name': name of the planet

        epoch : astropy.time.core.Time
            'pl_tranmid': the mid time of the next transit in barycentric frame

        period : astropy.units.quantity.Quantity
            'pl_orbper': orbital period of the planet around its host star in u.day

        period_err : astropy.units.quantity.Quantity
            'pl_orbpererr1': measured error of orbital period of the planet around its host star 

        transit_duration : astropy.units.quantity.Quantity
            'pl_trandur': duration of a transit in u.hour

        Coordinates : astropy.coordinates.sky_coordinate.SkyCoord (ICRS)
            'sky_coord': right ascension and azimuth of host star in degrees

        eccentricity : float
            'pl_eccen': eccentricity of the orbit of the planet around its host star

        star_Teff : astropy.units.quantity.Quantity
            'st_Teff': Effective temperature of the host star in Kelvin

        star_jmag : float
            'st_j': Magnitude of the host star in the J-band
            
        Max_Delta_days : int
            Days for which the eclipses get computed

        -------------------

        Other parameters get initialized empty and get assigned later.
        More parameters can be added manually. The parameters all come from the
        'NasaExoplanetArchive.query_planet(name, all_columns=True)' function from
        astroquery. For the filtering of the Exoplanet candidates refer to
        'Request_Table_NasaExoplanetArchive.py' and the file used to filter the Archive:
        'Nasa_Archive_Selection.txt', in this file you may find information to look up keywords
        that can be used to add additional information to a instance of Eclipses.

        -------------------

        Furthermore initializing an instance of this class calls astroplan.EclipsingSystem and creates an instance of
        EclipsingSystem, using the already initialized data. The instance is stored under self.Planets_eclipse.
        Additionally the number of eclipses possible in the evaluated timespan is computed and stored in self.num_eclipses.

    """

    @help_fun_logger
    def __init__(self, Max_Delta_days, planet=None):
        try:
            if planet == None:
                empty_class = True
            else:
                empty_class = False

        except Exception:

            if planet.all() == None:
                empty_class = True
            else:
                empty_class = False

        if empty_class:

            """ Create Empty class to load data from file """

            self.name = None
            self.epoch = None
            self.period = None
            self.period_err = None
            self.transit_duration = None
            self.eccentricity = None
            self.star_Teff = None
            self.star_jmag = None
            # self.pl_mass = None
            self.Planets_eclipse = None
            self.num_eclipses = None

            self.target_observable = []
            self.eclipse_observable = []

        else:

            """ Initialize Eclipse instance from Nasa query_planet object """

            self.name = planet["pl_name"]
            self.epoch = Time(
                planet["pl_tranmid"],
                format="jd",
                scale="utc",
                location=paranal.location,
            )
            self.period = planet["pl_orbper"] * u.day
            self.period_err = planet["pl_orbpererr1"]
            self.transit_duration = (
                planet["pl_trandur"] * u.hour
            )  # CAREFUL in nexa_old the pl_trandur is given in u.day
            self.eccentricity = planet["pl_orbeccen"]
            self.star_Teff = planet["st_teff"]
            self.star_jmag = planet["sy_jmag"]
            # self.pl_mass = planet['pl_bmassj']

            self.target_observable = []
            self.eclipse_observable = []

            """
                Planet_next_eclipse : Class object with which the future transits of the planet can be calculated.
                                      For documentation review astroplan.EclipsingSystem documentation.
    
                WARNING:
                There are currently two major caveats in the implementation of
                ''EclipsingSystem''. The secondary eclipse time approximation is
                only accurate when the orbital eccentricity is small, and the eclipse
                times are computed without any barycentric corrections. The current
                implementation should only be used forapproximate mid-eclipse times for
                low eccentricity orbits, with event durations longer than the
                barycentric correction error (<=16 minutes). 
                
                This is corrected for with the barycentric correction, which can be found in the classmethod Observability
    
                From EclipsingSystem.__doc__
    
            """

            Planet_next_eclipse = astroplan.EclipsingSystem(
                primary_eclipse_time=self.epoch,
                orbital_period=self.period,
                duration=self.transit_duration,
            )

            self.Planets_eclipse = Planet_next_eclipse  # time in barycentric frame

            # number of eclipses occurring during Max_Delta_days.

            self.num_eclipses = int(np.floor(Max_Delta_days / (self.period / u.day)))

            """ coordinates of the object in IRCS """
            try:
                self.Coordinates = FixedTarget.from_name(self.name)
            except Exception:
                sky_coord = SkyCoord(ra=planet["ra"] * u.deg, dec=planet["dec"] * u.deg)
                self.Coordinates = FixedTarget(name=self.name, coord=sky_coord)

    ##########################################################################################################

    @help_fun_logger
    def Observability(
        self,
        obs_time,
        Nights,
        constraints,
        check_eclipse,
        check_target=False,
        delta_midnight=None,
    ):
        """
            Calculates if the Transit and the target are observable for each date during the given timespan in ''Nights'' under
            the given ''constraints'' and writes it as dict objects into ''~self.eclipse_observable'' or ''~self.target_observable''.

            Parameters
            ----------
            obs_time : astropy.time.Time
                Contains the datetime as Time format after which the possible observations should be found.

            Nights : class Nights
                Containing night data of paranal, see Nights documentation.

            constraints : class astroplan.Constraint
                Constraints under which the observational events should get constrained.

            check_eclipse : bool
                If ''check_eclipse'' =  1, checks if transits/eclipses are observable.

            check_target : bool, optional
                If ''check_target'' = 1, checks if target is observable during the given nights. The default is 0.

            delta_midnight : numpy.linspace, optional
                array containing a grid of timesteps for which the nights datetimes should get computed. Default is None

        """

        print(self.name + " is getting processed")
        result_eclipse, result_target = [], []

        eclipse_times = self.Planets_eclipse.next_primary_eclipse_time(
            obs_time, n_eclipses=self.num_eclipses
        )
        eclipse_times.location = paranal.location
        ltt_bary = eclipse_times.light_travel_time(self.Coordinates.coord)
        eclipse_times -= ltt_bary

        n = np.arange(len(eclipse_times))

        if check_eclipse:
            eclipse_begin = eclipse_times - self.transit_duration / 2
            eclipse_end = eclipse_times + self.transit_duration / 2
            eclipse_error = (n + 1) * self.period_err

            # Beginning and end of transit
            eclipse_night = Time(
                np.concatenate([eclipse_begin, eclipse_times, eclipse_end])
            )
            # Needs to be reshaped to (N, 2) shape
            eclipse_night = eclipse_night.reshape((3, -1)).T
            eclipse_obs = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.Coordinates,
                times_ingress_egress=eclipse_night[:, (0, 2)],
            )
            eclipse_obs = eclipse_obs[0]
            for i in np.where(eclipse_obs)[0]:
                print(f"{self.name} total Eclipse is observable")
                eclipse_triplet = eclipse_night[i]
                moon_sep, moon_phase, airmass, altazs = fun.airmass_moon_sep_obj_altaz(
                    self, eclipse_triplet
                )

                eclipse_dict = {
                        "name": self.name,
                        "obs_time": eclipse_times[i],
                        "obs_time_error": eclipse_error[i],
                        "is_primary_eclipse_observable": True,
                        "transit_duration": self.transit_duration,
                        "stellar_effective_temperature": self.star_Teff,
                        "magnitude_j": self.star_jmag,
                    }
                eclipse_dict.update({
                    f"eclipse_{step}": {
                        "time": eclipse_triplet[j],
                        "airmass": airmass[j],
                        "moon sep": moon_sep[0][j],
                        "moon phase": moon_phase[j],
                        "az": altazs[j].az,
                        "alt": altazs[j].alt,
                    }
                    for j, step in enumerate(("begin", "mid", "end"))
                })
                result_eclipse += [eclipse_dict]
        if check_target:
            """ Check if Nights object has attribute nights to calculate observability of target """
            if hasattr(Nights, "night"):
                k = list.index(Nights.date, eclipse_times)
                night = Nights.night[k]
            else:
                Nights.Calculate_nights_paranal(delta_midnight)
                night = Nights.night[0]
            """ Check if target observable independent of Transit, can not be turned on through menu yet. """
            tar_obs = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.Coordinates,
                times=night,
            )
            tar_obs = tar_obs[0].reshape((3, -1)).T
            tar_obs = np.all(tar_obs, axis=1)
            print("{self.name} Target is observable without any primary eclipse")
            for i in np.where(tar_obs)[0]:
                moon_sep, moon_phase, airmass, altazs = fun.airmass_moon_sep_obj_altaz(
                    self, night[i]
                )
                result_target += {
                    "Name": self.name,
                    "Effective Temperature": self.star_Teff,
                    "J-magnitude": self.star_jmag,
                    "Object w. o. primary eclipse observable?": True,
                    "Obs Data": {
                        "time": night[i],
                        "airmass": airmass,
                        "moon sep": moon_sep[0],
                        "moon phase": moon_phase,
                        "az": altazs.az,
                        "alt": altazs.alt,
                    },
                }

        self.eclipse_observable = result_eclipse
        self.target_observable = result_target

        return self.eclipse_observable, self.target_observable


##########################################################################################################


@help_fun_logger
def load_Eclipses_from_file(filename, Max_Delta_days):
    """
        Loads Eclipses class objects from pickled file with ''filename''.

        Parameters
        ----------
        filename : str
            Name of file containing the pickled Eclipses class object data.

        Max_Delta_days : int
            Days for which the eclipses got computed, necessary to initialize Eclipses class.

        Returns
        -------
        Eclipses_List : list
            Contains all Eclipses class objects that have been loaded from the file.

    """
    Eclipses_List = []

    def att_identifier(att):
        if type(att) == int:
            Planet.num_eclipses = att
        elif type(att) == FixedTarget:
            Planet.Coordinates = att
        elif type(att) == astroplan.periodic.EclipsingSystem:
            Planet.Planets_eclipse = att
        else:
            return att

    planet = fun.pickled_items(filename)
    att = None
    while True:
        Planet = Eclipses(Max_Delta_days)
        try:
            if att != None:
                print(att)
                Planet.name = att
            else:
                Planet.name = next(planet)

            Planet.epoch = next(planet)
            Planet.period = next(planet)
            Planet.period_err = next(planet)
            Planet.transit_duration = next(planet)
            Planet.eccentricity = next(planet)
            Planet.star_Teff = next(planet)
            Planet.star_jmag = next(planet)
            # Planet.pl_mass = next(planet)
            att = None
            while att == None:
                att = next(planet)
                att = att_identifier(att)

            target_observable = att

            eclipse_observable = next(planet)

            att = None
            while att == None:
                att = next(planet)
                att = att_identifier(att)

            Planet.eclipse_observable.extend(eclipse_observable)
            Planet.target_observable.extend(target_observable)
            Eclipses_List.append(Planet)
        except StopIteration:
            print("Eclipses_List has been loaded from {}".format(filename))
            logging.info("Eclipses_List has bin loaded from {}".format(filename))
            Planet.eclipse_observable.extend(eclipse_observable)
            Planet.target_observable.extend(target_observable)
            Eclipses_List.append(Planet)
            break

    return Eclipses_List


##########################################################################################################


class Targets:
    """
        Initialization of Target class. For a star or other object the necessary data for observations get initialized here.
        Use NasaExoplanetArchive.query_star to initialize or do manually:

        target = 
        name = target['st_name']
        star_Teff = target['st_teff']
        star_jmag = target['st_j']

        Parameters:
        -------------------
        name : string
            'name': name of the object

        Coordinates : astropy.coordinates.sky_coordinate.FixedTarget (ICRS)
            'sky_coord': right ascension and azimuth of host star in degrees

        star_Teff : astropy.units.quantity.Quantity
            'st_Teff': Effective temperature of the host star in u.K (Kelvin)

        star_jmag : float
            'st_j': Magnitude of the host star in the J-band

        -------------------

        Other parameters get initialized as empty lists and get assigned later.
        More parameters can be added manually. The parameters all come from the
        'NasaExoplanetArchive.query_star(name, all_columns=True)' function from
        astroquery.

        -------------------

    """

    @help_fun_logger
    def __init__(self, name, star_Teff, star_jmag, Coordinates=None):
        """ Initialize Eclipse instance from Nasa query_planet object """
        self.name = name
        self.star_Teff = star_Teff
        self.star_jmag = star_jmag

        self.target_observable = []

        try:
            self.Coordinates = FixedTarget.from_name(self.name)
        except Exception:
            self.Coordinates = FixedTarget(Coordinates, name=self.name)

    ##########################################################################################################

    def target_observable(self, Nights, constraints, delta_midnight=None):
        """
            Calculates for which times during the time span of Nights, the target is observable under the given constraints.
            LATER : Could include plotting of target observability.

            Parameters
            ----------
            Nights : class
                Nights at Paranal for which to compute if the target is observable.
            constraints : list
                list of Astroplan constraints to constrain the observability.
            delta_midnight : numpy.linspace, Obtional
                grid of timesteps within 24 hours for which the observation should be calculated.

        """
        if delta_midnight == None:
            # defines number of timesteps per 24 hours
            delta_midnight = np.linspace(-12, 12, 1000) * u.hour

        print(self.name + " is getting processed")
        for date in Nights.date:
            """ Check if Nights object has attribute nights to calculate observability of target """
            if hasattr(Nights, "night"):
                k = list.index(Nights.date, date)
                night = Nights.night[k]
            else:
                Nights.Calculate_nights_paranal(delta_midnight)
                night = Nights.night[k]

            """ Check if target observable """
            tar_obs = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.Coordinates,
                times=night,
            )
            if any(tar_obs[0] == True):
                print(
                    "{} Target is observable without any primary eclipse".format(
                        self.name
                    )
                )
                for n, tar in enumerate(tar_obs[0]):
                    if tar == True:
                        (
                            moon_target_sep,
                            moon_phase,
                            airmass,
                            obs_altazs,
                        ) = fun.airmass_moon_sep_obj_altaz(self, night[n])
                        self.target_observable.append(
                            {
                                "Name": self.name,
                                "Effective Temperature": self.star_Teff,
                                "J-magnitude": self.star_jmag,
                                "Object observable?": tar,
                                "Obs Data": {
                                    "time": night[n],
                                    "airmass": airmass,
                                    "moon sep": moon_target_sep[0],
                                    "moon phase": moon_phase,
                                    "az": obs_altazs.az,
                                    "alt": obs_altazs.alt,
                                },
                            }
                        )


##########################################################################################################
