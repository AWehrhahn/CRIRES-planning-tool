#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:49:49 2020

Classes for Transit_List observability

@author: jonaszbinden
"""

import datetime
import logging
import pickle
from os.path import dirname, join

import astroplan
import astropy
import astropy.units as u
import numpy as np
import pandas as pd
from astroplan import FixedTarget, Observer
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
from tqdm import tqdm

import classes_methods.Helper_fun as fun

""" Location and UTC offset Paranal """
paranal = Observer.at_site("paranal", timezone="Chile/Continental")
logger = logging.getLogger(__name__)


class Exoplanets:
    """
        Stores list of available planets for processing. The list of candidates gets requested with the file
    """

    def __init__(self):
        """
            Initialize lists to sort the planets retrieved from Nasa Exoplanet Archive according
            to the available data to compute Transits.

            data : contains all planets found in the Nasa Exoplanet Archive with names from a name-list.

            data_missing : contains all the planet names for which not all the necessary datas and infos about the planet were available.

            Parse_Planet_Nasa : contains all planet tables for which the necessary data about the planet have been found.
                                This list gets used for later processing of the planets to asses their observability.

        """
        #:pd.DataFrame: contains all planets found in the Nasa Exoplanet Archive with names from a name-list.
        self.data = []
        #:List: same as data, but in a list
        self.data_list = []
        #:pd.DataFrame: contains all the planet names for which not all the necessary datas and infos about the planet were available.
        self.data_missing = []
        #:pd.DataFrame: contains all planet tables for which the necessary data about the planet have been found. This list gets used for later processing of the planets to asses their observability.
        self.data_complete = []
        #:str: The catalog that was used as the data source
        self.catalog = None

    def load_catalog(self, catalog):
        """ 
        
            Checking if Planet can be found in Nasa Exoplanet Archive 
            
            Parameters
        	----------
        	name : str
                name of planet, loaded from the file PlanetList.csv or any other file 		
                containing all the accepted planets from the Nasa Exoplanet Archive
        
        
        """
        filenames = {
            "nexa_new": "PlanetList.csv",
            "nexa_old": "PlanetList_old.csv",
            "custom": "PlanetList_edit.csv",
        }
        # Use a standard catalog if available
        # Otherwise interpret it as a filename
        if catalog in filenames.keys():
            filename = join(dirname(__file__), "../csv_files", filenames[catalog])
        else:
            filename = catalog

        candidate_list = pd.read_csv(filename)

        self.data = candidate_list
        self.data_list = [p for _, p in candidate_list.iterrows()]
        self.catalog = catalog
        return self.data

    ##########################################################################################################

    def filter_data(self):
        """
        Checks if the planet tables have all the necessary data for later processing.

        Parameters
        ----------
        catalog: str
            T

        """
        if self.catalog == "nexa_old":

            for planet in self.data:
                flag = None
                logger.info(f"Checking Planet {planet['pl_name']} for Transit Data")
                if np.ma.is_masked(planet["pl_tranmid"][0]) is True:
                    logger.warning(
                        f"Planet {planet['pl_name'][0]} has no data about transit mid"
                    )
                    flag = True
                if np.ma.is_masked(planet["pl_orbper"][0]) is True:
                    logger.warning(
                        f"Planet {planet['pl_name'][0]} has no data about orbit period"
                    )
                    flag = True
                if np.ma.is_masked(planet["pl_trandur"][0]) is True:
                    logger.warning(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit duration\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["ra"][0]) is True:
                    logger.warning(
                        "Planet "
                        + planet["pl_name"][0]
                        + " has no data about transit duration\n"
                    )
                    flag = True
                if np.ma.is_masked(planet["dec"][0]) is True:
                    logger.warning(
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
                        logger.warning(
                            "No Sky coordinates found for " + planet["pl_name"][0]
                        )
                    except:
                        flag = True
                if not sky_coords:
                    flag = True

                if not flag:
                    self.data_complete.append(planet)
                    logger.info(
                        "Planet " + planet["pl_name"] + " added to Transit_data_avail\n"
                    )
                else:
                    self.data_missing.append(planet["pl_name"])
                    logger.info(
                        "Planet " + planet["pl_name"] + " added to data_missing\n"
                    )

        elif self.catalog == "nexa_new":
            select = (
                np.isnan(self.data["ra"])
                | np.isnan(self.data["pl_trandur"])
                | np.isnan(self.data["pl_orbper"])
                | np.isnan(self.data["pl_tranmid"])
                | np.isnan(self.data["dec"])
            )
            self.data_missing = self.data[select]
            self.data_complete = self.data[~select]
        return self.data_complete


class Nights(object):
    def __init__(self, d, max_delta_days, load_from_pickle=0):
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
        d_end = d + dt * max_delta_days
        self.max_delta_days = max_delta_days
        if load_from_pickle == 1:
            """ Check if there exist pkl files with night data for the preferred range, if not: initializes empty lists for Nights instance """
            try:
                d_str = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(
                    d_str, self.max_delta_days
                )
                nights = fun.pickled_items(filename)
                self.start = nights.__next__()
                self.end = nights.__next__()
                self.max_delta_days = nights.__next__()
                self.date = nights.__next__()
                self.night = nights.__next__()
                self.loaded = 1
            except Exception as e:
                logger.error(e)
                self.loaded = 0
                d_str = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(
                    d_str, self.max_delta_days
                )
                logger.warning(
                    "No Night data found for {}, computing nights...".format(filename)
                )

        self.start = d
        self.end = d_end
        self.max_delta_days = (d_end - d).days
        # self.date = []
        # self.night = []
        self.loaded = 0

        # list with datetime objects for the midnights in question
        self.date = [self.start + dt * k for k in range(self.max_delta_days)]

    ##########################################################################################################

    def calculate_nights_paranal(
        self, delta_midnight, observatory=paranal, write_to_pickle=0
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
            logger.info(
                "Nights loaded from file, continueing with processing planets for observability"
            )
        else:

            """ All times are in UTC respectively """
            logger.info(
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

            if write_to_pickle == 1:
                """Write Nights_paranal_table to file"""

                d = self.start
                d = (
                    d.isoformat()
                )  # start date from which the nights in paranal are calculated
                filename = "Nights_paranal_{}_{}d.pkl".format(d, self.max_delta_days)
                fun.pickle_dumper_objects(filename, self)


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

    def __init__(self, max_delta_days, planet=None):
        #:str: Name of the Planet
        self.name = None
        #:Time: Time of the middle of a Transit
        self.epoch = None
        #:Quantity: Orbital Period
        self.period = None
        #:float: Uncertainty of the orbital period
        self.period_err = None
        #:Quantity: Duration of the transit
        self.transit_duration = None
        #:float: Eccentricity of the planet orbit
        self.eccentricity = None
        #:Quantity: Effective temperature of the host star
        self.star_teff = None
        #:Quantity: J band magnitude of the host star
        self.star_jmag = None
        #:EclipsingSystem: astroplan representation of the planet system
        self.planets_eclipse = None
        #:SkyCoord: Sky Coordinates of the host star system
        self.coordinates = None
        #:int: number of eclipses occurring during max_delta_days
        self.num_eclipses = None
        #:list: targets that are observable during max_delta_days
        self.target_observable = []
        #:list: eclipses that are observable during max_delta_days
        self.eclipse_observable = []

        if planet is not None:
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
            # CAREFUL in nexa_old the pl_trandur is given in u.day
            self.transit_duration = planet["pl_trandur"] * u.hour
            self.eccentricity = planet["pl_orbeccen"]
            self.star_teff = planet["st_teff"] * u.K
            self.star_jmag = planet["sy_jmag"] * u.mag

            # time in barycentric frame
            self.planets_eclipse = astroplan.EclipsingSystem(
                primary_eclipse_time=self.epoch,
                orbital_period=self.period,
                duration=self.transit_duration,
            )

            self.num_eclipses = int(np.floor(max_delta_days / (self.period / u.day)))

            """ coordinates of the object in IRCS """
            sky_coord = SkyCoord(ra=planet["ra"], dec=planet["dec"], unit=u.deg)
            self.coordinates = FixedTarget(name=self.name, coord=sky_coord)

    ##########################################################################################################

    def observability(
        self,
        obs_time,
        nights,
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

        logger.info(f"{self.name} is getting processed")
        result_eclipse, result_target = [], []

        eclipse_times = self.planets_eclipse.next_primary_eclipse_time(
            obs_time, n_eclipses=self.num_eclipses
        )
        eclipse_times.location = paranal.location
        ltt_bary = eclipse_times.light_travel_time(self.coordinates.coord)
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
            is_observable = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.coordinates,
                times_ingress_egress=eclipse_night[:, (0, 2)],
            )
            is_observable = is_observable[0]

            if np.any(is_observable):
                eclipse_night = eclipse_night[is_observable]
                eclipse_error = eclipse_error[is_observable]

                moon_sep, moon_phase, airmass, altazs = fun.airmass_moon_sep_obj_altaz(
                    self, eclipse_night
                )

                for i in range(len(eclipse_night)):
                    eclipse_triplet = eclipse_night[i]

                    eclipse_dict = {
                        "name": self.name,
                        "obs_time": eclipse_triplet[2],
                        "obs_time_error": eclipse_error[i],
                        "is_primary_eclipse_observable": True,
                        "transit_duration": self.transit_duration,
                        "stellar_effective_temperature": self.star_teff,
                        "magnitude_j": self.star_jmag,
                    }
                    eclipse_dict.update(
                        {
                            f"eclipse_{step}": {
                                "time": eclipse_triplet[j],
                                "airmass": airmass[i][j],
                                "moon sep": moon_sep[0][i][j],
                                "moon phase": moon_phase[i][j],
                                "az": altazs[i][j].az,
                                "alt": altazs[i][j].alt,
                            }
                            for j, step in enumerate(("begin", "mid", "end"))
                        }
                    )
                    result_eclipse += [eclipse_dict]
        if check_target:
            """ Check if Nights object has attribute nights to calculate observability of target """
            if hasattr(nights, "night"):
                k = list.index(nights.date, eclipse_times)
                night = nights.night[k]
            else:
                nights.Calculate_nights_paranal(delta_midnight)
                night = nights.night[0]
            """ Check if target observable independent of Transit, can not be turned on through menu yet. """
            tar_obs = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.coordinates,
                times=night,
            )
            tar_obs = tar_obs[0].reshape((3, -1)).T
            tar_obs = np.all(tar_obs, axis=1)
            logger.info("{self.name} target is observable without any primary eclipse")
            for i in np.where(tar_obs)[0]:
                moon_sep, moon_phase, airmass, altazs = fun.airmass_moon_sep_obj_altaz(
                    self, night[i]
                )
                result_target += {
                    "Name": self.name,
                    "Effective Temperature": self.star_teff,
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


def load_eclipses_from_file(filename, max_delta_days):
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
    eclipses_list = []

    def att_identifier(att):
        if type(att) == int:
            planet.num_eclipses = att
        elif type(att) == FixedTarget:
            planet.coordinates = att
        elif type(att) == astroplan.periodic.EclipsingSystem:
            planet.planets_eclipse = att
        else:
            return att

    planet = fun.pickled_items(filename)
    att = None
    while True:
        planet = Eclipses(max_delta_days)
        try:
            if att != None:
                planet.name = att
            else:
                planet.name = next(planet)

            planet.epoch = next(planet)
            planet.period = next(planet)
            planet.period_err = next(planet)
            planet.transit_duration = next(planet)
            planet.eccentricity = next(planet)
            planet.star_teff = next(planet)
            planet.star_jmag = next(planet)
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

            planet.eclipse_observable.extend(eclipse_observable)
            planet.target_observable.extend(target_observable)
            eclipses_list.append(planet)
        except StopIteration:
            logger.info("Eclipses_List has bin loaded from {}".format(filename))
            planet.eclipse_observable.extend(eclipse_observable)
            planet.target_observable.extend(target_observable)
            eclipses_list.append(planet)
            break

    return eclipses_list


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

    def __init__(self, name, star_teff, star_jmag, coordinates=None):
        """ Initialize Eclipse instance from Nasa query_planet object """
        self.name = name
        self.star_teff = star_teff
        self.star_jmag = star_jmag

        self.target_observable = []

        try:
            self.coordinates = FixedTarget.from_name(self.name)
        except Exception:
            self.coordinates = FixedTarget(coordinates, name=self.name)

    ##########################################################################################################

    def target_observable(self, nights, constraints, delta_midnight=None):
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

        logger.info(self.name + " is getting processed")
        for date in nights.date:
            """ Check if Nights object has attribute nights to calculate observability of target """
            if hasattr(nights, "night"):
                k = list.index(nights.date, date)
                night = nights.night[k]
            else:
                nights.Calculate_nights_paranal(delta_midnight)
                night = nights.night[k]

            """ Check if target observable """
            tar_obs = astroplan.is_event_observable(
                constraints=constraints,
                observer=paranal,
                target=self.coordinates,
                times=night,
            )
            if any(tar_obs[0] == True):
                logger.info(
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
                                "Effective Temperature": self.star_teff,
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
