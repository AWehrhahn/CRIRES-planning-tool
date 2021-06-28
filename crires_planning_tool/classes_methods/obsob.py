import pandas as pd
import numpy as np
import datetime as dt
from astropy.time import Time
from astropy import coordinates
from astropy import units as u
import astroplan as ap
# from astropy.constants import h, c, kb
import wget
import os.path
from tqdm import tqdm

def jdcnv(date):
    return Time(date).jd

def eq2hor(date, sky_coords, earth_location):
    altaz = coordinates.AltAz(obstime=date, location=earth_location)
    return sky_coords.transform_to(altaz)

def helcorr(skycoord, earth_location, obstime):
    return skycoord.radial_velocity_correction("heliocentric", obstime=obstime, location=earth_location)

def planck(teff, wave):
    h = constants.h
    c = constants.c
    kb = constants.kB
    B = 2 * np.pi * h * c**2 / wave**5 * (np.exp((h*c)/(wave*kB*teff)) - 1)**(-1)
    return B

def ComputeSNR(m, t, airm, method="Crires", mband="K"):
    """
    Computes SNR of stellar spectrum from magnitude, exposure time (in seconds) and airmass
    """
    # SNR with airmass = 1

    # This is from old Carmenes documentation, factor 1.1774 so that it agrees better with Ansgars result
    if method == "Carmenes":
        if mband == "J":
            SNR_1 = 1.1774 * 100 / np.sqrt(40 * 10 ** ((m - 4.2) / 2.5)) * np.sqrt(t)
        else:
            print("Use Jmag for calculation of Carmenes SNR.")
            SNR_1 = np.nan
    elif method == "Crires":
        if mband == "K":
            SNR_1 = 449.4241 * np.sqrt(10 ** (-m / 2.5)) * np.sqrt(t) - 6.3144
        else:
            print("Use Kmag for calculation of Crires SNR.")
            SNR_1 = np.nan

    else:
        print("Method not recognized. Use Crires or Carmenes.")
        SNR_1 = np.nan

    # Scale to airmass = airm
    extcof = 0.05  # extinction coefficient, see Lombardi et al., 2018
    SNR_airm = SNR_1 * 10 ** (extcof / 5 * (1 - airm))

    return SNR_airm

def snr_estimate_nexposures(eclipse, planet):
    mag = planet.star_jmag.to_value(u.mag)
    exptime = planet.transit_duration.to_value(u.second)
    airmass = [eclipse[f"eclipse_{s}"]["airmass"].to_value(1) for s in ["begin", "mid", "end"]]
    airmass = np.array(airmass)

    snr_data = ComputeSNR(mag, exptime, airmass)

    median_snr = np.median(snr_data)
    min_snr = np.min(snr_data)
    max_snr = np.max(snr_data)

    # we don't estimate the number of exposures here
    num_exp_possible = 1
    # estimates the number of exposures possible according to the transit duration and the maximum exposure time
    eclipse["n_exposures_possible"] = num_exp_possible
    eclipse["snr_median"] = median_snr
    eclipse["snr_minimum"] = min_snr
    eclipse["snr_maximum"] = max_snr
    eclipse["average_exposure_time"] = exptime
    return eclipse

def estimate_snr(eclipses_list, minimum_snr=100):
    for i, planet in tqdm(enumerate(eclipses_list), total=len(eclipses_list), desc="Planets"):
        for j, eclipse in tqdm(enumerate(planet.eclipse_observable), total=len(planet.eclipse_observable), leave=False, desc="Transits"):
            eclipses_list[i].eclipse_observable[j] = snr_estimate_nexposures(eclipse, planet)
    return eclipses_list

def GetPlanetDataNexa(planet):
    """
    Loads data for given planet from nexa file and returns dictionary with all information.
    Downloads nexa catalog if file at filepath does not exist.
    """
    filepath = "../NexaCatalog.csv"

    if os.path.exists(filepath) == False:
        # download data if no file exists
        DownloadDataNexa(filepath)
    nexaData = pd.read_csv(filepath)

    pl_index = np.where(nexaData["pl_name"] == planet)[0]

    if pl_index.size == 0:
        # print(planet + ' not in the Nexa database')
        return None
    else:
        pl_index = pl_index[0]

    planetData = {
        "ra": nexaData["ra"][pl_index],
        "dec": nexaData["dec"][pl_index],
        "T0": nexaData["pl_tranmid"][pl_index],
        "orbPer": nexaData["pl_orbper"][pl_index],
        "orbInc": nexaData["pl_orbincl"][pl_index],
        "orbEcc": nexaData["pl_orbeccen"][pl_index],
        "SMA": nexaData["pl_orbsmax"][pl_index],
        "RpJ": nexaData["pl_radj"][pl_index],
        "RsSun": nexaData["st_rad"][pl_index],
        "MpJ": nexaData["pl_bmassj"][pl_index],
        "MsSun": nexaData["st_mass"][pl_index],
        "Tdur": nexaData["pl_trandur"][pl_index] / 24,  # convert from hours to days
        "plName": nexaData["pl_name"][pl_index],
        "PerErr": nexaData["pl_orbpererr1"][pl_index],
        "TransitDepth": nexaData["pl_trandep"][pl_index],
        "Teq": nexaData["pl_eqt"][pl_index],
        "Teff": nexaData["st_teff"][pl_index],
        "Vmag": nexaData["sy_vmag"][pl_index],
        "Hmag": nexaData["sy_hmag"][pl_index],
        "Jmag": nexaData["sy_jmag"][pl_index],
        "Kmag": nexaData["sy_kmag"][pl_index],
    }

    return planetData


def DownloadDataNexa(filepath):
    """
    Downloads the columns in properties for all planets in Nexa to a file in filepath
    """
    properties = [
        "hostname",
        "pl_letter",
        "pl_name",
        "pl_orbper",
        "pl_orbsmax",
        "pl_radj",
        "pl_bmassj",
        "ra",
        "dec",
        "pl_orbincl",
        "pl_orbeccen",
        "pl_orbpererr1",
        "sy_vmag",
        "sy_hmag",
        "sy_jmag",
        "sy_kmag",
        "st_teff",
        "st_rad",
        "st_mass",
        "pl_eqt",
        "pl_trandep",
        "pl_trandur",
        "pl_tranmid",
    ]

    urlRoot = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
    select = "select+"
    select = select + ",".join(properties)
    select = select[:7] + select[8:]
    table = "+from+pscomppars"
    outformat = "&format=csv"

    url = "".join((urlRoot, select, table, outformat))

    wget.download(url, out=filepath)


def ListAllPlanetsNexa(maxKmag=12):
    """
    returns list of planets of Nexa with Vmag< maxVmag
    """
    filepath = "../NexaCatalog.csv"
    if os.path.exists(filepath) == False:
        # download data if no file exists
        DownloadDataNexa(filepath)
    nexaData = pd.read_csv(filepath)

    planetlist = []
    for planetnumber in range(len(nexaData["pl_name"])):
        if nexaData["sy_kmag"][planetnumber] < maxKmag:
            planetlist.append(nexaData["pl_name"][planetnumber])

    return planetlist


def GetPlanetDataCustom(planet):
    """
    Loads data for given planet from file with custom table
    """
    filepath = "../CustomCatalog.csv"

    if os.path.exists(filepath) == False:
        print("Custom Catalog not found!")
    customData = pd.read_csv(filepath)

    pl_index = np.where(customData["pl_name"] == planet)[0]

    if pl_index.size == 0:
        print(planet + " not in the custom database")
        return None
    else:
        pl_index = pl_index[0]

    planetData = {
        "ra": customData["ra"][pl_index],
        "dec": customData["dec"][pl_index],
        "T0": customData["pl_tranmid"][pl_index],
        "orbPer": customData["pl_orbper"][pl_index],
        "orbInc": customData["pl_orbincl"][pl_index],
        "orbEcc": customData["pl_orbecl"][pl_index],
        "SMA": customData["pl_orbsmax"][pl_index],
        "RpJ": customData["pl_radj"][pl_index],
        "RsSun": customData["st_rad"][pl_index],
        "MpJ": customData["pl_bmassj"][pl_index],
        "MsSun": customData["st_mass"][pl_index],
        "Tdur": customData["pl_trandur"][pl_index] / 24,  # convert from hours to days
        "plName": customData["pl_name"][pl_index],
        "PerErr": customData["pl_orbpererr1"][pl_index],
        "TransitDepth": customData["pl_trandep"][pl_index],
        "Teq": customData["pl_eqt"][pl_index],
        "Teff": customData["st_teff"][pl_index],
        "Vmag": customData["sy_vmag"][pl_index],
        "Hmag": customData["sy_hmag"][pl_index],
        "Jmag": customData["sy_jmag"][pl_index],
        "Kmag": customData["sy_kmag"][pl_index],
    }

    return planetData


def CheckRequiredParameters(planetdata):
    """
    Returns False if any of the required parameters is nan, otherwise returns True
    """
    required_parameters = [
        "ra",
        "dec",
        "T0",
        "orbPer",
        "orbInc",
        "SMA",
        "RpJ",
        "RsSun",
        "Tdur",
        "Jmag",
    ]

    for parameter in required_parameters:
        if pd.isna(planetdata[parameter]):
            return False

    return True


def TransitInformation(
    planet,
    catalog,
    observatory,
    d_start,
    d_end,
    observation_puffer=1 / 24,
    min_required_altitude=0,
    max_airm_good=2,
    max_sunalt_good=-20,
    SNR_method="Crires",
    mband="K",
    verbose=True,
):
    """
    For a given planet, returns a dataframe with details about each transit of the system
    The following is an outdated description...
    Parameters
    ----------
    planet : string
        Name of the system (as given in the respective Catalog)
    catalog: string
        Name of the catalog to be used. Choices: TEP, Nexa
    observation_puffer : float
        Time before and after transit that is observed (in days)
    min_required_altitude: float
        minimal altitude of the star during a transit such that the transit is considered as an option (degrees)
    d_start: datetime object
        start date of the observation window
    d_end: datetime object
        end date of the observation window    
    max_airm_good: float
        maximal airmass considered as 'good' observation conditions
    max_sunalt_good: float
        maximal solar altitude considered as 'good' observation conditions    
    
    Returns
    -------
    PlT: pd.Dataframe
        Dataframe containing these entries:
            'Tmid', 'Tmid_err', 'Obs_start', 'Obs_end', 'Trans_start', 'Trans_end', (all in HJD)
            'Sunalt_start','Sunalt_mid','Sunalt_end', (solar altitude at begin, center, and end of the observation)
            'Airm_start', 'Airm_mid', 'Airm_end', (airmass at begin, center, and end of the observation)
            'Moon_dist', (angular distance between target and moon)
            'v_bary', (barycentric velocity at center of transit)
            'GoodCondBefore', 'GoodCondDuring', 'GoodCondAfter', (time spent in region of 'good' obs. conditions (days), accurate to ~10min)
            'system', (name of the system)
            'transit_depth', (transit depth as given in the catalog)
            'stellar_SNR', (SNR of stellar spectrum computed from T_eff, Vmag and exposure time)
            'stellar_SNR_5min', (SNR for a fixed exposure time of 5 minutes)
            'stellar_SNR_transdur', (SNR for an exposure time equal to the transit duration)
            'deltav_5min', (shift in projected planet velocity over the timespan of 5 minutes [km/s])
            'H', (scale height of the planetary atmosphere in meters)
            'Delta_d', (change in transit depth due to existence of atmosphere)
                                          
    """

    if catalog == "TEP":
        physprop, hompar, hommes, obsplan = LoadTEP()
        planetData = GetPlanetDataTEP(planet)
    elif catalog == "Nexa":
        planetData = GetPlanetDataNexa(planet)
        if planetData == None:
            if verbose:
                print("Planet not in Nexa, trying to search in custom catalog.")
            planetData = GetPlanetDataCustom(planet)
    elif catalog == "Custom":
        planetData = GetPlanetDataCustom(planet)
    else:
        print("Requested catalog not recognized. Choices are TEP and Nexa")
        return None, None

    if planetData == None:
        if verbose:
            print(planet + " not contained in catalog!")
        return None, None

    if not CheckRequiredParameters(planetData):
        if verbose:
            print("At least one of the required parameters is missing!")
        return None, None

    jd_start = jdcnv(d_start)
    jd_end = jdcnv(d_end)

    ###
    ###First create dataframe with all transits
    ###

    t_min = jd_start
    t_max = jd_end

    T0 = planetData["T0"]
    period = planetData["orbPer"]
    period_err = planetData["PerErr"]
    Tdur = planetData["Tdur"]
    ra = planetData["ra"]
    dec = planetData["dec"]

    observatory_data = coordinates.EarthLocation.of_site(observatory)
    lon = observatory_data["longitude"]
    lat = observatory_data["latitude"]
    alt = observatory_data["altitude"]

    trnum_start = np.floor((t_min - T0) / period)
    trnum_end = np.ceil((t_max - T0) / period)
    # Relevant transit epochs
    tr = np.arange(trnum_start, trnum_end, 1)

    # Get list of all relevant times (start, mid and end of observation)
    t_list = []
    for epoch in tr:
        Tmid = T0 + float(epoch) * period
        T_before = Tmid - Tdur / 2 - observation_puffer
        T_after = Tmid + Tdur / 2 + observation_puffer

        if (Tmid < t_min) or (Tmid > t_max):
            # This may happen because the transit may occur in the first
            # relevant epoch but still before tmin. Likewise for tmax.
            continue

        t_list.extend([T_before, Tmid, T_after])

    t_list = np.array(t_list)

    altaz = eq2hor(
        t_list,
        (ra, dec),
        observatory_data
    )
    altitude = altaz[0]

    # remove when planet not visible from observation site
    alt_filtered = []
    t_filtered = []

    nan_list = [np.nan, np.nan, np.nan]
    for i in range(int(altitude.size / 3)):
        minalt = np.where(altitude[i * 3 : i * 3 + 3] >= min_required_altitude)[0]

        if len(minalt) == 3:
            # target visible
            alt_filtered.extend(altitude[i * 3 : i * 3 + 3])
            t_filtered.extend(t_list[i * 3 : i * 3 + 3])

    alt_filtered = np.array(alt_filtered)
    t_filtered = np.array(t_filtered)

    total_nights = int(t_filtered.size / 3)

    if total_nights == 0:
        if verbose:
            print("No transit found for " + planet)
        return None, None

    # calculate sun altitude and airmass at each time
    notnan = pd.notnull(t_filtered)
    
    sunpos_radec = coordinates.get_sun(t_filtered[notnan])
    sunpos_altaz = eq2hor(
        t_filtered[notnan],
        (sunpos_radec[1][0],
        sunpos_radec[2][0]),
        observatory_data
    )
    sunalt = np.ones(t_filtered.shape) * np.nan
    sunalt[notnan] = sunpos_altaz[0]

    airm = coordinates.AltAz(az=90 - alt_filtered)

    # calculate distance to moon (only for midpoints before and after)
    midpoints = np.array(
        [t_filtered[i * 3 + 1] for i in range(int(t_filtered.size / 3))]
    )

    mpos = coordinates.get_moon(midpoints)
    moondist = mpos[0].seperation(ra, dec)

    ##Get the date of each night
    t_jd_mid = np.array([t_filtered[i * 3 + 1] for i in range(total_nights)])

    # rough determination of timezone
    timezone = ((int(lon / 15) + 12) % 24) - 12

    t_jd_mid += timezone / 24

    # Subtract 1 day if night starts at day before (if time is before 8am)
    t_jd_mid[t_jd_mid % 1 > 0.5] -= 1
    dates_complete = Time(t_jd_mid, format="jd").iso
    night = [date[:10] for date in dates_complete]

    # Calculate barycentric velocity
    baryc_vel = []
    for tmid in [t_filtered[i * 3 + 1] for i in range(total_nights)]:
        baryc_vel.append(helcorr((ra, dec),  (lon, lat, alt), tmid)[0])
    data_dict = {
        "System": [planet for i in range(total_nights)],
        "Night": night,
        "Tmid": [t_filtered[i * 3 + 1] for i in range(total_nights)],
        "Tmid_err": [
            (t_filtered[i * 3 + 1] - T0) / period * period_err
            for i in range(total_nights)
        ],
        "Obs_start": [t_filtered[i * 3] for i in range(total_nights)],
        "Obs_end": [t_filtered[i * 3 + 2] for i in range(total_nights)],
        "Trans_start": [
            t_filtered[i * 3] + observation_puffer for i in range(total_nights)
        ],
        "Trans_end": [
            t_filtered[i * 3 + 2] - observation_puffer for i in range(total_nights)
        ],
        "Sunalt_start": [sunalt[i * 3] for i in range(total_nights)],
        "Sunalt_mid": [sunalt[i * 3 + 1] for i in range(total_nights)],
        "Sunalt_end": [sunalt[i * 3 + 2] for i in range(total_nights)],
        "Airm_start": np.array([airm[i * 3] for i in range(total_nights)]),
        "Airm_mid": np.array([airm[i * 3 + 1] for i in range(total_nights)]),
        "Airm_end": np.array([airm[i * 3 + 2] for i in range(total_nights)]),
        "Moon_dist": moondist,
        "V_bary": baryc_vel,
    }

    data_dict["GoodCond"] = (
        (np.array(data_dict["Airm_start"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_start"]) < max_sunalt_good)
        & (np.array(data_dict["Airm_end"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_end"]) < max_sunalt_good)
    )

    ###Add system data that is the same for all transits of a planet
    system_dict = {}

    mag = mband + "mag"
    # system name
    system_dict["system"] = planet
    system_dict["Mstar"] = planetData["MsSun"]
    system_dict["Rstar"] = planetData["RsSun"]
    system_dict["Mpl"] = planetData["MpJ"]
    system_dict["Rpl"] = planetData["RpJ"]
    system_dict["Period"] = planetData["orbPer"]
    system_dict["Teff"] = planetData["Teff"]
    system_dict[mag] = planetData[mag]

    # transit depth
    transit_depth = planetData["TransitDepth"]
    if np.isnan(transit_depth):
        # calculate depth if not in catalog
        rjup2rsun = 0.1005
        transit_depth = (planetData["RpJ"] * rjup2rsun / planetData["RsSun"]) ** 2

    system_dict["transit_depth"] = transit_depth

    system_dict["transit_dur"] = data_dict["Trans_end"][0] - data_dict["Trans_start"][0]

    # stellar SNR
    # here units are km and s everywhere
    # first calculate maximum exposure so that lines do not shift between pixels
    pixel_size = 1  # pixel size in km/s
    K_p = (
        2
        * np.pi
        * (planetData["SMA"] * 1.496e8)
        / (planetData["orbPer"] * 24 * 60 * 60)
    )
    # Acceleration at time of transit is almost constant, the maximum of the acceleration sine curve
    acceleration = K_p * 2 * np.pi / (planetData["orbPer"] * 24 * 60 * 60)
    exposure_time = pixel_size / acceleration

    SNR = ComputeSNR(
        planetData[mag], exposure_time, data_dict["Airm_mid"].min(), SNR_method, mband
    )
    system_dict["SNR_nosmear"] = SNR
    system_dict["t_nosmear"] = exposure_time
    system_dict["acc"] = acceleration
    system_dict["K_p"] = K_p

    # SNR and planet vel shift for fixed exposure time of 5 minutes
    SNR_5min = ComputeSNR(
        planetData[mag], 5 * 60, data_dict["Airm_mid"].min(), SNR_method, mband
    )
    system_dict["SNR_5min"] = SNR_5min
    # deltav_5min = acceleration * 5 * 60
    # data_dict['deltav_5min'] = ones_array * deltav_5min

    # SNR with exposure time of transit duration
    transit_dur = (
        (data_dict["Trans_end"][0] - data_dict["Trans_start"][0]) * 24 * 60 * 60
    )
    SNR_transdur = ComputeSNR(
        planetData[mag], transit_dur, data_dict["Airm_mid"].min(), SNR_method, mband
    )
    system_dict["SNR_transdur"] = SNR_transdur

    # scale height
    grav_constant = 6.67430 * 10 ** -11
    k_b = 1.380649 * 10 ** -23
    m_H = 1.6735575 * 10 ** -27  # mass of hydrogen atom [kg]
    surface_gravity = (
        grav_constant
        * planetData["MpJ"]
        * 1.898
        * 10 ** 27
        / (planetData["RpJ"] * 69911000) ** 2
    )

    if planetData["RpJ"] * 11.2095 < 1.5:
        # water atmosphere for small planets with R < 1.5 R_earth
        mean_mol_weight = 18
    else:
        # for hot jupiter atmosphere
        mean_mol_weight = 2.3

    Teq = planetData["Teq"]
    if np.isnan(planetData["Teq"]):
        # if no value for Teq, calculate it using formula from Kempton 2018
        Teq = (
            planetData["Teff"]
            * np.sqrt(planetData["RsSun"] * 696340 / (planetData["SMA"] * 1.496e8))
            * 0.25 ** 0.25
        )
    H = k_b * Teq / (mean_mol_weight * m_H * surface_gravity)
    system_dict["Teq"] = Teq
    system_dict["H"] = H

    # Delta_d = change in transit depth bc of atmosphere
    Delta_d = (
        2 * planetData["RpJ"] * 69911000 * H / (planetData["RsSun"] * 696340000) ** 2
    )
    system_dict["Delta_d"] = Delta_d

    # number of transits with good conditions
    moondist_okay_mask = np.array(data_dict["Moon_dist"]) > 30
    goodtransit_mask = data_dict["GoodCond"] & moondist_okay_mask

    system_dict["N_GoodTransit"] = np.where(goodtransit_mask)[0].size

    # number of transits required so that a SNR of SNR_desired is achieved
    SNR_desired = 1000
    t_required = (SNR_desired / system_dict["SNR_5min"]) ** 2 * 5 * 60
    N_Trans_SNR1000 = t_required / (system_dict["transit_dur"] * 24 * 3600)
    system_dict["N_Trans_SNR1000"] = N_Trans_SNR1000

    # number of transits required for detection (scaled by value of HD 18 bc 1 transit was enough for that target)
    N_trans_req = (1 / (SNR_transdur * Delta_d)) ** 2 / 208.264983
    system_dict["N_trans_req"] = N_trans_req

    # number of good transits in a certain time interval
    timeinterval_start = dt.datetime(2021, 1, 1, 00)
    timeinterval_end = dt.datetime(2021, 8, 1, 00)

    jd_timeinterval_start = jdcnv(timeinterval_start)
    jd_timeinterval_end = jdcnv(timeinterval_end)

    inInterval = (np.array(data_dict["Tmid"]) > jd_timeinterval_start) & (
        np.array(data_dict["Tmid"]) < jd_timeinterval_end
    )
    system_dict["N_GoodTransit_Semester"] = np.where(goodtransit_mask & inInterval)[
        0
    ].size

    # calculate score to rank the systems
    system_dict["systemscore"] = (
        Delta_d * system_dict["N_GoodTransit"] / N_Trans_SNR1000 * 1000
    )

    planet_df = pd.DataFrame(data_dict)
    system_df = pd.DataFrame.from_records([system_dict])

    return planet_df, system_df


def RequestFromList_Transit(
    planet_list,
    catalog,
    observatory,
    d_start,
    d_end,
    observation_puffer=1 / 24,
    min_required_altitude=0,
    max_airm_good=2,
    max_sunalt_good=-20,
    SNR_method="Crires",
    mband="K",
    verbose=True,
):

    planet_transit_list = []
    system_list = []

    for planet in tqdm(planet_list):
        if verbose:
            print("Finding transits for " + planet)
        planet_df, system_df = TransitInformation(
            planet,
            catalog,
            observatory,
            d_start,
            d_end,
            observation_puffer,
            min_required_altitude,
            max_airm_good,
            max_sunalt_good,
            SNR_method,
            mband,
            verbose,
        )
        planet_transit_list.append(planet_df)
        system_list.append(system_df)

    alltrans = pd.concat(planet_transit_list, ignore_index=True)
    systemdat = pd.concat(system_list, ignore_index=True)

    systemdat_ranked = systemdat.sort_values(
        "systemscore", ascending=False, ignore_index=True
    )

    rounding_dict = {
        "SNR_nosmear": 0,
        "Mpl": 3,
        "Rpl": 3,
        "Period": 2,
        "t_nosmear": 1,
        "transit_depth": 3,
        "transit_dur": 3,
        "acc": 6,
        "K_p": 1,
        "SNR_5min": 0,
        "SNR_transdur": 0,
        "H": 0,
        "Delta_d": 6,
        "N_Trans_SNR1000": 3,
        "N_trans_req": 2,
        "systemscore": 2,
        "Sunalt_start": 2,
        "Sunalt_mid": 2,
        "Sunalt_end": 2,
        "Airm_start": 2,
        "Airm_mid": 2,
        "Airm_end": 2,
        "Moon_dist": 2,
        "V_bary": 2,
        "Sunalt_start_b": 2,
        "Sunalt_mid_b": 2,
        "Sunalt_end_b": 2,
        "Sunalt_start_a": 2,
        "Sunalt_mid_a": 2,
        "Sunalt_end_a": 2,
        "Airm_start_b": 2,
        "Airm_mid_b": 2,
        "Airm_end_b": 2,
        "Airm_start_a": 2,
        "Airm_mid_a": 2,
        "Moon_dist_b": 2,
        "Moon_dist_a": 2,
    }

    return alltrans.round(rounding_dict), systemdat_ranked.round(rounding_dict)


def EmissionInformation(
    planet,
    catalog,
    obs_time,
    time_in_eclipse,
    min_required_altitude,
    d_start,
    d_end,
    max_airm_good,
    max_sunalt_good,
    observatory,
    SNR_method="Crires",
    mband="K",
    verbose=True,
):

    """
    Everything with _b is before eclipse, with _a is after eclipse
    """
    if catalog == "TEP":
        physprop, hompar, hommes, obsplan = LoadTEP()
        planetData = GetPlanetDataTEP(planet)
    elif catalog == "Nexa":
        planetData = GetPlanetDataNexa(planet)
        if planetData == None:
            if verbose:
                print("Try to search in custom catalog.")
            planetData = GetPlanetDataCustom(planet)
    elif catalog == "Custom":
        planetData = GetPlanetDataCustom(planet)
    else:
        print("Requested catalog not recognized. Choices are TEP and Nexa")
        return None, None

    if planetData == None:
        if verbose:
            print(planet + " not contained in catalog!")
        return None, None

    if not CheckRequiredParameters(planetData):
        if verbose:
            print("At least one of the required parameters is missing!")
        return None, None

    jd_start = jdcnv(d_start)
    jd_end = jdcnv(d_end)

    ###
    ###First create dataframe with all observation opportunities
    ###
    t_min = jd_start
    t_max = jd_end

    T0 = planetData["T0"]
    period = planetData["orbPer"]
    Tdur = planetData["Tdur"]
    ra = planetData["ra"]
    dec = planetData["dec"]

    observatory_data = coordinates.EarthLocation.of_site(observatory)
    lon = observatory_data["longitude"]
    lat = observatory_data["latitude"]
    alt = observatory_data["altitude"]

    trnum_start = np.floor((t_min - T0) / period)
    trnum_end = np.ceil((t_max - T0) / period)
    # Relevant transit epochs
    tr = np.arange(trnum_start, trnum_end, 1)

    # Get list of all relevant times (start, mid and end for both before and after)
    t_list = []
    for epoch in tr:
        Tmid = T0 + float(epoch) * period + period / 2
        T_before = Tmid - Tdur / 2 + time_in_eclipse - obs_time / 2
        T_after = Tmid + Tdur / 2 - time_in_eclipse + obs_time / 2

        if (Tmid < t_min) or (Tmid > t_max):
            # This may happen because the transit may occur in the first
            # relevant epoch but still before tmin. Likewise for tmax.
            continue

        t_list.extend(
            [
                T_before - obs_time / 2,
                T_before,
                T_before + obs_time / 2,
                T_after - obs_time / 2,
                T_after,
                T_after + obs_time / 2,
            ]
        )

    t_list = np.array(t_list)

    altaz = eq2hor(
        t_list,
        (ra, dec),
        (lon, lat, alt),
    )
    altitude = altaz[0]

    # remove when planet not visible from observation site
    alt_filtered = []
    t_filtered = []

    nan_list = [np.nan, np.nan, np.nan]
    for i in range(int(altitude.size / 6)):
        minalt_before = np.where(altitude[i * 6 : i * 6 + 3] >= min_required_altitude)[
            0
        ]
        minalt_after = np.where(
            altitude[i * 6 + 3 : i * 6 + 6] >= min_required_altitude
        )[0]

        if (len(minalt_before) < 3) & (len(minalt_after) < 3):
            # target not visible at all -> skip
            continue
        if len(minalt_before) < 3:
            alt_filtered.extend(nan_list)
            t_filtered.extend(nan_list)
        else:
            alt_filtered.extend(altitude[i * 6 : i * 6 + 3])
            t_filtered.extend(t_list[i * 6 : i * 6 + 3])
        if len(minalt_after) < 3:
            alt_filtered.extend(nan_list)
            t_filtered.extend(nan_list)
        else:
            alt_filtered.extend(altitude[i * 6 + 3 : i * 6 + 6])
            t_filtered.extend(t_list[i * 6 + 3 : i * 6 + 6])
    alt_filtered = np.array(alt_filtered)
    t_filtered = np.array(t_filtered)

    total_nights = int(t_filtered.size / 6)

    if total_nights == 0:
        if verbose:
            print("No eclipse opportunity found for " + planet)

        return None, None

    # calculate sun altitude and airmass at each time
    notnan = pd.notnull(t_filtered)
    sunpos_radec = coordinates.get_sun(t_filtered[notnan])
    sunpos_altaz = eq2hor(
        t_filtered[notnan],
        (sunpos_radec[1][0],
        sunpos_radec[2][0]),
        (lon,
        lat,
        alt),
    )
    sunalt = np.ones(t_filtered.shape) * np.nan
    sunalt[notnan] = sunpos_altaz[0]

    airm = coordinates.AltAz(az=90 - alt_filtered)

    # calculate distance to moon (only for midpoints before and after)
    moondist = np.ones(int(t_filtered.size / 3)) * np.nan
    midpoints = np.array(
        [t_filtered[i * 3 + 1] for i in range(int(t_filtered.size / 3))]
    )

    mpos = coordinates.get_moon(midpoints[pd.notnull(midpoints)])
    moondist[pd.notnull(midpoints)] = mpos.seperation((ra, dec))

    ##Get the date of each night
    t_jd_b = np.array([t_filtered[i * 6 + 1] for i in range(total_nights)])
    t_jd_a = np.array([t_filtered[i * 6 + 4] for i in range(total_nights)])

    # rough determination of timezone
    timezone = ((int(lon / 15) + 12) % 24) - 12

    t_jd_b += timezone / 24
    t_jd_a += timezone / 24

    t_notnan_b = t_jd_b[pd.notnull(t_jd_b)]
    t_notnan_a = t_jd_a[pd.notnull(t_jd_a)]

    # Subtract 1 day if night starts at day before (if time is before 8am)
    t_notnan_b[t_notnan_b % 1 > 0.5] -= 1
    dates_complete_b = Time(t_notnan_b, format="jd").iso
    dates_b = [date[:10] for date in dates_complete_b]

    t_notnan_a[t_notnan_a % 1 > 0.5] -= 1
    dates_complete_a = Time(t_notnan_a, format="jd").iso
    dates_a = [date[:10] for date in dates_complete_a]

    night = [np.nan for i in range(total_nights)]
    # use night of b if available, because that is earlier than a
    for i in range(len(dates_a)):
        night_ind = np.where(pd.notnull(t_jd_a))[0][i]
        night[night_ind] = dates_a[i]
    for i in range(len(dates_b)):
        night_ind = np.where(pd.notnull(t_jd_b))[0][i]
        night[night_ind] = dates_b[i]

    data_dict = {
        "System": [planet for i in range(total_nights)],
        "Night": night,
        "Tmid_b": [t_filtered[i * 6 + 1] for i in range(total_nights)],
        "Tmid_a": [t_filtered[i * 6 + 4] for i in range(total_nights)],
        "Obs_start_b": [t_filtered[i * 6] for i in range(total_nights)],
        "Obs_end_b": [t_filtered[i * 6 + 2] for i in range(total_nights)],
        "Obs_start_a": [t_filtered[i * 6 + 3] for i in range(total_nights)],
        "Obs_end_a": [t_filtered[i * 6 + 5] for i in range(total_nights)],
        "Sunalt_start_b": [sunalt[i * 6] for i in range(total_nights)],
        "Sunalt_mid_b": [sunalt[i * 6 + 1] for i in range(total_nights)],
        "Sunalt_end_b": [sunalt[i * 6 + 2] for i in range(total_nights)],
        "Sunalt_start_a": [sunalt[i * 6 + 3] for i in range(total_nights)],
        "Sunalt_mid_a": [sunalt[i * 6 + 4] for i in range(total_nights)],
        "Sunalt_end_a": [sunalt[i * 6 + 5] for i in range(total_nights)],
        "Airm_start_b": [airm[i * 6] for i in range(total_nights)],
        "Airm_mid_b": [airm[i * 6 + 1] for i in range(total_nights)],
        "Airm_end_b": [airm[i * 6 + 2] for i in range(total_nights)],
        "Airm_start_a": [airm[i * 6 + 3] for i in range(total_nights)],
        "Airm_mid_a": [airm[i * 6 + 4] for i in range(total_nights)],
        "Airm_end_a": [airm[i * 6 + 5] for i in range(total_nights)],
        "Moon_dist_b": [moondist[i * 2] for i in range(int(moondist.size / 2))],
        "Moon_dist_a": [moondist[i * 2 + 1] for i in range(int(moondist.size / 2))],
    }

    data_dict["GoodCond_b"] = (
        (np.array(data_dict["Airm_start_b"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_start_b"]) < max_sunalt_good)
        & (np.array(data_dict["Airm_end_b"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_end_b"]) < max_sunalt_good)
    )

    data_dict["GoodCond_a"] = (
        (np.array(data_dict["Airm_start_a"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_start_a"]) < max_sunalt_good)
        & (np.array(data_dict["Airm_end_a"]) < max_airm_good)
        & (np.array(data_dict["Sunalt_end_a"]) < max_sunalt_good)
    )

    obs_df = pd.DataFrame(data_dict)

    ###
    ###Create df with system data that is the same for all transits of a planet
    ###
    system_dict = {}

    mag = mband + "mag"
    # system name
    system_dict["system"] = planet
    system_dict["Mstar"] = planetData["MsSun"]
    system_dict["Rstar"] = planetData["RsSun"]
    system_dict["Mpl"] = planetData["MpJ"]
    system_dict["Rpl"] = planetData["RpJ"]
    system_dict["Period"] = planetData["orbPer"]
    system_dict["Ecc"] = planetData["orbEcc"]
    system_dict["Teff"] = planetData["Teff"]
    system_dict[mag] = planetData[mag]

    system_dict["transit_dur"] = planetData["Tdur"]

    # stellar SNR
    # here units are km and s everywhere
    # first calculate maximum exposure so that lines do not shift between pixels
    pixel_size = 1  # pixel size in km/s
    K_p = (
        2
        * np.pi
        * (planetData["SMA"] * 1.496e8)
        / (planetData["orbPer"] * 24 * 60 * 60)
    )
    # Acceleration at time of transit is almost constant, the maximum of the acceleration sine curve
    acceleration = K_p * 2 * np.pi / (planetData["orbPer"] * 24 * 60 * 60)
    exposure_time = pixel_size / acceleration

    SNR = ComputeSNR(
        planetData[mag],
        exposure_time,
        np.nanmin(data_dict["Airm_mid_b"]),
        SNR_method,
        mband,
    )
    system_dict["SNR_nosmear"] = SNR
    system_dict["t_nosmear"] = exposure_time
    system_dict["acc"] = acceleration
    system_dict["K_p"] = K_p

    # SNR for fixed exposure time of 5 minutes
    SNR_5min = ComputeSNR(
        planetData[mag], 5 * 60, np.nanmin(data_dict["Airm_mid_b"]), SNR_method, mband
    )
    system_dict["SNR_5min"] = SNR_5min

    # scale height
    grav_constant = 6.67430 * 10 ** -11
    k_b = 1.380649 * 10 ** -23
    m_H = 1.6735575 * 10 ** -27  # mass of hydrogen atom [kg]
    surface_gravity = (
        grav_constant
        * planetData["MpJ"]
        * 1.898
        * 10 ** 27
        / (planetData["RpJ"] * 69911000) ** 2
    )
    mean_mol_weight = 2.3  # for hot jupiter atmosphere

    Teq = planetData["Teq"]
    if np.isnan(planetData["Teq"]):
        # if no value for Teq, calculate it using formula from Kempton 2018
        Teq = (
            planetData["Teff"]
            * np.sqrt(planetData["RsSun"] * 696340 / (planetData["SMA"] * 1.496e8))
            * 0.25 ** 0.25
        )
    H = k_b * Teq / (mean_mol_weight * m_H * surface_gravity)
    system_dict["Teq"] = Teq
    system_dict["H"] = H

    # number of transits with good conditions
    moondist_okay_mask_b = np.array(data_dict["Moon_dist_b"]) > 30
    goodtransit_mask_b = (data_dict["GoodCond_b"] == 1) & moondist_okay_mask_b
    moondist_okay_mask_a = np.array(data_dict["Moon_dist_a"]) > 30
    goodtransit_mask_a = (data_dict["GoodCond_a"] == 1) & moondist_okay_mask_a

    system_dict["N_GoodBefore"] = np.where(goodtransit_mask_b)[0].size
    system_dict["N_GoodAfter"] = np.where(goodtransit_mask_a)[0].size
    system_dict["N_GoodBoth"] = np.where(goodtransit_mask_b & goodtransit_mask_a)[
        0
    ].size

    # calculate score to rank the systems
    # reference wavelength to calculate ratio of planck spectra
    if mband == "K":
        lambda_ref = 2.2 * 10 ** -6
    elif mband == "J":
        lambda_ref = 1.25 * 10 ** -6
    else:
        lambda_ref = 2.2 * 10 ** -6
    planet_signal = (
        planck(Teq, lambda_ref) * (system_dict["Rpl"] * 69911000) ** 2
    )
    star_signal = (
        planck(system_dict["Teff"], lambda_ref)
        * (system_dict["Rstar"] * 696340000) ** 2
    )
    signal_strength = planet_signal / star_signal

    scoreNumberofNights = 1
    if (system_dict["N_GoodBefore"] == 0) & (system_dict["N_GoodAfter"] == 0):
        scoreNumberofNights = 0
    elif (system_dict["N_GoodBefore"] > 3) & (system_dict["N_GoodAfter"] > 3):
        scoreNumberofNights = 1.5

    score = SNR_5min * scoreNumberofNights * signal_strength * 100
    system_dict["systemscore"] = score

    system_df = pd.DataFrame.from_records([system_dict])

    return obs_df, system_df


def RequestFromList_Emission(
    planet_list,
    catalog,
    observatory,
    d_start,
    d_end,
    obs_time,
    time_in_eclipse,
    min_required_altitude=0,
    max_airm_good=2,
    max_sunalt_good=-20,
    SNR_method="Crires",
    mband="K",
    verbose=True,
):

    planet_event_list = []
    system_list = []

    for planet in tqdm(planet_list):
        if verbose:
            print("Finding emission observation opportunities for " + planet)
        planet_df, system_df = EmissionInformation(
            planet,
            catalog,
            obs_time,
            time_in_eclipse,
            min_required_altitude,
            d_start,
            d_end,
            max_airm_good,
            max_sunalt_good,
            observatory,
            SNR_method,
            mband,
            verbose,
        )
        planet_event_list.append(planet_df)
        system_list.append(system_df)

    allevents = pd.concat(planet_event_list, ignore_index=True)
    systemdat = pd.concat(system_list, ignore_index=True)

    systemdat_ranked = systemdat.sort_values(
        "systemscore", ascending=False, ignore_index=True
    )

    rounding_dict = {
        "SNR_nosmear": 0,
        "Mpl": 3,
        "Rpl": 3,
        "Period": 2,
        "t_nosmear": 1,
        "transit_depth": 3,
        "transit_dur": 3,
        "acc": 6,
        "K_p": 1,
        "SNR_5min": 0,
        "SNR_transdur": 0,
        "H": 0,
        "Delta_d": 6,
        "N_Trans_SNR1000": 3,
        "N_trans_req": 2,
        "systemscore": 2,
        "Sunalt_start": 2,
        "Sunalt_mid": 2,
        "Sunalt_end": 2,
        "Airm_start": 2,
        "Airm_mid": 2,
        "Airm_end": 2,
        "Moon_dist": 2,
        "V_bary": 2,
        "Sunalt_start_b": 2,
        "Sunalt_mid_b": 2,
        "Sunalt_end_b": 2,
        "Sunalt_start_a": 2,
        "Sunalt_mid_a": 2,
        "Sunalt_end_a": 2,
        "Airm_start_b": 2,
        "Airm_mid_b": 2,
        "Airm_end_b": 2,
        "Airm_start_a": 2,
        "Airm_mid_a": 2,
        "Moon_dist_b": 2,
        "Moon_dist_a": 2,
    }

    return allevents.round(rounding_dict), systemdat_ranked.round(rounding_dict)
