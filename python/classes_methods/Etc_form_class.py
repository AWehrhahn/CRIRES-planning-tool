#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:30:00 2020

This file contains the class etc_form to read in, change and update the input file for the ETC calculator 'etc-form.json'.
IMPORTANT: Do not change the files 'etc-form-default-snr.json' or 'etc-form-default-ndir.json' 
except if necessary due to updates on the etc-cli side. .


@author: jonaszbinden
GitHub: jonaszubindu
"""
import argparse
import json
import logging
import os

from astropy import units as u

# import pandas as pd
import time
from functools import wraps
from json import JSONDecodeError, JSONEncoder
from os.path import dirname, join

import numpy as np
import requests

from classes_methods import etc_cli, misc

# from requests import Session
from requests_futures.sessions import FuturesSession

# TODO! workaround until we figure put how to handle ssl certificate correctly
import warnings
warnings.filterwarnings('ignore', message='Unverified HTTPS request')

##########################################################################################################

""" Warnings """
JSONDecodeWarning = Warning(
    "Something went wrong processing the etc-form file... I will try to find the error for you"
)
NDITWarning = Warning("NDIT not available from output file")

logger = logging.getLogger(__name__)

def DecodeWarning(key, value):
    DecodeWarning = FutureWarning(
        f"the error is related to the present {key} input value: {value.__str__()}"
    )
    return DecodeWarning


ErrorNotFoundWarning = DeprecationWarning(
    "Sorry, I cannot find the error, check the file etc-format.json and try to run it manually on the ETC calculator webpage. \n Maybe she can help you find the error..."
)

dit_std = 10  # default value for DIT
nditSTD = 1  # default value for NDIT

##########################################################################################################


class FormEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


##########################################################################################################


class etc_form:
    """
    Include ETC constraints here as a different mode to compute
    additional independent constraints

    This can be advanced by any method to change some input parameter of
    'etc-form.json' for any type of targets.

    WARNING: If the general structure of the file changes due to for instance
             change from inputtype "Spectrum" to "Emission Line", this must be regarded
             when adding methods to alter 'etc-form.json'. Might conflict with other methods!
    """

    def __init__(self, inputtype):
        """
            Initializes 'etc-form-default.json' via json to a namespace object according to inputtype ''ndit'' or ''snr''.

            Parameters:
            -----------
            inputtype : string
                specify if the ETC-calculator should be run in S/N mode ''snr-Templ'' or in
                NDIT mode ''ndit-Templ'' if spectral templetas from MARCS catalogue are used, if the spectrum is assumed to
                be blackbody with an effective Temperature Teff, then use ''snr-Teff'' and ''ndit-Teff'' respectively.

        """
        filenames = {
            "snr-Teff": "etc-form-default-snr-Teff.json",
            "ndit-Teff": "etc-form-default-ndit-Teff.json",
            "snr-Templ": "etc-form-default-ndit-Templ.json",
            "ndit-Templ": "etc-form-default-ndit-Templ.json",
        }
        filename = join(dirname(__file__), "../json_files", filenames[inputtype])

        with open(filename) as args:
            etc_obj = json.load(args)
        self.input = etc_obj

    ##########################################################################################################

    def update_etc_form(self, **kwargs):
        """
            changes input values in 'etc-form.json'

            Parameters:
            -----------
            Keyword arguments recognized by update_etc_form:

            airmass : float

            moon_target_sep : list
                Two values, first value is moon_target_separation in degrees, second value is moon_alt in degrees.

            moon_phase : float
                moon_sun_separation in degrees.

            snr : int or float
                Minimum signal to noise ratio S/N.

            dit : int or float
                DIT exposure time for single exposure.

            ndit : int
                NDIT number of single exposures for one single observation.

                NDIT*DIT = Texp total exposure time for one single observation.

            inputtype : string
                snr or ndit depending on ETC calculator should calculate the NDIT for a certain minimum S/N
                or S/N for a certain NDIT.

            temperature : float
                Effective temperature of the target object.

            brightness : float
                Object brightness, standard is J-band magnitude, system: AB.
            
            gsmag : float 
                Brightness of the guide star. 
            
            sptype : string
                Spectral type of guide star.

            others can be added:...

        """

        if "airmass" in kwargs:
            self.input["sky"]["airmass"] = kwargs["airmass"].to_value(1)
            # self.input.sky.airmass = 12 # Chabis Test

        if "moon_target_sep" in kwargs:
            self.input["sky"]["moon_target_sep"] = kwargs["moon_target_sep"][0].to_value(u.deg)
            self.input["sky"]["moon_alt"] = kwargs["moon_target_sep"][1].to_value(u.deg)

        if "moon_phase" in kwargs:
            self.input["sky"]["moon_sun_sep"] = kwargs["moon_phase"].to_value(u.deg)

        self.input["timesnr"]["snr"]["snr"] = kwargs.get("snr", 100)
        if "dit" in kwargs:
            self.input["timesnr"]["dit"] = kwargs.get("dit")

        if self.input["target"]["sed"]["spectrum"]["spectrumtype"] == "template":
            if "temperature" in kwargs:
                temperature = kwargs["temperature"].to_value(u.K)
                if temperature > 8000:
                    logger.warning(
                        f"WARNING : Temperature exceeds MARCS spT catalog levels! Teff = {temperature}, taking T = 8000 K"
                    )
                    temperature = 8000
                elif temperature < 4000:
                    logger.warning(
                        f"WARNING : Temperature does not reach lower MARCS spT catalog levels! Teff = {temperature}, taking T = 4000 K"
                    )
                    temperature = 4000
                else:
                    temperature = int(np.round(temperature / 500) * 500) 
                self.input["target"]["sed"]["spectrum"]["params"]["spectype"] = f"p{temperature}:g+4.0:m0.0:t02:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00"

        elif self.input["target"]["sed"]["spectrum"]["spectrumtype"] == "blackbody":
            if "temperature" in kwargs:
                self.input["target"]["sed"]["spectrum"]["params"]["temperature"] = kwargs["temperature"].to_value(u.K)

        if "gsmag" in kwargs:
            self.input["seeingiqao"]["params"]["gsmag"] = kwargs["gsmag"].to_value(u.mag)

        if "sptype" in kwargs:
            self.input["seeingiqao"]["params"]["sptype"] = kwargs["sptype"]

        if "brightness" in kwargs:
            self.input["target"]["brightness"]["params"]["mag"] = kwargs.get("brightness").to_value(u.mag)

        if "inputtype" in kwargs:
            self.input["timesnr"]["inputtype"] = kwargs.get("inputtype")

        if self.input["timesnr"]["inputtype"] == "snr":
            self.input["timesnr"]["dit"] = kwargs.get("dit", dit_std)

        elif self.input["timesnr"]["inputtype"] == "ndit":
            self.input["timesnr"]["dit"] = kwargs.get("dit", dit_std)
            self.input["timesnr"]["ndit"] = kwargs.get("ndit", dit_std)
        return self.input

    ##########################################################################################################


    def write_etc_format_file(self, input_filename=None):
        """
            Writes self.etc to a new JSON file named 'etc-form.json' such
            that it can be interpreted by the ETC online-calculator.
        """
        if input_filename is None:
            path = join(dirname(__file__), "../json_files")
            input_filename = join(path, "etc-form.json")

        with open(input_filename, "w") as dumpfile:
            json.dump(self.input, dumpfile, indent=2, cls=FormEncoder)

    ##########################################################################################################

    def run_etc_calculator(self, input_data, output_filename=None):
        """
            Runs ETC calculator through commandline and asks for output data file

            Parameters
            ----------
            name : str
                Name of the target for which the ETC shoudl calculate S/N.

            tim : datetime.datetime or astropy.time.Time
                Observation time for which the ETC should calculate S/N.

            Returns
            -------
            NDIT : int
                Number of single exposures with DIT to reach
                signal to noise S/N as defined in 'etc-form.json'.

            output : pandas DataFrame
                DataFrame object containing the output data from the
                ETC calculator

        """
        # Try to access the ETC server, give it a few tries if necessary
        future = call_ETC("crires", input_data, output_filename)
        output = future.result().data
        return output

    def get_ndit(self, output):
        if self.input["timesnr"]["inputtype"] == "ndit":
            ndit = self.input["timesnr"]["ndit"]
        else:
            ndit = output["data"]["time"]["ndit"]
        return ndit

    ##########################################################################################################

    # def etc_debugger(
    #     self,
    #     inputtype,
    #     name,
    #     tim,
    #     temperature,
    #     brightness,
    #     airmass,
    #     moon_phase,
    #     moon_target_sep,
    #     gsmag,
    #     snr=None,
    # ):
    #     """
    #         This tries to find the error in the etc-format file. As soon as the ETC calculator gets updated with better input error handling
    #         this function must be updated or replaced by additional error handling in the functions running the ETC calculator.

    #         Parameters
    #         ----------
    #         JSONDecodeError : Exception
    #             Handle of the JSONDecodeError that occurred while running the ETC calculator.

    #         temperature : float
    #             Temperature input parameter that was used.

    #         brightness : float
    #             Brightness input parameter that was used.

    #         airmass : float
    #             Airmass input parameter that was used.

    #         moon_phase : float
    #             Illumination of the moon, also known as moon_sun_separation.

    #         moon_target_sep : list
    #             Two values, first value is moon_target_separation, second value is moon_alt, altitude above horizon of the moon.

    #         gsmag : float 
    #             Brightness of the guide star. 
                
    #         Raises
    #         ------
    #         JSONDecodeError
    #             If the errornous parameter was found, raises the JSONDecodeError and reviels the faulty parameter.

    #         Returns
    #         -------
    #         None. If no raises occur, the etc_debugger tells the user that it has not found the error and gives the problem
    #         back to the user

    #     """

    #     path = join(dirname(__file__), "../json_files")

    #     print(
    #         "Something went wrong processing the etc-form file... I will try to fix it for you"
    #     )
    #     print(name, tim, temperature, brightness, airmass, moon_phase, moon_target_sep)
    #     os.system(f"cp {path}/etc-form.json {path}/etc-form-copy.json")

    #     cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(temperature=temperature)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("temperature", temperature)  # Warning
    #     cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(brightness=brightness)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("brightness", brightness)  # Warning
    #     cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(airmass=airmass)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("airmass", airmass)
    #     cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(moon_phase=moon_phase)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("moon_phase", moon_phase)
    #         cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(moon_target_sep=moon_target_sep)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("moon_target_sep", moon_target_sep)
    #         cls = type(self)
    #     ETC = cls.__new__(cls)
    #     ETC.__init__(inputtype)
    #     ETC.update_etc_form(gsmag=gsmag)
    #     ETC.write_etc_format_file()
    #     try:
    #         NDIT, output = ETC.run_etc_calculator(name, tim)
    #     except JSONDecodeError:
    #         raise DecodeWarning("gsmag", gsmag)

    #     # others...

    #     print("I will continue with the next planet for now...")
    #     raise ErrorNotFoundWarning  # Warning


##########################################################################################################

def response_hook(resp, *args, **kwargs):
    # parse the json storing the result on the response object
    try:
        resp.data = resp.json()
    except JSONDecodeError:
        pass

def send_request(post_data, url, uploadfile=None):
    session = FuturesSession()
    session.hooks['response'] = response_hook

    kwargs = {"data": json.dumps(post_data), "verify": False}

    if uploadfile is None:
        kwargs["headers"] = {'Content-Type': 'application/json'}
    else:
        kwargs["files"] = {"target": open(uploadfile, 'rb')}
    future = session.post(url, **kwargs)
    return future

def call_ETC(instrument, input_data, output_filename=None):
    """ Call the exposure time calculator API

    Parameters
    ----------
    instrument : str
        Name of the instrument, here "crires"
    input_filename : str
        Name of the json document with the input settings
    output_filename : str, optional
        Name of the json document to save the results in, by default None

    Returns
    -------
    data: dict
        ETC retrieved data
    """    

    """ Temporary this is in use, does not require VPN """
    baseurl = "https://etctestpub.eso.org/observing/etc/etcapi/"
    # baseurl = 'https://etc.eso.org/observing/etc/etcapi/'
    # baseurl = 'http://localhost:8000/observing/etc/etcapi/'
    # baseurl = 'https://etctest.hq.eso.org/observing/etc/etcapi/'

    etcName = etc_cli.getEtcUrl(instrument)
    url = baseurl + etc_cli.getEtcUrl(etcName)
    future = send_request(input_data, url)
    return future
    # ------------------------------------------------------------------------------------------------


##########################################################################################################

