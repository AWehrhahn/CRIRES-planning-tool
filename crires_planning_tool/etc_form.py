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
import json
import logging
import warnings
from json import JSONDecodeError
from os.path import dirname, join

import numpy as np
from astropy import units as u
from requests_futures.sessions import FuturesSession

# TODO! workaround until we figure put how to handle ssl certificate correctly
warnings.filterwarnings("ignore", message="Unverified HTTPS request")

logger = logging.getLogger(__name__)

dit_std = 10  # default value for DIT
nditSTD = 1  # default value for NDIT


class EtcForm:
    """
    Include ETC constraints here as a different mode to compute
    additional independent constraints

    This can be advanced by any method to change some input parameter of
    'etc-form.json' for any type of targets.

    WARNING: If the general structure of the file changes due to for instance
             change from inputtype "Spectrum" to "Emission Line", this must be regarded
             when adding methods to alter 'etc-form.json'. Might conflict with other methods!
    """

    def __init__(self, filename, instrument="crires"):
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
        if filename in filenames.keys():
            filename = join(dirname(__file__), "json_files", filenames[filename])
        with open(filename) as args:
            etc_obj = json.load(args)

        self.input = etc_obj
        self.output = None
        self.timeout = 5
        self.instrument = instrument
        self.baseurl = "https://etctestpub.eso.org/observing/etc/etcapi/"

    def __getitem__(self, key):
        return self.input[key]

    def __setitem__(self, key, value):
        self.input[key] = value

    def update_form(self, **kwargs):
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
            self.input["sky"]["moon_target_sep"] = kwargs["moon_target_sep"][
                0
            ].to_value(u.deg)
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
                self.input["target"]["sed"]["spectrum"]["params"][
                    "spectype"
                ] = f"p{temperature}:g+4.0:m0.0:t02:st:z+0.00:a+0.00:c+0.00:n+0.00:o+0.00"

        elif self.input["target"]["sed"]["spectrum"]["spectrumtype"] == "blackbody":
            if "temperature" in kwargs:
                self.input["target"]["sed"]["spectrum"]["params"][
                    "temperature"
                ] = kwargs["temperature"].to_value(u.K)

        if "gsmag" in kwargs:
            self.input["seeingiqao"]["params"]["gsmag"] = kwargs["gsmag"].to_value(
                u.mag
            )

        if "sptype" in kwargs:
            self.input["seeingiqao"]["params"]["sptype"] = kwargs["sptype"]

        if "brightness" in kwargs:
            self.input["target"]["brightness"]["params"]["mag"] = kwargs.get(
                "brightness"
            ).to_value(u.mag)

        if "inputtype" in kwargs:
            self.input["timesnr"]["inputtype"] = kwargs.get("inputtype")

        if self.input["timesnr"]["inputtype"] == "snr":
            self.input["timesnr"]["dit"] = kwargs.get("dit", dit_std)

        elif self.input["timesnr"]["inputtype"] == "ndit":
            self.input["timesnr"]["dit"] = kwargs.get("dit", dit_std)
            self.input["timesnr"]["ndit"] = kwargs.get("ndit", dit_std)
        return self.input

    def write(self, filename):
        """
            Writes self.etc to a new JSON file named 'etc-form.json' such
            that it can be interpreted by the ETC online-calculator.
        """
        with open(filename, "w") as dumpfile:
            json.dump(self.input, dumpfile, indent=2)

    @classmethod
    def read(cls, filename):
        form = EtcForm(filename)
        return form

    def response_hook(self, resp, *args, **kwargs):
        # parse the json storing the result on the response object
        try:
            resp.data = resp.json()
        except JSONDecodeError:
            pass

    def send_request(self, post_data, uploadfile=None):
        """ Send a request to the ETC with a set of input data

        Parameters
        ----------
        post_data : dict
            input data for the etc
        uploadfile : [type], optional
            [description], by default None

        Returns
        -------
        [type]
            [description]
        """
        session = FuturesSession()
        session.hooks["response"] = self.response_hook

        kwargs = {"data": json.dumps(post_data), "verify": False}

        if uploadfile is None:
            kwargs["headers"] = {"Content-Type": "application/json"}
        else:
            kwargs["files"] = {"target": open(uploadfile, "rb")}
        
        url = self.get_etc_url()
        future = session.post(url, **kwargs)
        return future

    def get_etc_url(self):
        """ Create the url for the ETC and this instrument

        Returns
        -------
        url : str
            Url of the ETC
        """
        if('4most' in self.instrument.lower()):
            return self.baseurl + '4Most/'
        elif('crires' in self.instrument.lower()):
            return self.baseurl + 'Crires2/'
        # elif('harps' in etcname.lower()):  # not yet working
        #     return 'HarpsNirps/'
        else:
            print("error: no match for etcname: " + self.instrument)

    def call(self):
        """ 
        Call the exposure time calculator API
        The input settings are defined in self.input

        Returns
        -------
        data: dict
            ETC retrieved data
        """
      
        response = None
        n_attempts = max(self.timeout, 1)
        for i in range(n_attempts):
            future = self.send_request(self.input)
            try:
                response = future.result()
                response = response.data
            except ConnectionError as ex:
                logger.warning(str(ex))
            

            # check for success
            if response is not None and response["success"]:
                break
            else:
                # unfortunately the etc does not tell us why it failed
                # so we just try again
                response = None
                logger.warning("ETC failed, trying again...")

        if response is None:
            raise ConnectionError("Could not connect to the ETC")

        self.output = response
        return response

