#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Handles the data access to the Nasa Exoplanet Archive
Especially since it does not use TAP in  astroquery (yet)

@author: Ansgar Wehrhahn
"""
import json
import logging
import urllib
from tempfile import NamedTemporaryFile

from astropy.table import Table
from astropy.utils.data import (clear_download_cache, download_file,
                                import_file_to_cache)
from astroquery.nasa_exoplanet_archive import \
    NasaExoplanetArchive as _NasaExoplanetArchive
from astroquery.utils.tap.core import TapPlus

logger = logging.getLogger(__name__)


class NasaExoplanetsArchive:
    def __init__(self, timeout=10, verbose=0):
        #:int: the number of tries before timing out the query
        self.timeout = timeout
        #:bool: verbose TAP connection
        self.verbose = verbose
        #:str: the url of the TAP interface
        self.url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
        #:str: the table to query, should be "pscomppars" or "ps"
        self.table = "pscomppars"
        #:list: a list of citations to use
        self.citation = [
            "This research has made use of the NASA Exoplanet "
            "Archive, which is operated by the California Institute of "
            "Technology, under contract with the National Aeronautics and "
            "Space Administration under the Exoplanet Exploration Program."
        ]

    @property
    def url(self):
        return self._url

    @url.setter
    def url(self, value):
        self._archive = TapPlus(url=value, verbose=self.verbose >= 2)
        self._url = value

    @property
    def archive(self):
        return self._archive

    @archive.setter
    def archive(self, value):
        raise AttributeError("Set the archive via the url property")

    def reset_cache(self):
        clear_download_cache(pkgname="crires-planning-tool")

    def get_filename_from_cache(self, filename):
        try:
            url = "file:" + urllib.parse.quote(filename)
            fname = download_file(
                url, cache=True, pkgname="crires-planning-tool", show_progress=False
            )
        except:
            fname = None
        return fname

    def import_file_to_cache(self, key, filename):
        url = "file:" + urllib.parse.quote(key)
        import_file_to_cache(
            url, filename, pkgname="crires-planning-tool", replace=True
        )

    def from_cache_table(self, asql_query):
        try:
            fname = self.get_filename_from_cache(asql_query)
            data = Table.read(fname, format="pandas.csv")
        except (KeyError, urllib.error.URLError):
            data = None
        except:
            data = None
        return data

    def to_cache_table(self, asql_query, data):
        with NamedTemporaryFile("w", suffix=".csv") as tempfile:
            data.write(tempfile, format="pandas.csv")
            tempfile.flush()
            self.import_file_to_cache(asql_query, tempfile.name)

    def from_cache_json(self, fname):
        fname = self.get_filename_from_cache("aliases.json")
        try:
            with open(fname) as f:
                data = json.load(f)
        except:
            data = None
        return data

    def to_cache_json(self, fname, data):
        with NamedTemporaryFile("w", suffix=".json") as tempfile:
            json.dump(data, tempfile)
            tempfile.flush()
            self.import_file_to_cache(fname, tempfile.name)

    def tap(self, asql_query):
        """Get data from the TAP interface

        Parameters
        ----------
        asql_query : str
            The query to send to the TAP interface in Astronomy SQL format.

        Returns
        -------
        data : astropy.Table
            the returnded data

        Raises
        ------
        ConnectionError
            Could not establish a connection
        RuntimeError
            Could not get results with the given query
        """

        data = self.from_cache_table(asql_query)
        if data is not None:
            if self.verbose >= 1:
                logger.info("Retrieving catalog data from local cache")
            return data

        for i in range(self.timeout):
            try:
                # data = NasaExoplanetArchive.query_object(name)
                query = self.archive.launch_job(asql_query, verbose=self.verbose >= 1)
                data = query.results
                break
            except ConnectionError:
                if self.verbose >= 1:
                    logger.warning(f"Connection failed, attempt {i} of {self.timeout}")
                continue

        if data is None and i == self.timeout - 1:
            raise ConnectionError("Connection to the Nasa Exoplanet Archive timed out")
        elif data is None:
            raise RuntimeError(f"Invalid query: {asql_query}")

        # Store it in the cache for later
        self.to_cache_table(asql_query, data)

        return data

    def query_object(self, name, regularize=True):
        """Get data for a specific target

        Parameters
        ----------
        name : str
            Name of the target in the format "star planet"
        regularize : bool, optional
            Whether to check the name against the known aliases in the database.
            The database otherwise only accepts exact matches. By default True

        Returns
        -------
        data : astropy.Table
            Table with the results for the target
        """
        return self.query_object([name], regularize=regularize)

    def query_objects(self, names, regularize=True):
        """
        Request information for a set of planets or stars

        Parameters
        ----------
        names : list
            list of planet names
        regularize : bool, optional
            whether to check the database for the correct names, by default True.
            This creates overhead but is more convenient.

        Returns
        -------
        [type]
            [description]
        """
        original = {}
        if regularize:
            aliases = self.from_cache_json("aliases.json")
            aliases = {} if aliases is None else aliases
            hasChanged = False
            for i, name in enumerate(names):
                if name in aliases.keys():
                    alias = aliases[name]
                else:
                    alias = _NasaExoplanetArchive._regularize_object_name(name)
                    hasChanged = True
                    aliases[name] = alias
                original[alias] = name
                names[i] = alias
            if hasChanged:
                self.to_cache_json("aliases.json", aliases)

        asql_query = f"SELECT * FROM {self.table} WHERE "
        asql_query += " OR ".join([f"pl_name like '%{name}%'" for name in names])

        data = self.tap(asql_query)

        if regularize:
            for i, d in enumerate(data):
                try:
                    name = d["pl_name"]
                    name = original[name]
                except KeyError:
                    name = d["hostname"]
                    name = original.get(name, name)
                    name = f"{name} {d['pl_letter']}"
                data[i]["pl_name"] = name

        return data

    def query_criteria(self, criteria):
        asql_query = f"SELECT * FROM {self.table} WHERE "
        asql_query += " AND ".join(criteria)
        data = self.tap(asql_query)
        return data
