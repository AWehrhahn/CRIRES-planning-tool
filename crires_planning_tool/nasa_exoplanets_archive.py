from astroquery.utils.tap.core import TapPlus
from astroquery.nasa_exoplanet_archive import (
    NasaExoplanetArchive as _NasaExoplanetArchive,
)


class NasaExoplanetsArchive:
    def __init__(self):
        #:int: the number of tries before timing out the query
        self.timeout = 10
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
        self._archive = TapPlus(url=value)
        self._url = value

    @property
    def archive(self):
        return self._archive

    @archive.setter
    def archive(self, value):
        raise AttributeError("Set the archive via the url property")

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
        data = None
        for i in range(self.timeout):
            try:
                # data = NasaExoplanetArchive.query_object(name)
                query = self.archive.launch_job(asql_query)
                data = query.results
                break
            except ConnectionError:
                print(f"Connection failed, attempt {i} of {self.timeout}")
                continue

        if data is None and i == self.timeout - 1:
            raise ConnectionError("Connection to the Nasa Exoplanet Archive timed out")
        elif data is None:
            raise RuntimeError(f"Invalid query: {asql_query}")

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
        if regularize:
            name = _NasaExoplanetArchive._regularize_object_name(name)
        asql_query = f"SELECT top 100 * FROM {self.table} WHERE pl_name='{name}'"
        data = self.tap(asql_query)
        return data

    def query_objects(self, names, regularize=True):
        aliases = {}
        if regularize:
            for i, name in enumerate(names):
                alias = _NasaExoplanetArchive._regularize_object_name(name)
                aliases[alias] = name
                names[i] = alias
        else:
            aliases = {n:n for n in names}

        asql_query = f"SELECT top 100 * FROM {self.table} WHERE "
        asql_query += " OR ".join([f"pl_name='{name}'" for name in names])

        data = self.tap(asql_query)
        return data, aliases
