from astroquery.utils.tap.core import TapPlus
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

class NasaExoplanetsArchive:
    def __init__(self):
        #:int: the number of tries before timing out the query
        self.timeout = 10
        #:str: the url of the TAP interface
        self.url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
        #:str: the table to query, should be "pscomppars" or "ps"
        self.table = "pscomppars"
        #:list: a list of citations to use
        self.citation = ["This research has made use of the NASA Exoplanet "
            "Archive, which is operated by the California Institute of "
            "Technology, under contract with the National Aeronautics and "
            "Space Administration under the Exoplanet Exploration Program."]

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
        """ Get data using the TAP interface """
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
        """ Get data for a specific object """
        if regularize:
            name = NasaExoplanetArchive._regularize_object_name(name)
        asql_query = f"SELECT top 100 * FROM {self.table} WHERE hostname='{name}'"
        data = self.tap(asql_query)
        return data
