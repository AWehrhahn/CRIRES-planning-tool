

# CRIRES-planning-tool documentation

Version 0.1.0 2021

#### Introduction:

The CRIRES-planning-tool is intended to be used to plan transit observations of exoplanets for CRIRES+, the new cross-dispersed high-resolution infrared spectrograph for the ESO VLT [CRIRES+](https://www.eso.org/sci/facilities/develop/instruments/crires_up.html). Observaton of exoplanets can be planned in two ways. Single candidate by name in a given timespan or by using constraints for observable candidates by CRIRES+, which can be loaded from the file: Nasa_Archive_Selection.txt (see section: **Constraints for Candidates**). The known exoplanets fulfilling these constraints are downloaded from NASA Exoplanet Archive and each candidate is checked for its observability from Cerro Paranal, Chile for during a given time frame. Each observable candidate is checked for a minimum signal-to-noise ratio (S/N)≥100 during 20 exposures. This constrain will likely be updated in the future to distinguish between observations of super-Earth's or giant planets, and sensitive to the host star spectral type. Each exposure is related to its total exposure time, calculated from the detector integration time multiplied with the number of detector integrations: (TEXP = DIT x NDIT) and NDIT is optimized to be within 16≤NDIT≤32 for each exposure (see section: **Exposure Time Calculator**). Candidates with observability of the complete transit are added to the observation list and further information can be found in the output excel files (see section: **Result files**). The tool uses two ways to calculate the number of exposures possible during a single transit. The details are described in my master thesis: *Planning observations of terrestrial Exoplanets around M type Stars with CRIRES+*, section *3.3 Signal-to-noise ratio and the Exposure Time Calculator*. The tool comes with plotting tools and a commandline window to access its functionalities. This document shall give an overview about the functionalities, accessibility, structure,  installation, and further development possibilities of the CRIRES-planning-tool. Code documentation can be found in **Code documentation** and a dependency tree is presented in **Dependencies**. The methods used for astronomical calculations are used from the astropy and astroplan library, about which documentation can be found here: [astroplan](https://astroplan.readthedocs.io/en/latest/), [astropy](https://docs.astropy.org/en/stable/).  



#### Installation:

1. Navigate to your chosen directory to install the CRIRES-planning-tool.

2. Download github repository: 
   `git clone https://github.com/AWehrhahn/CRIRES-planning-tool`

3. (optional) Setup a virtual environment to install the correct packages to run the planning tool

   `cd CRIRES-planning-tool/python`

   `pip install virtualenv`

   `virtualenv --python python3.7 [name of your venv]`

   activate your virtual environment:

   `source [name of your venv]/bin/activate`

   after you are done with running CRIRES-planning-tool use 

   `deactivate`

   to deactivate the virtual environment. 

4. install the requirements to run CRIRES-planning-tool stored in requirements.txt

   `pip  install -r requirements.txt`

5. install the package itself

    `pip install -e .`

   
To update the tool in case of changes you can use `git pull` and run the tool in your virtual environment, by activating it and installing the requirements with pip. Just as if it was freshly installed.

#### Commandline Menu
The command line interface of the planning tool can be accessed via `python -m crires-planning-tool`.
It takes at least three parameters, first and second are the starting and end date. This is followed by a list of planet or star names to run the tool for. Use the `-h` option to get detailed descriptions of all options. An example usecase is the following:
`python -m crires-planning-tool 01-01-2010 01-01-2020 "TRAPPIST-1"`
Note that you need to use the `--output=<filename>` option to generate output files, and you can use the `--plot` option to create visualization plots of the results.

#### Constraints for Candidates

By default the possible observations are constrained by four measures as described in `transit_planner.get_default_constraints()`. These can be replaced by passing the constraints parameter to `transit_planner.transit_calculation()`. The default constraints are: Altitude above 30 deg, airmass less than 1.7, 
astronomical twilight at the observing location, and at least 45 deg seperation from the moon.


#### Exposure Time Calculator

For regular operation we use an approximation of the ETC to determine the SNR much faster. For more accurate values one should call the ETC manually using the `EtcForm` class in `etc_form.py` as described below.

The exposure time calculator is called through a client and requires a json input file containing all the input data to compute the exposure time or the signal-to-noise ratio for a  particular observation. The exposure time calculator is provided by ESO and maintained by [Jakob Vinther](j.vinther@eso.org). The theory behind the ETC can be looked up in my thesis: *Planning observations of terrestrial Exoplanetsaround M type Stars with CRIRES+*, section *3.3 Signal-to-noise ratio and the Exposure Time Calculator*. The public interface can be accessed [here](https://etctestpub.eso.org/observing/etc/crires). Any updates of the etc conflicting with the CRIRES-planning tool should be checked in correspondence with Jakob Vinther. Here are a few reasons why CRIRES-planning-tool might not be able to access the ETC anymore and strategies to solve it:

1. The baseurl to call the ETC with the cli has changed. You can change the baseurl in the file: etc_form.py

2. The structure of the input json file has changed. There are several ways to fix this. The easiest way is by accessing the api version of the etc and plugging in standard inputs. Then replace the correct file in json_files/ 


one can download the jsonfile and depending on the desired input method as one of the following. Before you store you store the file, make sure that you make a copy of the old json file(s). 

​	calculating S/N, using spectral templates -> store file as: 

> etc-form-default-snr-Templ.json

​	calculating S/N, using effective temperature of the target -> store file as: 

> etc-form-default-snr-Teff.json

​	or calculating exposure time for minimum S/N, using spectral templates: -> store file as: 

> etc-form-default-ndit-Templ.json

​	or calculating exposure time for minimum S/N, using effective temperature of the target: -> store file as: 

>etc-form-default-ndit-Teff.json

Check for differences between the old and the new file. Check in the script Etc_form_class.py if the function update_etc_form is still following the right structure to write input data into the replaced json file and adjust the structure adequately. To test the structure of the input json file, navigate to CRIRES-planning-tool/python, open an iPython console and type the following:

```python
from classes_methods import Etc_form_class
etc_form = Etc_form_class.etc_form('[jsonfile-type]')
```

and write the desired type of json input file at [jsonfile-type]: snr-Templ, edit-Templ, snr-Teff, ndit-Teff. Now you can investigate the structure of etc_form by writing [etc_form] + [.] + [Tab] and navigate through the file... 

If none of these two strategies solve the problem, you need to contact Jakob Vinther.


#### Result files

The results are stored in csv files, each column has a header that describes its contents.
