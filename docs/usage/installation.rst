Installation
============


1. Navigate to your chosen directory to install the CRIRES-planning-tool.

2. Download github repository: 
   `git clone https://github.com/AWehrhahn/CRIRES-planning-tool`

3. Enter the new directory
   `cd CRIRES-planning-tool`

3. (optional) Setup a virtual environment for the planning tool
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

To update the tool in case of changes you can use `git pull` and run the tool
in your virtual environment, by activating it and installing the requirements with pip, 
just as if it was freshly installed.
