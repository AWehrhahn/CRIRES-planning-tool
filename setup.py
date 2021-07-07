from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="crires-planning-tool",
    description="Plans observations for exoplanet transits for the CRIRES+ instrument",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Ansgar Wehrhahn, Linn Boldt-Christmas, Jonas Zubindu",
    author_email="ansgar.wehrhahn@physics.uu.se",
    url="https://github.com/awehrhahn/CRIRES-planning-tool",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Beta",
    ],
)
