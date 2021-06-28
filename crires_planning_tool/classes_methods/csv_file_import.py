#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:00:27 2020

Initializes .csv files with lists of planets, could be extended to also process other types of .csv list.

@author: jonaszbinden
"""

import os
import sys
from os.path import dirname, join

import pandas as pd

path = join(dirname(__file__), "../csv_files")
default_file = "PlanetList"

if len(sys.argv) > 1:
    file_name = sys.argv[1]
else:
    file_name = default_file


if len(file_name) == 0:
    file_name = default_file
    print("Default file " + default_file + " is being used")
elif os.path.exists(join(path, file_name + ".csv")) == False:
    message = f"Error:file {file_name} does not exist in current directory"
    raise FileNotFoundError(message)


def main():
    df = pd.read_csv(join(path, file_name + ".csv"))
    df_names = df["pl_name"]
    return df_names


if __name__ == "__main__":
    main()

