#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 12:20:40 2021

@author: linnboldtc
"""
import pandas as pd
import plotly.express as px

# --- Import CSV data of given period
df = pd.read_csv('/Users/linnboldtc/CRIRES-planning-tool/crires_planning_tool/csv_files/Eclipse_events_processed_2020-12-01_15d.csv')

# --- Plotting the length of each planet's transit (from start to end time) across the given period
fig = px.timeline(df, x_start="time_begin", x_end="time_end", y="name", color="stellar_effective_temperature", title='Length of each planet transit')
fig.update_yaxes(autorange="reversed") # otherwise tasks are listed from the bottom up
fig.show()

# --- Since ranking/exposure time/transit time seems a bit wonky at the moment, this is showing what the SNR at each transit is as some indicator of how possible it is
fig = px.scatter(df, x = 'time_mid', y = 'snr_median', color="stellar_effective_temperature", text='name', title='SNR at each transit (mid) time')
fig.update_traces(textposition='top center')
fig.update_layout(hovermode="x unified")
fig.show()