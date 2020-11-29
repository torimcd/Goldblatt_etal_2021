#!/usr/bin/env python3

"""
Author: Victoria McDonald
email: vmcd@atmos.washington.edu
website: http://torimcd.github.com
license: BSD

"""
import matplotlib
matplotlib.use("Agg")

import numpy
import matplotlib
import matplotlib.pyplot as plt



# ------------------------------------------------------------------------
# This script plots a table of CO2 and S/So values run in each model
# ------------------------------------------------------------------------

table_values = [['0.700', '130000', ' '],
				['0.725', '92000', ' '],
				['0.750', '65500', ' '],
				['0.775', '45119', ' '],
				['0.800', '31688'. ' '],
				['0.825', '20809', ' '],
				['0.850', '13562', ' '],
				['0.875', '8692', ' '],
				['0.900', '5254', '6000'],
				['0.925', '3108', '2700'],
				['0.950', '1616', '1400'],
				['0.975', '813', '700'],
				['1.000', '368.9', '284.7'],
				['1.025', '139.2', '120'],
				['1.050', '51.4', '35'],
				['1.075', '15.1', ' '],
				['1.100', '4.5', ' '],]

row_labels = [r'$\mathsf{S/S_0}$', r'$\mathsf{CAM4 CO_2 (ppmv)}$', r'$\mathsf{CAM5 CO_2 (ppmv)}$']

#create plot
fig = plt.figure()

plt.table(table_values, row_labels=[])


plt.show()
fig.savefig("ED_figure1.pdf", bbox_inches='tight')

