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

table_values = ['_______________________________________________',
				' ',
				r'$\mathit{S/S_0}$' + '      ' +   r'$\mathrm{CAM4}$' + ' ' +  r'$\mathrm{CO_2}$' + ' ' + r'$\mathrm{(ppmv)}$' + '      '+   r'$\mathrm{CAM5}$' + ' ' +  r'$\mathrm{CO_2}$' + ' ' + r'$\mathrm{(ppmv)}$',
				'_______________________________________________',

				'0.700                   130,000 ',
				'0.725                     92,000 ',
				'0.750                     65,500 ',
				'0.775                     45,119 ',
				'0.800                     31,688 ',
				'0.825                     20,809 ',
				'0.850                     13,562 ',
				'0.875                       8,692 ',
				'0.900                       5,254                            6,000',
				'0.925                       3,108                            2,700',
				'0.950                       1,616                            1,400',
				'0.975                          813                               700',
				'1.000                          368.9                            284.7',
				'1.025                          139.2                            120',
				'1.050                            51.4                              35',
				'1.075                            15.1  ',
				'1.100                              4.5  ',
				'_______________________________________________']


#create plot
fig = plt.figure(figsize=(3.46457, 4))

space = 0
for row in table_values:
	if row == table_values[3] or row == table_values[2]:
		fs = 10
		space = space - 0.01

	else:
		fs = 8

	plt.text(0, 1-space, row, fontsize=fs)

	space = space + 0.05

plt.axis('off')

plt.show()
fig.savefig("ED_figure1.eps", format='eps', bbox_inches='tight')

