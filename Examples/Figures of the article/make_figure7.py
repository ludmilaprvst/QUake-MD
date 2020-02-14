#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 14:32:46 2020

@author: PROVOLU
"""

import sys
sys.path.append('../Source')
import PlotEvtObject as peo
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches



def get_pdf_mag(nom_fichier):
    PDF = np.genfromtxt(nom_fichier,
                        skip_header=5,
                        dtype=[('Depth', 'f8'),
                               ('Mag', 'f8'),
                               ('Weight', 'f8')])
    return PDF


Evtfilename = 'Data/Evt2016b.test.txt'
Obsfilename = 'Data/Obs2016b.test.txt'
evid = 650009

# Loading intensity data
files = peo.FilesForPlot(Evtfilename, Obsfilename)
evt = peo.PlotEvt(files)
evt.build(evid)
IDPs = evt.Obsevid[evt.Obsevid['Iobs']>=1]
Obsbin = pd.read_csv('Sorties/650009_Hmax25km/IDP_binning.txt')

# Loading IPEs solutions and predictions
IPEs_results_25 = pd.read_csv('Sorties/650009_Hmax25km/All_IPEs_classical_results.txt')
IPEs_results_11 = pd.read_csv('Sorties/650009_Hmax11km/All_IPEs_classical_results.txt')
Depi_plot = np.logspace(-1, 3, 1000)

# Loading the PDF
PDF_25 = get_pdf_mag('Sorties/650009_Hmax25km/HM.txt')
PDF_11 = get_pdf_mag('Sorties/650009_Hmax11km/HM.txt')

xi,yi = np.linspace(1, 25, 100), np.linspace(5.5, 8, 100)
xi,yi = np.meshgrid(xi, yi)
zi_25 = ml.griddata(PDF_25['Depth'], PDF_25['Mag'], PDF_25['Weight'], xi, yi, interp='linear')
zi_11 = ml.griddata(PDF_11['Depth'], PDF_11['Mag'], PDF_11['Weight'], xi, yi, interp='linear')

# Plot the figure 5
axis_font = {'fontname' : 'Arial', 'size': '12'}
fig_intensity = plt.figure(figsize=(10, 6))
gs = mpl.gridspec.GridSpec(4, 3, width_ratios=[1.2, 1, 0.1],
                                 height_ratios=[0.1, 0.1, 0.8, 0.5])
ax_PDF = fig_intensity.add_subplot(gs[2:, 0])
ax_Ifit = fig_intensity.add_subplot(gs[:3, 1])
ax_HM = fig_intensity.add_subplot(gs[3, 1])
ax_cb25 = fig_intensity.add_subplot(gs[:3, 2])
ax_cbPDF = fig_intensity.add_subplot(gs[1, 0])
ax_cbPDF2 = fig_intensity.add_subplot(gs[0, 0])


ax_PDF.imshow(zi_25, vmin=PDF_25['Weight'].min(), vmax=PDF_25['Weight'].max(), origin='lower', 
                              extent=[1, 25, 5.5, 8], aspect='auto',
                              interpolation='nearest', cmap=plt.cm.get_cmap('winter_r'))
ax_PDF.imshow(zi_11, vmin=PDF_11['Weight'].min(), vmax=PDF_11['Weight'].max(), origin='lower', 
                              extent=[1, 25, 5.5, 8], aspect='auto',
                              interpolation='nearest', cmap=plt.cm.get_cmap('copper_r'))
ax_PDF.set_xlabel('Depth [km]', **axis_font)
ax_PDF.set_ylabel('Magnitude', **axis_font)

ax_Ifit.fill_between([0.1, 1000], 7.5, 9.5, color='PaleVioletRed', alpha=0.1)
ax_Ifit.semilogx(IDPs['Depi'], IDPs['Iobs'], '.', color='Grey')
ax_Ifit.semilogx(Obsbin['Depi'].values, Obsbin['I'].values,'d', color=None,
                 markeredgecolor='k')
ax_Ifit.errorbar(Obsbin['Depi'], Obsbin['I'].values, yerr=Obsbin['StdI'].values,
                                fmt='d', color='k', mec='k')

normcb25 = mpl.colors.Normalize(vmin=1, vmax=25)
cmapcb25 = mpl.cm.get_cmap('viridis_r')


normcbPDF = mpl.colors.Normalize(vmin=PDF_25['Weight'].min(), vmax=PDF_25['Weight'].max())
cmapcbPDF = mpl.cm.get_cmap('winter_r')
cmapcbPDF2  = mpl.cm.get_cmap('copper_r')



for ind, row in IPEs_results_11.iterrows():
    Hypo = np.sqrt(Depi_plot**2 + row['H']**2)
    Ipred = row['C1'] + row['C2']*row['Mag']+row['Beta']*np.log10(Hypo) +row['Gamma']*Hypo
    couleur_depth = cmapcb25(normcb25(row['H']))
    ax_Ifit.semilogx(Depi_plot, Ipred, '--', color=couleur_depth,  alpha=0.5)
    ax_HM.plot(row['H'], row['Mag'], '^', color=couleur_depth)
for ind, row in IPEs_results_25.iterrows():
    Hypo = np.sqrt(Depi_plot**2 + row['H']**2)
    Ipred = row['C1'] + row['C2']*row['Mag']+row['Beta']*np.log10(Hypo) +row['Gamma']*Hypo
    couleur_depth = cmapcb25(normcb25(row['H']))
    ax_Ifit.semilogx(Depi_plot, Ipred, color=couleur_depth, alpha=0.6)
    ax_HM.plot(row['H'], row['Mag'], 'o', color=couleur_depth)


   
ax_Ifit.set_xlim([1, 500])
ax_Ifit.set_ylim([2, 10])
ax_Ifit.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
locmaj = mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2, 0.5, 1, 2, 5, 10 ))
ax_Ifit.xaxis.set_major_locator(locmaj)
ax_Ifit.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax_Ifit.set_xlabel('Epicentral distance [km]', **axis_font)
ax_Ifit.set_ylabel('Intensity', **axis_font)

ax_HM.set_xlim([0, 25+2]) 
ax_HM.set_ylim([5.5, 8]) 
ax_HM.set_xlabel('Depth [km]', **axis_font)
ax_HM.set_ylabel('Magnitude', **axis_font)
ax_HM.grid(which='both')
plt.tight_layout(pad=0.1, w_pad=0.5, h_pad=0)

cb = mpl.colorbar.ColorbarBase(ax_cb25, cmap=cmapcb25, norm=normcb25, orientation='vertical')
cb.set_label('Depth [km]', size=12)
cb.set_ticks(np.arange(0, 25+5, 5))
cb.ax.invert_yaxis()

cbPDF = mpl.colorbar.ColorbarBase(ax_cbPDF, cmap=cmapcbPDF, norm=normcbPDF, orientation='horizontal')
cbPDF.set_ticks([])
cbPDF.set_label('Space of solution with maximal depth 25 km')
ax_cbPDF.text(0.035, 0.03, 'Increasing weight', ha='center',
          va='center', rotation=0, color='White', weight='bold')
cbPDF2 = mpl.colorbar.ColorbarBase(ax_cbPDF2, cmap=cmapcbPDF2, norm=normcbPDF, orientation='horizontal')
cbPDF2.set_ticks([])
cbPDF2.set_label('Space of solution with maximal depth 11 km')
ax_cbPDF2.text(0.035, 0.03, 'Increasing weight', ha='center',
          va='center', rotation=0, color='White', weight='bold')
Io_uncertainties = mpatches.Patch(color='PaleVioletRed', alpha=0.3, label='I0 uncertainties')
Ibin_plt, = ax_Ifit.plot([], [], 'd', markerfacecolor='k', markeredgecolor='k', label='Binned intensity with RAVG method')
Iobs_plt, = ax_Ifit.plot([], [], '.', color='Gray', label='Observed intensities')
ax_Ifit.legend(handles=[Io_uncertainties, Ibin_plt, Iobs_plt], loc=4)

plt.savefig('Figure7.png', dpi=200, bbox_inches='tight')



