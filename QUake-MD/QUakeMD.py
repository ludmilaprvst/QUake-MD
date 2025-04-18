#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 12:07:33 2019

@author: baize-funck-ame
"""
# Importation des modules internes à Python

import queue as que
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
#import matplotlib.mlab as ml
from scipy.interpolate import griddata
import os.path
import copy
import sys
from mpl_toolkits.mplot3d import axes3d
import datetime
try:
    import mpl_toolkits.basemap.pyproj as pyproj
except:
    import pyproj
    
print('Import of the general python modules: success')
import library as ale
print('Import of the library module: success')
import PlotEvtObject as peo
print('Import of the PlotEvt object modules: success')
#from Modules_QUakeMD import *
from Modules_QUakeMD import add_I02obsbin, SearchBestStartDepth, inversion_MHI0, read_empe2, weight_percentile
from Modules_QUakeMD import PDF
print('Import of the tool module: success')
import WLSIC
print('Import of the least square module: success')

Std ={'A':0.5,'B':0.5,'C':0.5,'E':0.750}

class QUakeMD():
    """
    Core of QUage-MD GUI: compute a space of M/H solutions based on IDPs and IPEs
    
    Described in [ref biblio]
    """
    def __init__(self, evts, Evtname, Obsname, 
                       listVarEq, listVarCoeff, listeIbin, Parametername="", Ic=False, output_folder=False,
                       I0_option=True, imposed_depth=False, tag ='',
                       depth_min_ref=False, depth_max_ref=False,
                       LimitForSamplingInStd=2):
        
        self.evts = evts
        if not (len(listVarEq) == len(listVarCoeff)):
            return
        self.listVarEq = listVarEq
        self.listVarCoeff = listVarCoeff
        self.listeIbin = listeIbin
        self.Ic = Ic
        self.output_folder = output_folder
        # Create an outut folder if not given
        if not self.output_folder:
            self.output_folder = '{date:%Y-%m-%d}'.format(date=datetime.datetime.now())
            print('Output foldername:')
            print(self.output_folder)
            if not os.path.isdir(self.output_folder +'/'):
                os.makedirs(self.output_folder + '/')
        
        self.I0_option = I0_option
        self.imposed_depth = imposed_depth
        self.tag = tag
        if depth_min_ref == False:
            self.depth_min_ref = 1
        else:
            self.depth_min_ref = depth_min_ref
        if depth_max_ref == False:
            self.depth_max_ref = 25
        else:
            self.depth_max_ref = depth_max_ref
        self.LimitForSamplingInStd = LimitForSamplingInStd
         
        # Initialization of the logfile
        logfilename = '{date:%Y-%m-%d_%H-%M-%S.txt}'.format(date=datetime.datetime.now())
        self.logfile = open(self.output_folder + '/' + logfilename,'w')
        self.writeOnLogFile("Output_folder: " + self.output_folder)
        self.writeOnLogFile("Number of events: " + str(self.evts.qsize()))
        self.writeOnLogFile("Event File: " + Evtname)
        self.writeOnLogFile("Observation File: " + Obsname)
        #self.writeOnLogFile("Parameter File: " + Parametername)
        
        
        self.nbre_prof_test = 25
        self.NSamples = 20
        
        # Initialization of PDF parameters
        
        
        #Definition standard deviations for I0 based on quality factor
        self.Std = {'A': 0.5, 'B': 0.5, 'C': 0.5, 'E': 0.750}
        
        # Initialization of python outputs (for future)
        self.index_result = 0
        self.result_by_EMPE = pd.DataFrame(columns=['NumEvt','Bin_method', 'C1', 'C2', 'Beta',
                                               'Gamma', 'Mag', 'StdM', 'H', 'StdH', 'Io'])
    
        while self.evts.qsize() > 0:
            evt = self.evts.get()
            self.algorithm_QUakeMD(evt)
        self.logfile.close()
        
    def writeOnLogFile(self, s):
        self.logfile.write(s + '\n')
        print(s)
        
    
    def save_isoseist(self, evt, Ic):
        unique_metbin = np.unique(self.listeIbin)
        for metbin in unique_metbin:
            evt.Binning_Obs(1, Ic, method_bin=metbin)
            ObsBin = evt.ObsBinn
            ObsBin_save = copy.deepcopy(ObsBin)
            ObsBin_save = ObsBin_save[['EVID', 'Depi', 'I','StdI', 'StdLogR', 'Ndata']]
            ObsBin_save.to_csv(self.output_folder+'/'+str(evt.evid)+'/IDP_binning_' + metbin + '.txt')
            
    def cleaning_Ieq0(self, evt):
        evt.Obsevid.loc[evt.Obsevid['Iobs'] == 0, 'Iobs'] = 1
        evt.Obsevid.loc[evt.Obsevid['Iobs'] == 0, 'QIobs'] = 'C'
        evt.Obsevid = evt.Obsevid[evt.Obsevid['Iobs'] > 0]
        
    def init_figintensity(self):
        # Initialization of figure with results of inversion of intensity data
        fig_intensity = plt.figure(figsize=(6, 7))
        gs = mpl.gridspec.GridSpec(2, 2, width_ratios=[1, 0.1],
                               height_ratios=[1, 0.5])
        # Creation of subplot of fig_intensity
        ax = fig_intensity.add_subplot(gs[0])
        axMH_IPE = fig_intensity.add_subplot(gs[2])
        axcb = fig_intensity.add_subplot(gs[1])
        ax.set_xlabel('Epicentral distance [km]', size=12)
        ax.set_ylabel('Intensity', size=12)
        axMH_IPE.set_xlabel('Depth [km]', size=12)
        axMH_IPE.set_ylabel('Magnitude [Mw]', size=12)
        axMH_IPE.set_ylim([3, 7.5])
        return fig_intensity, ax, axMH_IPE, axcb

    def plot_figintensity(self, ax, axMH_IPE, axcb, evt, DataObs, depth, mag, I0, 
                          methode_bin, empe, C1, C2, Beta, gamma, comptlaw):
        depth_max_cb_plot = 25
        depth_min_cb_plot = 1
        normcb = mpl.colors.Normalize(vmin=self.depth_min_ref, vmax=self.depth_max_ref)
        normcb = mpl.colors.Normalize(vmin=depth_min_cb_plot, vmax=depth_max_cb_plot)
        cmapcb = mpl.cm.get_cmap('viridis_r')
        couleur_depth = cmapcb(normcb(depth))
        maxdepi = DataObs['Depi'].max()
        
        evt.Binning_Obs(depth, evt.Ic, method_bin=methode_bin)
        ObsBin = evt.ObsBinn
        
        titre_subplt = os.path.basename(empe)[:-4]
        
        axMH_IPE.plot(depth, mag, 'o', color = couleur_depth)
        ax.semilogx(DataObs['Depi'].values, DataObs['Iobs'].values,
                    '.', color='Gray', markerfacecolor='None',
                    markersize=2)
        Depi_pour_plot = np.logspace(-1, 3, 1000)
        Hypo_pour_plot = np.sqrt(Depi_pour_plot**2 + depth**2)
        
        Ipred = C1 + C2*mag + Beta*np.log10(Hypo_pour_plot) + gamma*Hypo_pour_plot
        Depi_bin = np.sqrt(ObsBin['Hypo']**2 - depth**2)
        
        ax.semilogx(Depi_bin, ObsBin['I'].values,
                    'd', color = couleur_depth, markeredgecolor='k')
        
        ax.errorbar(Depi_bin, ObsBin['I'].values, yerr=ObsBin['StdI'].values,
                    fmt='d', color=couleur_depth, mec='k')
        self.writeOnLogFile('I0 inversion, I0 predit')
        self.writeOnLogFile(str(I0) + ', '+str( C1 + C2*mag + Beta*np.log10(depth) + gamma*depth))
        ax.semilogx(0.1, I0, 's', color='Red')

        if comptlaw == 1:
            ax.fill_between([0.1, 1000], evt.Io_inf, evt.Io_sup,
                            color='PaleVioletRed', alpha=0.1)

        ax.semilogx(Depi_pour_plot, Ipred, color=couleur_depth)
        
        ax.set_title(titre_subplt)
        ax.set_xlim([1, maxdepi+100])
        ax.set_ylim([2, evt.Io_sup+1])
       
        axMH_IPE.set_title(titre_subplt)
        axMH_IPE.set_xlim([self.depth_min_ref-1, self.depth_max_ref+2])
        axMH_IPE.set_xlim([depth_min_cb_plot-1, depth_max_cb_plot+2])
        
        
        cb = mpl.colorbar.ColorbarBase(axcb,
                                       cmap=cmapcb,
                                       norm=normcb,
                                       orientation='vertical')
        cb.set_label('Depth [km]', size=12)
        cb.set_ticks(np.arange(self.depth_min_ref, self.depth_max_ref + 5, 5))
        cb.ax.invert_yaxis()
        ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        locmaj = mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2, 0.5, 1, 2, 5, 10 ))
        ax.xaxis.set_major_locator(locmaj)
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        return ax, axMH_IPE, axcb
    
    def togoodformat_Std(self, std):
        try:
            std = std[0]
        except TypeError:
            pass
        return std
    
    def fill_tripletsolution_byempe(self, evt, StdM_fin, StdH_fin, mag, depth,
                                    C1, C2, Beta, gamma, empe):
        # Initialization of the variable which will contain the space of solutions (after I0 filtering)
        Triplets = {'Magnitude': [],
                    'Profondeur': [],
                    'Io': [],
                    'Weights': []}
        # I0 filtering of the gaussian space of solutions
        #dM = self.LimitForSamplingInStd * StdM_fin / float(self.NSamples)
        dH = self.LimitForSamplingInStd * StdH_fin / float(self.NSamples)
        if dH == 0:
            dH = 0.25
        try:
            Mexplore = np.linspace(mag - self.LimitForSamplingInStd * StdM_fin, mag + self.LimitForSamplingInStd * StdM_fin,  self.NSamples)
            Hexp_min = max([self.depth_min_ref, depth - self.LimitForSamplingInStd * StdH_fin])
            Hexp_max = min([self.depth_max_ref, depth + self.LimitForSamplingInStd * StdH_fin])
            Hexplore = np.linspace(Hexp_min, Hexp_max, self.NSamples)
            
            sommewM = 0
            ntest_explore = len(Mexplore)*len(Hexplore)
            
            for Mm in Mexplore:
                wM = np.exp(-0.5*((Mm-mag)/StdM_fin)**2) 
                sommewM += wM
                sommewH = 0
                for Hm in Hexplore:                    
                    wH = np.exp(-0.5*((Hm-depth)/StdH_fin)**2)
                    sommewH += wH
                    if Hm > self.depth_max_ref:
                        Hm = self.depth_max_ref
                    if Hm < self.depth_min_ref:
                        Hm = self.depth_min_ref

                    Io_test = C1 + C2 * Mm + Beta * np.log10(Hm) + gamma * Hm
                    if (Io_test >= evt.Io_inf) and (Io_test <= evt.Io_sup):
                        Triplets['Magnitude'].append(float(Mm))
                        Triplets['Profondeur'].append(float(Hm))
                        Triplets['Io'].append(float(Io_test))
                        try:                        
                            Triplets['Weights'].append(wH[0]*wM[0])
                        except IndexError:
                            produit = wH*wM 
                            try:
                                Triplets['Weights'].append(produit[0])
                            except IndexError:
                                Triplets['Weights'].append(produit)
            # MAJ of logfile
            n_in_explore = len(Triplets['Magnitude'])
            pourcent_in = n_in_explore/ntest_explore*100.
            self.writeOnLogFile('Pourcent of solutions compatible with I0: ' +str(pourcent_in) +'%')
            
        except ZeroDivisionError:
            with open('No_start_depth'+ self.tag+'.txt','a') as no_startdepth:
                no_startdepth.write(str(evt.evid) + '\t' + empe + '\n')
        except:
            print('Unknown error in solutions space constitution')
        return Triplets
    

    def control_StdHvalue(self, depth, StdH_fin):
        if abs(depth - self.depth_max_ref) < 0.001:
            StdH_fin = max([StdH_fin, 5. / self.LimitForSamplingInStd])
        if abs(depth - self.depth_min_ref) < 0.001:
            StdH_fin = max([StdH_fin, 1. / self.LimitForSamplingInStd])
        return StdH_fin
    
    def algorithm_QUakeMD(self, evt):
        self.writeOnLogFile("\n")
        self.writeOnLogFile("Id of the event : " + str(evt.evid))
        self.writeOnLogFile("StdI_0 = " + str(evt.QI0))
        self.writeOnLogFile("Minimal depth limit: " + str(self.depth_min_ref))
        self.writeOnLogFile("Maximal depth limit: " + str(self.depth_max_ref))
        self.writeOnLogFile("Sigma sampling :" + str(self.LimitForSamplingInStd))
        
        # # Initialization of PDF
        self.all_PDF = PDF(self.depth_min_ref, self.depth_max_ref)
        
        # Initialization of barycenter
        (big_somme,  Barycentre_Mag, Barycentre_LogH,
         Barycentre_Io, poids_manquants, poids_presents)= np.zeros(6)
    
        # Creation of output folder for PDF and figures by event
        foldername = self.output_folder + '/' + str(evt.evid)
        if not os.path.isdir(foldername):
            os.makedirs(foldername)
        
        self.writeOnLogFile("Output folder for individual results created")

        DataObs = copy.deepcopy(evt.Obsevid)
        DataObs_ref = evt.Obsevid.copy()
        
        # Dealing with the IDP 0 values
        self.cleaning_Ieq0(evt)
        # Dealing with the IDP with I value superior to epicentral intensity
        I0_ini = evt.Io_ini
        evt.Obsevid.loc[DataObs_ref['Iobs'] >= I0_ini, 'Iobs'] = I0_ini
        # Saving the different isoseismal
        if self.Ic == False:
            Ic_ref = copy.deepcopy(evt.Ic)
        else:
            Ic_ref = self.Ic
        self.save_isoseist(evt, Ic_ref)
        # Application of the different IPEs stored in .txt files (for loop on the .txt files)
        for index in range(len(self.listVarEq)):
            Poids_branche = self.listVarCoeff[index]
            methode_bin = self.listeIbin[index]
            empe = self.listVarEq[index]
            
            # Initialization of Io and Ic
            StdI0 = evt.QI0
            QI0_inv = StdI0
            I0 = I0_ini
            Ic = Ic_ref
            
            # Preparation of the data
            evt.Binning_Obs(8, Ic, method_bin=methode_bin) # Abritary 8 km depth
            #print(evt.ObsBinn)

            
            # Initialization of figure with results of inversion of intensity data
            fig_intensity, ax, axMH_IPE, axcb = self.init_figintensity()

            # Initialization of variables
            (nloi, comptlaw, sumpoidsEMPE) =  np.zeros(3)
            # Initialization of barycenter
            (Barycentre_MagLaw, Barycentre_LogHLaw, Barycentre_IoLaw) =  np.zeros(3)
            # Weight for barycenter calculus
            (poids_manquants_law, poids_presents_law) = np.zeros(2)
            # Initialization of PDF
            self.all_PDF.initPDFlaw()
            
            self.writeOnLogFile('Ic = ' + str(Ic))
            inversion = 'BasedOnIobs'
                
            # Update logfile
            self.writeOnLogFile('IPE: ' + str(empe))
            # Loading the IPE .txt files
            try:
                beta_liste, c1_liste, c2_liste, poidsEMPE_liste, gamma_liste = read_empe2(empe)
            except:
                self.writeOnLogFile('Error reading ' + str(empe))
            
            # Computing the gaussian space of solution for each IPE
            for Beta, C1, C2, gamma, w_empe in zip(beta_liste, c1_liste, c2_liste, gamma_liste, poidsEMPE_liste):
                Singular = False
                nloi += 1
                comptlaw += 1
                # Weight control
                sumpoidsEMPE += w_empe
                
                # Initialization of the inverted parameters
                I0_ini = evt.Io_ini
                #evt.Io_inv = I0_ini
                evt.add_invMHI0_variables(QI0_inv, I0_ini, Ic)
                print()
                start_depth, start_mag = SearchBestStartDepth(evt, Beta, C1, C2, gamma,
                                                              self.depth_min_ref, self.depth_max_ref, self.nbre_prof_test)
                #print(start_depth, start_mag)
                if not start_depth:
                    Singular = True
                    start_depth = 10
                    start_mag = (I0_ini - (C1 + Beta * np.log10(start_depth) + gamma * start_depth)) / C2
                    #inversion = 'BasedOnIo'
                    self.writeOnLogFile("No start depth, change inversion")
                
                # Update logfile
                self.writeOnLogFile('Inversion: ' + str(inversion))
                self.writeOnLogFile('I = ' + str(C1) + ' + ' + str(C2) + 'M ' + str(Beta) + 'log10(Dhypo) + ' + str(gamma) + 'Dhypo')
                # Inversion of M and H (in Modules_QUakeMD)
                (mag, depth, I0, StdM_fin,
                 StdH_fin, evt) = inversion_MHI0(evt, 
                                                 C1, C2, Beta, gamma,
                                                 start_depth, start_mag, 
                                                 self.depth_min_ref, self.depth_max_ref,
                                                 self.imposed_depth, self.I0_option) 
                StdM_fin = self.togoodformat_Std(StdM_fin)
                StdH_fin = self.togoodformat_Std(StdH_fin)
                
                self.index_result += 1
                # Update logfile
                self.writeOnLogFile('M = %.2f; H = %.2f; I0=%.f' %(mag, depth, I0))
                self.writeOnLogFile('StdM = %.2f; StdH = %.2f' %(StdM_fin, StdH_fin))
                
                # Control figure of the IPE fit to the intensity decrease
                ax, axMH_IPE, axcb = self.plot_figintensity(ax, axMH_IPE, axcb,  evt, DataObs, depth, mag, I0, 
                                      methode_bin, empe, C1, C2, Beta, gamma, comptlaw)
                
                # Avoiding too small values of H sigma
                StdH_fin = self.control_StdHvalue(depth, StdH_fin)
                # Storage of the IPEs central values and associated sigmas
                self.result_by_EMPE.loc[self.index_result] = [evt.evid, methode_bin, C1, C2, Beta, gamma, mag,
                                   StdM_fin, depth, StdH_fin, I0]
                Triplets = self.fill_tripletsolution_byempe(evt, StdM_fin, StdH_fin, mag, depth,
                                                C1, C2, Beta, gamma, empe)
                
                                 
#%%
                if len(Triplets['Magnitude']) == 0:
                    self.writeOnLogFile('I='+str(C1)+'+'+str(C2)+'M'+str(Beta)+'log(hypo)')
                    poids_manquants_law+=w_empe
                    poids_manquants+=w_empe*Poids_branche
                    self.writeOnLogFile('LimitForSamplingInStd')
                    self.writeOnLogFile(str(self.LimitForSamplingInStd))
                    self.writeOnLogFile('No solutions compatible with I0 found for this equation\n')    
                else:
                    # Normalization of the weights
                    Triplets['Weights'] = np.array(Triplets['Weights'])
                    Triplets['Weights'] = Triplets['Weights']/Triplets['Weights'].sum()
                    
                    poids_presents+=w_empe*Poids_branche
                    poids_presents_law+=w_empe
                    
                    # Attributing IPE's weight
                    Triplets['Weights'] = Triplets['Weights']*w_empe
                    # Barycenter computation for the IPEs in the .txt file
                    Barycentre_MagLaw += np.sum(np.array(Triplets['Magnitude'])*np.array(Triplets['Weights']))
                    Barycentre_LogHLaw += np.sum(np.log10(np.array(Triplets['Profondeur']))*np.array(Triplets['Weights']))    
                    Barycentre_IoLaw += np.sum(np.array(Triplets['Io'])*np.array(Triplets['Weights']))
                    
                    # Attributing IPE's .txt file  weight
                    Triplets['Weights'] = Triplets['Weights']*Poids_branche

                    # Barycenter computation
                    Barycentre_Mag += np.sum(np.array(Triplets['Magnitude'])*np.array(Triplets['Weights']))
                    Barycentre_LogH += np.sum(np.log10(np.array(Triplets['Profondeur']))*np.array(Triplets['Weights']))    
                    Barycentre_Io += np.sum(np.array(Triplets['Io'])*np.array(Triplets['Weights']))
                    # Storage of the space of solutions in a matrix for all IPEs and for each .txt files
                    self.all_PDF.update_PDF(Triplets, Poids_branche)
                        

            # Write PDF file by law
            try:
                # Control on barycenter 
                Barycentre_MagLaw = Barycentre_MagLaw/poids_presents_law
                Barycentre_LogHLaw = Barycentre_LogHLaw/poids_presents_law
                Barycentre_IoLaw = Barycentre_IoLaw/poids_presents_law
                self.all_PDF.PDF_HMLaw = self.all_PDF.PDF_HMLaw/self.all_PDF.PDF_HMLaw.sum()
                # MAJ du logfile
                self.writeOnLogFile('Barycenter of the group of IPE:')
                self.writeOnLogFile('M = %.2f; H = %.2f; I0 = %.2f' % (Barycentre_MagLaw, 10**Barycentre_LogHLaw, Barycentre_IoLaw))
                #print (str(evt.year))
                
                # End of the control figure of the IPE fit to the intensity decrease
                plt.tight_layout()      
                Io_uncertainties = mpatches.Patch(color='PaleVioletRed', alpha=0.3, label='I0 uncertainties')
                Ibin_plt, = ax.plot([], [], 'd', markerfacecolor='w', markeredgecolor='k', label='Binned intensity with'+ methode_bin +' method')
                Iobs_plt, = ax.plot([], [], '.', color='Gray', label='Observed intensities')
                ax.legend(handles=[Io_uncertainties, Ibin_plt, Iobs_plt])
                
                fig_intensity.savefig(foldername+'/'+str(evt.evid)+'_fit_intensity_Law_'+str(index)+'_'+ methode_bin+'.jpeg', dpi=100,
                                        bbox_inches='tight')
                # Save the IPE's .txt group space of solutions
                fileLaw = open(self.output_folder+'/'+str(evt.evid)+'/Law_'+str(index)+'_'+ methode_bin+ '_HM.txt', 'w')
                self.write_header_PDF(fileLaw, evt, Barycentre_IoLaw, Barycentre_MagLaw, Barycentre_LogHLaw)
                self.write_HMPDF_files(fileLaw, self.all_PDF.PDF_HMLaw)
                fileLaw.close()
            except ZeroDivisionError:
                self.all_PDF.PDF_HMLaw = self.all_PDF.PDF_HMLaw/self.all_PDF.PDF_HMLaw.sum()
                self.writeOnLogFile('Line of protocole:')
                self.writeOnLogFile(str(Poids_branche) +", " + str(methode_bin))
                self.writeOnLogFile('could not find an Io compatible result')
                self.writeOnLogFile(str(Barycentre_MagLaw) + ", "+str(10**Barycentre_LogHLaw)+ ", "+str(Barycentre_IoLaw))
              
        # Save the space of solutions
        try:
            # Control on barycenter
            if poids_manquants!=0:
                Barycentre_Mag = Barycentre_Mag/poids_presents
                Barycentre_LogH = Barycentre_LogH/poids_presents
                Barycentre_Io = Barycentre_Io/poids_presents 
                self.all_PDF.normalize_PDF()
        except ZeroDivisionError:
            self.writeOnLogFile('No result compatible with I0 could be found')
            self.writeOnLogFile('Last solution computed:')
            self.writeOnLogFile(str(mag)+ ", "+str(depth)+ ", "+str(I0))
            self.writeOnLogFile(str(evt.Io_inf)+ ", "+str(evt.Io_sup))
            with open('Noresult.txt', 'a') as no_result:
                no_result.write(str(evt.evid)+'\t'+str(mag)+'\t'+str(depth)+'\t'+str(I0)+'\n')
            return
        # MAJ du logfile
        self.writeOnLogFile('Final barycenter :\n')
        self.writeOnLogFile('M=%.2f;H=%.2f;I0=%.2f\n' % (Barycentre_Mag, 10**Barycentre_LogH, Barycentre_Io))
        # Finalisation et sauvegarde de la figure de controle des solutions avec les intensites
       
        # Save and plot space of solution in HM space
        filePDF = open(self.output_folder+'/'+str(evt.evid)+'/HM.txt','w')
        self.write_header_PDF(filePDF, evt, Barycentre_Io, Barycentre_Mag, Barycentre_LogH)
        prof_plot, mag_plot, poids_plot = self.write_HMPDF_files(filePDF, self.all_PDF.PDF_HM)
        filePDF.close()
        self.plot_HM(evt, 1, 25, prof_plot, mag_plot, poids_plot,
                     Barycentre_Mag, empe)
        
        # Weighted percentiles
        poids_plot = np.array(poids_plot)
        percentile16_prof = weight_percentile(prof_plot, poids_plot, 0.16)
        percentile84_prof = weight_percentile(prof_plot, poids_plot, 0.84)

        percentile16_mag = weight_percentile(mag_plot, poids_plot, 0.16)
        percentile84_mag = weight_percentile(mag_plot, poids_plot, 0.84)
        
        # Enregistrement et plot de la PDF HMIo
        filePDF = open(self.output_folder+'/'+str(evt.evid)+'/HMIo.txt','w')
        self.write_header_PDF(filePDF, evt, Barycentre_Io, Barycentre_Mag, Barycentre_LogH)
        prof_plot, mag_plot, io_plot, poids_plot = self.write_HMI0_files(filePDF, self.all_PDF.PDF_HMIo)
        filePDF.close()
        self.plot_HMI0(evt, prof_plot, mag_plot, io_plot, poids_plot)
        
        # Plot et sauvegarde de la PDF HIo 
        filePDF = open(self.output_folder+'/'+str(evt.evid)+'/HIo.txt','w')
        self.write_header_PDF(filePDF, evt, Barycentre_Io, Barycentre_Mag, Barycentre_LogH)
        prof_plot, io_plot, poids_plot = self.write_HI0_files(filePDF, self.all_PDF.PDF_HIo)
        filePDF.close()
        self.plot_HI0(evt, io_plot, prof_plot, poids_plot, 1, 25)
        
        # Calcul des percentiles
        poids_plot = np.array(poids_plot)
        percentile16_io = weight_percentile(io_plot, poids_plot, 0.16)
        percentile84_io = weight_percentile(io_plot, poids_plot, 0.84)

        if not os.path.isfile(self.output_folder+'/'+'file_temp_'+ self.tag +'.txt'):
            file_temp = open(self.output_folder+'/'+'file_temp_'+ self.tag +'.txt', 'w')
            file_temp.write('EVID\tI0 cat.\tQI0\tIc\tMbary\tM16th\tM84th\tHbary\tH16th\tH84th\tI0bary\tI016th\tI084th\n')
            file_temp.close()
        with open(self.output_folder+'/'+'file_temp_'+ self.tag +'.txt', 'a') as file_temp:
            file_temp.write('%d\t%0.1f\t%s\t%0.1f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n' 
                            % (evt.evid, evt.Io_ini, evt.QI0name, evt.Ic,
                                 Barycentre_Mag, percentile16_mag, percentile84_mag,
                                 10**Barycentre_LogH, percentile16_prof, percentile84_prof,
                                 Barycentre_Io, percentile16_io, percentile84_io))
        
        self.result_by_EMPE.to_csv(self.output_folder+'/'+str(evt.evid)+'/All_IPEs_classical_results.txt')
        return self.result_by_EMPE
    
    def plot_HM(self, evt, depth_min, depth_max,
                prof_plot, mag_plot, poids_plot,
                Barycentre_Mag, empe):
        
        # Plot MH space of solutions
        # Preparation des donnees pour le plot (gridding)
        mag_lim_min = min(mag_plot)-0.5
        mag_lim_max = max(mag_plot)+0.75
        depth_min = 0.2
        depth_max = 25
        #Normalisation pour plot
        poids_plot = poids_plot/max(poids_plot)
        xi, yi = np.mgrid[depth_min:depth_max:100j, mag_lim_min:mag_lim_max:100j]
        points_HM = []
        for xx in range(len(prof_plot)):
            points_HM.append([prof_plot[xx], mag_plot[xx]])
        points_HM = np.array(points_HM)

        try:
            #zi = ml.griddata(prof_plot, mag_plot, poids_plot, xi, yi, interp='linear')
            zi = griddata(points_HM, poids_plot, (xi, yi), method='linear')
            not_plot = False
        except RuntimeError:
            self.writeOnLogFile('Magnitude too small for the plot:')
            self.writeOnLogFile(str(Barycentre_Mag))
            self.writeOnLogFile('Number of points for the grid:')
            self.writeOnLogFile(str(len(mag_plot))+ ", "+str(len(prof_plot)))
            self.writeOnLogFile(str(min(prof_plot))+ ", "+str(max(prof_plot))+ ", "+str(depth_min)+ ", "+str(depth_max))
            self.writeOnLogFile(str(min(mag_plot))+ ", "+str(max(mag_plot)))
            self.writeOnLogFile('NumEvt:' + str(evt.evid))
            self.writeOnLogFile(str(empe))
            not_plot = True
            return
        except AttributeError:
            not_plot = True
            print("Mettre a jour pour la nouvelle version de matplotlib")
            return
           
        plt.figure(figsize=(5,5))

        plt.imshow(zi.T, vmin=poids_plot.min(), vmax=poids_plot.max(), origin='lower', 
                   extent=[depth_min, depth_max, mag_lim_min, mag_lim_max], aspect='auto',
                   interpolation='nearest', cmap=plt.cm.get_cmap('winter_r'))
        cbar = plt.colorbar()
        cbar.ax.text(.6, (2*1.5+1)/8., 'Increasing weight', ha='center', va='center', rotation=90,
                     color='White', weight='bold')
        cbar.set_ticks([])
        plt.xlabel('Depth [km]')
        plt.ylabel('Magnitude')
        plt.title(str(evt.evid))

        plt.savefig(self.output_folder+'/'+str(evt.evid)+'/HM.png')  
        plt.show()
        
    def plot_HMI0(self, evt, prof_plot, mag_plot, io_plot, poids_plot):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(prof_plot, mag_plot, io_plot, c=poids_plot, cmap='winter_r')
        ax.set_xlim([1, 25])
        ax.set_ylim([min(mag_plot)-1, max(mag_plot)+1])
        ax.set_zlim([4, 10])
        ax.set_xlabel('Depth [km]')
        ax.set_ylabel('Magnitude')
        ax.set_zlabel('Io')
        plt.savefig(self.output_folder+'/'+str(evt.evid)+'/HMIo.png') 
        plt.show()  
        
    def plot_HI0(self, evt, io_plot, prof_plot, poids_plot,
                 depth_min, depth_max):
        # Figure PDF HIo
        plt.figure(figsize=(5,5))
        
        poids_plot = np.array(poids_plot)   
        #Normalisation
        poids_plot = poids_plot/max(poids_plot)
#        max_io_lim =max(io_plot)+2
#        min_io_lim =min(io_plot)-2
        
        xi, yi = np.mgrid[depth_min:depth_max:100j, 2:12:100j]
        points_HI0 = []
        for xx in range(len(prof_plot)):
            points_HI0.append([prof_plot[xx], io_plot[xx]])
        points_HI0 = np.array(points_HI0)
        zi = griddata(points_HI0, poids_plot, (xi, yi), method='linear')

        plt.imshow(zi.T, vmin=poids_plot.min(), vmax=poids_plot.max(), origin='lower', 
                  extent=[depth_min, depth_max, 2, 12], aspect='auto',
                  interpolation='nearest', cmap=plt.cm.get_cmap('winter_r'))
        cbar = plt.colorbar()
       
        cbar.ax.text(.6, (2*1.5+1)/8., 'Increasing weight', ha='center', va='center', rotation=90,
                     color='White', weight='bold')
        cbar.set_ticks([])
        plt.xlabel('Depth [km]')
        plt.ylabel('Io')
        plt.title(str(evt.evid))
        plt.savefig(self.output_folder+'/'+str(evt.evid)+'/HIo.png')  
        plt.show()
        
    def write_header_PDF(self, fileLaw, evt, Barycentre_Io, Barycentre_Mag, Barycentre_LogH):
        fileLaw.write('NumEvt: '+ str(evt.evid) + ', year=' + str(int(evt.year)))
        fileLaw.write(', I0 from catalogue = '+ str(round(evt.Io_ini, 1)) +'\n')
        fileLaw.write('Barycenter Io:' + str(round(Barycentre_Io,2))+'\n')
        fileLaw.write('Barycenter M:' + str(round(Barycentre_Mag,2))+'\n')
        fileLaw.write('Barycenter H:' + str(round(10**Barycentre_LogH,2))+'\n')
        
    def write_HMPDF_files(self, filePDF, PDF_HM):
        filePDF.write('H[km]\tMag\tPDF\n')
        prof_plot = []
        mag_plot = []
        poids_plot = []
        for jj in range(PDF_HM.shape[1]):    
            for ii in range(PDF_HM.shape[0]):
                magni = self.all_PDF.EchelleMRef[jj]
                h = self.all_PDF.EchelleLogHRef[ii]
                if PDF_HM[ii][jj]>10**-6:
                    prof_plot.append(h)
                    mag_plot.append(magni)
                    poids_plot.append(PDF_HM[ii][jj])
                    filePDF.write('%0.4f\t%0.4f\t%0.6f\n' % (h, magni, PDF_HM[ii][jj]))
        return prof_plot, mag_plot, poids_plot
    
    def write_HI0_files(self, filePDF, PDF_HIo):
        prof_plot = []
        io_plot = []
        poids_plot = []
        filePDF.write('H[km]\Io\tPDF\n')
        for ii in range(PDF_HIo.shape[0]):
            for jj in range(PDF_HIo.shape[1]):
                h = self.all_PDF.EchelleLogHRef[ii]
                Int0 = self.all_PDF.EchelleIoRef[jj]
                if PDF_HIo[ii][jj]>0.000001:
                    prof_plot.append(h)
                    io_plot.append(Int0)
                    poids_plot.append(PDF_HIo[ii][jj])
                    filePDF.write('%0.4f\t%0.4f\t%0.6f\n' % (h, Int0, PDF_HIo[ii][jj]))
        return prof_plot, io_plot, poids_plot
    
    def write_HMI0_files(self, filePDF, PDF_HMIo):
        filePDF.write('H[km]\tMag\tIo\tPDF\n')
        prof_plot = []
        mag_plot = []
        io_plot = []
        poids_plot = []
        for kk in range(PDF_HMIo.shape[2]):
            for jj in range(PDF_HMIo.shape[1]):    
                for ii in range(PDF_HMIo.shape[0]):
                    io = self.all_PDF.EchelleIoRef[kk]
                    magni = self.all_PDF.EchelleMRef[jj]
                    h = self.all_PDF.EchelleLogHRef[ii]
                    if PDF_HMIo[ii][jj][kk]>10**-6:
                        prof_plot.append(h)
                        mag_plot.append(magni)
                        io_plot.append(io)
                        poids_plot.append(PDF_HMIo[ii][jj][kk])
                        filePDF.write('%0.4f\t%0.4f\t%0.4f\t%0.6f\n' % (h, magni, io, PDF_HMIo[ii][jj][kk]))
        return prof_plot, mag_plot, io_plot, poids_plot
