#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:31:20 2018

@author: PROVOLU
"""
import numpy as np
import pandas as pd
import sys

from scipy.optimize import curve_fit
#import statsmodels.api as sm
import matplotlib.pyplot as plt


#import Fonctions_bin as a
#from Modules_Getbeta import CalcDist, read_obsfile, read_evtfile, read_critfile


class WLSIC():
    def __init__(self, ObsBinn, depth, Beta, Gamma, I0):
        self.beta = Beta
        self.Obsbin = ObsBinn
        self.depth = depth
        self.I0 = I0
        self.gamma = Gamma
        
    def EMIPE_H(self, Depi, H):
        I = self.I0 + self.beta*np.log10(np.sqrt(Depi**2+H**2)/H)+self.gamma*(np.sqrt(Depi**2+H**2)-H)
        return I

    def EMIPE_JACdH(self, Depi, H):
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        g = self.beta*(tmpValue**2-1.)/(H*np.log(10))+self.gamma*(tmpValue-1)
        return g.reshape(len(Depi),1)
    
    def EMIPE_I0(self, Depi, I0):
        I = I0 + self.beta*np.log10(np.sqrt(Depi**2+self.depth**2)/self.depth)+self.gamma*(np.sqrt(Depi**2+self.depth**2)-self.depth)
        return I

    def EMIPE_JACdI0(self, Depi, I0):
        g = np.ones(len(Depi))
        return g.reshape(len(Depi),1)

    def do_wlsic_depth(self, depth_inf, depth_sup):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0
        resH = curve_fit(self.EMIPE_H, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdH, bounds=(depth_inf, depth_sup),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH
    
    def do_wlsic_I0(self, I0_inf, I0_sup):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        
        resI0 = curve_fit(self.EMIPE_I0, Depi, Ibin, p0=self.I0,
                                  jac= self.EMIPE_JACdI0, bounds=(I0_inf, I0_sup),
                                  sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                  ftol=1e-2)
        return resI0

class WLSIC_M():
    def __init__(self, ObsBinn, depth, mag, Beta, Gamma, C1, C2):
        self.beta = Beta
        self.gamma = Gamma
        self.Obsbin = ObsBinn
        self.depth = depth
        self.C1 = C1
        self.C2 = C2
        self.mag = mag
        
    def EMIPE_HM(self, Depi, H):
        I = self.C1 + self.C2*self.mag+ self.beta*np.log10(np.sqrt(Depi**2+H**2))+self.gamma*np.sqrt(Depi**2+H**2)
        return I
    
    def EMIPE_M(self, Depi, M):
        I = self.C1 + self.C2*M+ self.beta*np.log10(np.sqrt(Depi**2+self.depth**2))+self.gamma*np.sqrt(Depi**2+self.depth**2)
        return I

    def EMIPE_JACdHM(self, Depi, H):                        
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        g = (tmpValue)*((self.beta/(np.log(10)*Hypo))+self.gamma)
        return g.reshape(len(Depi),1)

    def EMIPE_JACdM(self, Depi, H):                        
        g = self.C2*np.ones(len(Depi))
        return g.reshape(len(Depi),1)

    def do_wlsic_depthM(self, depth_inf, depth_sup):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resH = curve_fit(self.EMIPE_HM, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdHM, bounds=(depth_inf, depth_sup),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH
    
    def do_wlsic_depthM_std(self):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resH = curve_fit(self.EMIPE_HM, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdHM, bounds=(self.depth-0.01, self.depth+0.01),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH

    def do_wls_M(self):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resM = curve_fit(self.EMIPE_M, Depi, Ibin, p0=self.mag,
                                 jac= self.EMIPE_JACdM, 
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=1e-2)
        return resM
    
    def do_wls_M_std(self):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resM = curve_fit(self.EMIPE_M, Depi, Ibin, p0=self.mag, bounds=(self.mag-0.001, self.mag+0.001),
                                 jac= self.EMIPE_JACdM, 
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=1e-2)
        return resM


    
class WLS_SIGMA():
    def __init__(self, ObsBinn, depth, mag, Beta, Gamma, C1, C2, loss_fc):
        self.beta = Beta
        self.gamma = Gamma
        self.Obsbin = ObsBinn
        self.depth = depth
        self.C1 = C1
        self.C2 = C2
        self.mag = mag
        self.loss_fc = loss_fc

    def EMIPE_HM(self, Depi, M, H):
        I = self.C1 + self.C2*M+ self.beta*np.log10(np.sqrt(Depi**2+H**2))+self.gamma*np.sqrt(Depi**2+H**2)
        return I

    def EMIPE_JACdHM(self, Depi, M, H):                        
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        gH = (tmpValue)*((self.beta/(np.log(10)*Hypo))+self.gamma)
        gM = self.C2*np.ones(len(Depi))
        g = np.transpose([gM, gH])
        return g
    
    # a travailler
    def do_wls_HM(self):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resM = curve_fit(self.EMIPE_HM, Depi, Ibin, jac= self.EMIPE_JACdHM, 
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 bounds=([-np.inf, 0], [np.inf, np.inf]), loss=self.loss_fc,
                                 verbose=2,
                                 ftol=1e-4)
        return resM
    
class WLSIC_BH():
     def __init__(self, ObsBinn, depth, Beta, I0):
        self.beta = Beta
        self.Obsbin = ObsBinn
        self.depth = depth
        self.I0 = I0

     def EMIPE_HB(self, Depi, B, H):
        I = self.I0 + B*np.log10(np.sqrt(Depi**2+H**2)/H)
        return I

     def EMIPE_JACdHB(self, Depi, B, H):                        
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        gH = self.beta*(tmpValue**2-1.)/(H*np.log(10))
        gB = np.log10(np.sqrt(Depi**2+self.depth**2)/self.depth)
        g = np.transpose([gB, gH])
        return g
    
     # a travailler
     def do_wls_HB(self, depth_inf, depth_sup):
        Ibin = self.Obsbin['I'].values
        depi_carre = np.round(self.Obsbin['Hypo'].values**2 - self.depth**2,5)
        depi_carre = np.array(depi_carre)
        Depi = np.sqrt(depi_carre)
        ind_nan = np.isnan(Depi)
        Depi[ind_nan] = 0

        resM = curve_fit(self.EMIPE_HB, Depi, Ibin, jac= self.EMIPE_JACdHB, 
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 bounds=([-6, depth_inf,], [0, depth_sup]), 
                                 verbose=0,
                                 ftol=1e-5)
        return resM
    
"""
NumEvt = 880077
millesime = '2016'
binning = a.ROBS_c

ObsFileName = '../../../Data/Obs'+millesime+'.txt'
EvtFileName = '../../../Data/Evt'+millesime+'.txt'
nomFichierCritique = '../../../Data/IcDc.'+millesime+'.txt'


ObsFile = read_obsfile(ObsFileName)
EvtFile = read_evtfile(EvtFileName)
fichierCritique = read_critfile(nomFichierCritique)

EvtFile = EvtFile[EvtFile['EVID']==NumEvt]
ObsFile = ObsFile[ObsFile['EVID']==NumEvt]
Lon_evt = EvtFile['Lon'].values[0]
Lat_evt = EvtFile['Lat'].values[0]

ObsFile.loc[:,'Depi'] = ObsFile.apply(lambda row:CalcDist(row['Lon'],row['Lat'],Lon_evt,Lat_evt),axis=1) 

Std = {'I0':{'A':0.25,'B':0.375,'C':0.5,'E':0.75,'K':0.375},
           'Iobs':{b'A':0.5,b'B':0.577,b'C':0.710,b'D':1.0}}

Beta = -3.5
I0 = EvtFile['I0'].values[0]
Hini = 10

ObsFile.loc[:,'Hypo'] = ObsFile.apply(lambda row:np.sqrt(row['Depi']**2+Hini**2), axis=1) 
ObsFile.loc[:,'QIobs2'] = ObsFile.apply(lambda row:Std['Iobs'][row['QIobs']], axis=1)
Ic = {}
Dc = {}
Ic[NumEvt] = 3
Dc[NumEvt] = 130

ObsBin = binning(ObsFile['Iobs'].values,
                            ObsFile['Hypo'].values,
                             ObsFile['QIobs2'].values,
                             I0,0.25, 10, int(NumEvt), Ic[NumEvt], Dc[NumEvt], 30)
colonnes = ['EVID', 'Hypo', 'I', 'Io', 'QIo', 'StdI', 'StdlogR', 'Ndata']
ObsBinn = pd.DataFrame(data=ObsBin,columns=colonnes)
Depi = np.sqrt(ObsBinn['Hypo'].values**2 - Hini**2)
res3 = WLSIC(ObsBinn, Hini, Beta, I0).do_wlsic_depth(8,13)
H = res3[0][0]
ObsFile.loc[:,'Hypo'] = ObsFile.apply(lambda row:np.sqrt(row['Depi']**2+H**2), axis=1)
ObsBin = binning(ObsFile['Iobs'].values,
                            ObsFile['Hypo'].values,
                             ObsFile['QIobs2'].values,
                             I0,0.25, 10, int(NumEvt), Ic[NumEvt], Dc[NumEvt], 30)
colonnes = ['EVID', 'Hypo', 'I', 'Io', 'QIo', 'StdI', 'StdlogR', 'Ndata']
ObsBinn = pd.DataFrame(data=ObsBin,columns=colonnes)
res4 = WLSIC(ObsBinn, H, Beta, I0).do_wlsic_I0(5.5,7.5)

dData = ObsBinn['I']-(res4[0][0]+Beta*np.log10(ObsBinn['Hypo'].values/H))
G = np.log10(ObsBinn['Hypo'].values/H)
Cd = ObsBinn['StdI'].values**2
Wd = 1./Cd
Wd = Wd/sum(Wd)
mod = sm.WLS(dData, G, w=Wd)
resB = mod.fit()
"""
#Cd = ObsBinn['StdI'].values**2
#Wd = 1./Cd
#Wd = Wd/sum(Wd)
#
#def EMIPE_K(x, Depi):
#    I = x[1] + Beta*np.log10(np.sqrt(Depi**2+ x[0]**2)/ x[0])
#    return I
#
#
#def EMIPE_K_JACdH(x, Depi, Ibin):
#    Hypo = np.sqrt(Depi**2+ x[0]**2)
#    tmpValue =  x[0]/Hypo
#    G = np.empty((len(Depi), 2))
#    G[:,0] = Beta*(tmpValue**2-1.)/( x[0]*np.log(10))
#    G[:,1] = np.ones(len(Depi))
#    return G
#
#
#def func_inv(x, Depi, Ibin):
#    return EMIPE_K(x, Depi) - Ibin
#
#Ibin = ObsBinn['I'].values
#Hinv = np.ones(len(Ibin))*Hini
#
#
#
#res = least_squares(func_inv, np.array([Hini, 6.5]), 
#                    jac= EMIPE_K_JACdH, bounds=([8,5.5], [13,7.5]),
#                    args=(Depi, Ibin), verbose=1, ftol=5e-3)
#print(res.x)

#def EMIPE_K2(Depi, H, I0):
#    I = I0 + Beta*np.log10(np.sqrt(Depi**2+H**2)/H)
#    return I
#
#
#def EMIPE_K_JACdH2(Depi, H, I0):
#    Hypo = np.sqrt(Depi**2+H**2)
#    tmpValue = H/Hypo
#    #g = Beta*(tmpValue**2-1.)/(H*np.log(10))
#    G = np.empty((len(Depi),2))
#    G[:,0] = Beta*(tmpValue**2-1.)/(H*np.log(10))
#    G[:,1] = np.ones(len(Depi))
#    return G
#
#def EMIPE_K2(Depi, H):
#    I = I0 + Beta*np.log10(np.sqrt(Depi**2+H**2)/H)
#    return I
#
#
#def EMIPE_K_JACdH2(Depi, H):
#    Hypo = np.sqrt(Depi**2+H**2)
#    tmpValue = H/Hypo
#    g = Beta*(tmpValue**2-1.)/(H*np.log(10))
#    #G = np.empty((len(Depi),2))
#    #G[:,0] = Beta*(tmpValue**2-1.)/(H*np.log(10))
#    #G[:,1] = np.ones(len(Depi))
#    return g.reshape(len(Depi),1)


#res2 = curve_fit(EMIPE_K2, Depi, Ibin, p0=(Hini, 6.5),
#                 jac= EMIPE_K_JACdH2, bounds=([8,5.5], [13,7.5]),
#                 sigma=Cd, ftol=5e-2, verbose=1)



#plt.semilogx(ObsBinn['Hypo'].values, ObsBinn['I'].values, 'o')
#plt.semilogx(ObsBinn['Hypo'].values, EMIPE_K(np.array([Hini, I0]), Depi), 'Orange')
#plt.semilogx(ObsBinn['Hypo'].values, EMIPE_K(np.array([res.x[0], res.x[1]]), Depi), 'r')
#plt.semilogx(ObsBinn['Hypo'].values, EMIPE_K(res2[0][0], Depi), '--b')