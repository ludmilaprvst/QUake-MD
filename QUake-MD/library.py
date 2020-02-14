# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:06:09 2019

@author: Baize-Funck Amelie

Binning function by intensity level
"""

import numpy as np
import time
from math import *
from numba import jit


#@jit
def RAVG_c(IObs: np.ndarray, Hypo: np.ndarray, QIobs: np.ndarray, I0: np.float,
           QI0: np.float, depth: np.float, evid: np.int, Ic: np.float, Nbin: np.int):
    """
    Binning function by intensity level
    
    The function groups intensity data by intensity levels by the geometrical mean
    of the intensity data points (IDP) of the intensity level. The mean is weighted. 
    A standard deviation of intensity is computed for each intensity bin, based
    on the weighted standard deviation of the hypocentral distance mean.
    
    :param IObs: the intensities values of the IDP
    :param Hypo: the hypocentral distances of the IDP
    :param QIobs: quality factors converted into standard deviation associated to the IDP
    :param I0: epicentral intensity values
    :param QI0: standard deviation based on quality factor of the epicentral intensity value
    :param depth: depth of the hypocenter
    :param evid: id of the earthquake
    :param Ic: intensity f completeness
    :return: a table with the intensity bins. First column: id of the earthquake,
    second column: hypocentral distance of the bin, third column: intensity of the bin,
    fourth column: epicentral intensity, fifth column: standard deviation of the epicentral intensity,
    sixth column: intensity standard deviation of the bin, seventh column: hypocentral distance
    standard deviation of the bin, eighth column: number of data in the bin
    """
    
    ClassWidth = 0.5
    MinDataPerClass = 1.0
    dI = 0.5
    compt = 0
    IntenObsMin = min_ndarray_jit(IObs)
    IntenObsMax = max_ndarray_jit(IObs)
    
    if Ic < 12:
        IbinMin = max([IntenObsMin, Ic])
    else:
        IbinMin = max([IntenObsMin, 1])
    IbinMax =  max([IntenObsMax, I0])
    Intensity_range = np.arange(IbinMax, IbinMin-dI, -dI)  
    Sortie = np.zeros((Intensity_range.shape[0], 8))

    for ii in range(Intensity_range.shape[0]):
        Intensity = Intensity_range[ii]
        IMin = Intensity - ClassWidth/2.
        IMax = Intensity + ClassWidth/2.
        index_temp = where_ndarray_jit(IObs, IMin, IMax)
        if len(index_temp) >= MinDataPerClass:
            QIobstemp = 1./QIobs[index_temp]
            Hypotemp = Hypo[index_temp]
            Hypotemp[Hypotemp ==0] = 0.1
            Iobstemp = IObs[index_temp]
            poids = QIobstemp**2
            LogR = np.log10(Hypotemp)
            Iavg = average_ndarray_jit(Iobstemp, poids)            
            Ravg = average_ndarray_jit(LogR, poids)
            StdlogR = sum(poids*((LogR-Ravg)**2))
            StdlogR = np.sqrt((StdlogR/sum(poids)))
            Ravg = 10**Ravg
            MinStdI = 0.710/np.sqrt(len(index_temp))
            StdI = abs(-3.5)*StdlogR
            StdI = np.max([StdI,MinStdI])
            if StdI == 0:
                StdI = 0.0001
            if compt ==0:
                Sortie[compt] = [evid, Ravg, Iavg, I0, QI0, StdI, StdlogR, len(index_temp)]
            elif (compt > 0):
                if (abs(np.log10(Sortie[compt-1][1]/Ravg)) >= 0.02) or (abs(Sortie[compt-1][2] - Iavg) >= 0.1):
                    Sortie[compt] = [evid, Ravg, Iavg, I0, QI0, StdI, StdlogR, len(index_temp)]
            compt += 1
    return Sortie


@jit
def RAVG_c_without(IObs: np.ndarray, Hypo: np.ndarray, QIobs: np.ndarray, I0: np.float, QI0: np.float, depth: np.float, evid: np.int, Ic: np.float, Nbin: np.int):
    
    ClassWidth = 0.5
    MinDataPerClass = 1.0
    dI = 0.5
    compt = 0
    IntenObsMin = np.min(IObs)
    IntenObsMax = np.max(IObs)
    if Ic < 12:
        iStart = np.max([round(Ic/dI)+1,round(IntenObsMin/dI)+1])
    else:
        iStart = np.max([1,round(IntenObsMin/dI)+1])
    iEnd = np.max([round(IntenObsMax/dI)+1, round(I0/dI)+1])
    
    Sortie = np.zeros((len(np.arange(iStart,iEnd+1,1)),8))
    for ii in np.arange(iStart,iEnd+1,1) [::-1]:
        Intensity = (ii-1)*dI
        IMin = Intensity - ClassWidth/2.
        IMax = Intensity + ClassWidth/2.
        index_temp = np.where((IObs>=IMin) & (IObs<=IMax))
        if len(index_temp[0]) >= MinDataPerClass:
            QIobstemp = 1./QIobs[index_temp[0]]
            Hypotemp = Hypo[index_temp[0]]
            Hypotemp[Hypotemp ==0] = 0.1
            Iobstemp = IObs[index_temp[0]]
            poids = QIobstemp**2
            LogR = np.log10(Hypotemp)
            Iavg = np.average(Iobstemp,weights=poids)
            Ravg = np.average(LogR,weights=poids)
            StdlogR = sum(poids*((LogR-Ravg)**2))

            StdlogR = np.sqrt((StdlogR/sum(poids)))
            Ravg = 10**Ravg
            MinStdI = 0.710/np.sqrt(len(index_temp[0]))
            StdI = abs(-3.5)*StdlogR
            StdI = np.max([StdI,MinStdI])
            if StdI == 0:
                StdI = 0.0001
            if (abs(np.log10(Sortie[compt-1][1]/Ravg)) >= 0.02) or (abs(Sortie[compt-1][2] - Iavg) >= 0.1):
                Sortie[compt] = [evid, Ravg, Iavg, I0, QI0, StdI, StdlogR, len(index_temp[0])]
            compt += 1
    return Sortie




@jit
def aproxi(number):
    return round(number)

@jit
def somme_jit(narray: np.ndarray):
    return sum(narray)

@jit
def floaty(a):
    return float(a)

@jit
def taille(narray):
    return len(narray)

@jit
def absolute(p):
    return abs(p)

@jit 
def maxJ(narray: np.ndarray):
    return np.max(narray)

@jit
def Distance_c(depth: int, Depi: np.ndarray):
    Hypo = []
    i = 0
    while i < Depi.size:
        Hypo.append(np.sqrt(Depi[i]**2 + depth**2))
        i += 1
    return np.array(Hypo)

@jit
def min_int_jit(a: int, b: int):
    if a < b:
        return a
    return b
    
@jit
def max_int_jit(a: int, b: int):
    if a > b:
        return a
    return b

@jit
def min_ndarray_jit(narray: np.ndarray):
    if narray.size == 0:
        return np.NaN
    i = 0
    imax = narray[0]
    while i < narray.size:
        element = narray[i]
        if imax > element:
            imax = element
        i+=1
    return imax

@jit
def max_ndarray_jit(narray: np.ndarray):
    if narray.size == 0:
        return np.NaN
    i = 0
    imax = narray[0]
    while i < narray.size:
        element = narray[i]
        if imax < element:
            imax = element
        i+=1
    return imax

@jit
def where_ndarray_jit(narray: np.ndarray, IMin, IMax):
    new_array = []
    i = 0
    while i < narray.size:
        IObs = narray[i]
        if ((IObs>=IMin) & (IObs<=IMax)):
            new_array.append(i)
        i += 1
    new = np.array(new_array)
    return new

@jit
def average_ndarray_jit(narray: np.ndarray, weight: np.ndarray):
    summe=0
    scoef=0
    i=0
    if not narray.size == weight.size:
        return np.NaN
    while i < narray.size:
        summe += narray[i] * weight[i]
        scoef += weight[i]
        i += 1
    return summe / scoef


class Timer(object):
    def start(self):
        if hasattr(self, 'interval'):
            del self.interval
        self.start_time = time.time()
        
    def stop(self):
        if hasattr(self, 'start_time'):
            self.interval = time.time() - self.start_time
            del self.start_time















    
    
    
    
    
    
    
    
    
    
    
    
    
