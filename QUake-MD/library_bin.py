# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 16:23:37 2021

@author: PROVOST-LUD
"""

import pandas as pd
import numpy as np
import library_bin as libr


Stdobs = {'A':0.5,'B':0.577,'C':0.710,'D':1.0,'I':1.5,'K':2.0}
 


def ROBS(obsdata, depth, Ic, I0, QI0):
    """
    Function that bins intensiy data points in intensity class for one earthquake.
    For each intensity class, the epicentral distance of the intensity data points
    are meaned (geometric mean).
    Width of the intensity class : 0.
    Standard deviation associated to the intensity class is equal to the log. standard
    deviation of the meaned epicentral distances multiplied by the absolute value
    of a geometric attenuation coefficient equal to -3.5.
    
    :param obsdata: dataframe with the following columns : Iobs, value of the
                    intensity data point, Depi the associated epicentral intensity,
                    QIobs the associated quality (A, B or C) and EVID the earthquake ID.
    :param depth: hypocentral depth of the considered earthquake
    :param Ic: intensity of completeness of the macroseismic field
    :param I0: epicentral intensity of the earthquake
    :param QI0: epicentral intensity quality translated in standard deviation
    :type obsdata: pandas.DataFrame
    :type depth: float
    :type I0: float
    :type QI0: float
    
    :return: a pandas.DataFrame with the binned intensity. The columns are:
             'EVID'  the earthquake ID, 'I' the value of the binned intensity,
             'Depi' the associated epicentral distance, 'Hypo' the associated
             hypocentral distance, 'StdLogR' the log standard deviation of the
             epicentral geometric mean, 'StdI' the standard deviation associated
             to the binned intensity value, 'Io' the epicentral intensity,
             'Io_std' the epicentral intensity standard deviation and 'Ndata' the
             number of intensity data point used compute the intensity bin.
    """
    #print(obsdata.EVID.values[0])
    #print(obsdata.EVID.values[0])
    #print(Ic)
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    #print(obsdata.EVID.values[0])
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    obsdata.loc[:, 'LogDepi'] = obsdata.apply(lambda row: np.log10(row['Depi']), axis=1)
    
    grouped_byIbin = obsdata.groupby('Iobs')
    wm = lambda x: np.average(x, weights=obsdata.loc[x.index, "poids"])
    varwm = lambda x: np.sum(obsdata.loc[x.index, "poids"]*(np.average(x, weights=obsdata.loc[x.index, "poids"]) - x)**2)/np.sum(obsdata.loc[x.index, "poids"])
    obsbin = grouped_byIbin.agg({'LogDepi':[("RAVG", wm), ("VarLogR", varwm)],
                                'EVID': [('Ndata', 'count')]})
    obsbin.columns = obsbin.columns.get_level_values(1)
    obsbin.reset_index(inplace=True)
    obsbin.loc[:, 'EVID'] = obsdata.EVID.values[0]
    obsbin.loc[:, 'I0'] = I0
    obsbin.loc[:, 'QI0'] = QI0
    obsbin.loc[:, 'Depi'] = obsbin.apply(lambda row: 10**row['RAVG'], axis=1)
    obsbin.loc[:, 'Hypo'] = obsbin.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    obsbin.loc[:, 'StdLogR'] = obsbin.apply(lambda row: np.sqrt(row['VarLogR']), axis=1)
    obsbin.loc[:, 'StdI_data'] = obsbin.apply(lambda row: np.abs(-3.5)*row['StdLogR'], axis=1)
    obsbin.loc[:, 'StdI_min'] = obsbin.apply(lambda row: Stdobs['C']/np.sqrt(row['Ndata']), axis=1)
    obsbin.loc[:, 'StdI'] = obsbin.apply(lambda row: np.max([row['StdI_data'], row['StdI_min']]), axis=1)
    
    obsbin = obsbin[['EVID', 'Iobs', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'I0', 'QI0', 'Ndata']]
    obsbin.columns = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    return obsbin


def RAVG(obsdata, depth, Ic, I0, QI0):
    """
    Function that bins intensiy data points in intensity class for one earthquake.
    For each intensity class, the epicentral distance of the intensity data points
    are meaned (geometric mean).
    Width of the intensity class : 1.
    Standard deviation associated to the intensity class is equal to the log. standard
    deviation of the meaned epicentral distances multiplied by the absolute value
    of a geometric attenuation coefficient equal to -3.5.
    
    :param obsdata: dataframe with the following columns : Iobs, value of the
                    intensity data point, Depi the associated epicentral intensity,
                    QIobs the associated quality (A, B or C) and EVID the earthquake ID.
    :param depth: hypocentral depth of the considered earthquake
    :param Ic: intensity of completeness of the macroseismic field
    :param I0: epicentral intensity of the earthquake
    :param QI0: epicentral intensity quality translated in standard deviation
    :type obsdata: pandas.DataFrame
    :type depth: float
    :type I0: float
    :type QI0: float
    
    :return: a pandas.DataFrame with the binned intensity. The columns are:
             'EVID'  the earthquake ID, 'I' the value of the binned intensity,
             'Depi' the associated epicentral distance, 'Hypo' the associated
             hypocentral distance, 'StdLogR' the log standard deviation of the
             epicentral geometric mean, 'StdI' the standard deviation associated
             to the binned intensity value, 'Io' the epicentral intensity,
             'Io_std' the epicentral intensity standard deviation and 'Ndata' the
             number of intensity data point used compute the intensity bin.
    """
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    colonnes_bin = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    Iclass = np.unique(obsdata.Iobs.values)
    #Iclass = Iclass[Iclass>=Ic]
    class_width = 1
    obsbin = pd.DataFrame(columns=colonnes_bin)
    EVID = obsdata.EVID.values[0]
    for compt, Ibin in enumerate(Iclass):
        if compt == 0:
            condition = np.logical_and(obsdata.Iobs>=Ibin, obsdata.Iobs<=Ibin+class_width/2)
        else:
            condition = np.logical_and(obsdata.Iobs>=Ibin-class_width/2, obsdata.Iobs<=Ibin+class_width/2)
        tmp = obsdata[condition]
        Ndata = len(tmp)
        tmp.loc[:, 'LogDepi'] = tmp.apply(lambda row: np.log10(row['Depi']), axis=1)
        IAVG =  np.average(tmp.Iobs.values, weights=tmp.poids.values)
        RAVG = np.average(tmp.LogDepi.values, weights=tmp.poids.values)
        VarLogR = np.average((tmp.LogDepi.values-RAVG)**2, weights=tmp.poids.values)
        StdLogR = np.sqrt(VarLogR)
        StdI_data = np.abs(-3.5)*StdLogR
        StdI_min = Stdobs['C']/np.sqrt(Ndata)
        StdI = max([StdI_data, StdI_min])
        Depi = 10**RAVG
        Hypo = np.sqrt(Depi**2+depth**2)
        obsbin.loc[compt, :] = [EVID, IAVG, Depi, Hypo, StdLogR, StdI, I0, QI0, Ndata]
    return obsbin

def weighted_percentile(data, percents, weights=None):
    ''' percents in units of 1%
        weights specifies the frequency (count) of data.
    '''
    if weights is None:
        return np.percentile(data, percents)
    ind = np.argsort(data)
    d=data[ind]
    w=weights[ind]
    p=1.*w.cumsum()/w.sum()*100
    y=np.interp(percents, p, d)
    return y

def RP50(obsdata, depth, Ic, I0, QI0):
    
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    #print(obsdata.EVID.values[0])
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    obsdata.loc[:, 'LogDepi'] = obsdata.apply(lambda row: np.log10(row['Depi']), axis=1)
    
    grouped_byIbin = obsdata.groupby('Iobs')
    w50 = lambda x: weighted_percentile(x.values, 50, weights=obsdata.loc[x.index, "poids"].values)
    w16 = lambda x: weighted_percentile(x.values, 16, weights=obsdata.loc[x.index, "poids"].values)
    w84 = lambda x: weighted_percentile(x.values, 84, weights=obsdata.loc[x.index, "poids"].values)
    obsbin = grouped_byIbin.agg({'LogDepi':[("RAVG", w50), ("w16", w16), ("w84", w84)],
                                'EVID': [('Ndata', 'count')]})
    obsbin.columns = obsbin.columns.get_level_values(1)
    obsbin.reset_index(inplace=True)
    obsbin.loc[:, 'StdLogR'] = (obsbin.loc[:, 'w84'] - obsbin.loc[:, 'w16'])/2
    obsbin.loc[:, 'EVID'] = obsdata.EVID.values[0]
    obsbin.loc[:, 'I0'] = I0
    obsbin.loc[:, 'QI0'] = QI0
    obsbin.loc[:, 'Depi'] = obsbin.apply(lambda row: 10**row['RAVG'], axis=1)
    obsbin.loc[:, 'Hypo'] = obsbin.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    obsbin.loc[:, 'StdI_data'] = obsbin.apply(lambda row: np.abs(-3.5)*row['StdLogR'], axis=1)
    obsbin.loc[:, 'StdI_min'] = obsbin.apply(lambda row: Stdobs['C']/np.sqrt(row['Ndata']), axis=1)
    obsbin.loc[:, 'StdI'] = obsbin.apply(lambda row: np.max([row['StdI_data'], row['StdI_min']]), axis=1)
    
    obsbin = obsbin[['EVID', 'Iobs', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'I0', 'QI0', 'Ndata']]
    obsbin.columns = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    return obsbin


def RP84(obsdata, depth, Ic, I0, QI0):
    
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    #print(obsdata.EVID.values[0])
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    obsdata.loc[:, 'LogDepi'] = obsdata.apply(lambda row: np.log10(row['Depi']), axis=1)
    
    grouped_byIbin = obsdata.groupby('Iobs')
    w84 = lambda x: weighted_percentile(x.values, 84, weights=obsdata.loc[x.index, "poids"].values)
    w50 = lambda x: weighted_percentile(x.values, 50, weights=obsdata.loc[x.index, "poids"].values)
    w98 = lambda x: weighted_percentile(x.values, 98, weights=obsdata.loc[x.index, "poids"].values)
    obsbin = grouped_byIbin.agg({'LogDepi':[("RAVG", w84), ("w50", w50), ("w98", w98)],
                                'EVID': [('Ndata', 'count')]})
    obsbin.columns = obsbin.columns.get_level_values(1)
    obsbin.reset_index(inplace=True)
    obsbin.loc[:, 'StdLogR'] = (obsbin.loc[:, 'w98'] - obsbin.loc[:, 'w50'])/2
    obsbin.loc[:, 'EVID'] = obsdata.EVID.values[0]
    obsbin.loc[:, 'I0'] = I0
    obsbin.loc[:, 'QI0'] = QI0
    obsbin.loc[:, 'Depi'] = obsbin.apply(lambda row: 10**row['RAVG'], axis=1)
    obsbin.loc[:, 'Hypo'] = obsbin.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    obsbin.loc[:, 'StdI_data'] = obsbin.apply(lambda row: np.abs(-3.5)*row['StdLogR'], axis=1)
    obsbin.loc[:, 'StdI_min'] = obsbin.apply(lambda row: Stdobs['C']/np.sqrt(row['Ndata']), axis=1)
    obsbin.loc[:, 'StdI'] = obsbin.apply(lambda row: np.max([row['StdI_data'], row['StdI_min']]), axis=1)
    
    obsbin = obsbin[['EVID', 'Iobs', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'I0', 'QI0', 'Ndata']]
    obsbin.columns = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    return obsbin

def RF50(obsdata, depth, Ic, I0, QI0):
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    #print(obsdata.EVID.values[0])
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    obsdata.loc[:, 'poids2'] = obsdata.apply(lambda row: row['poids']**2, axis=1)
    obsdata.loc[:, 'LogDepi'] = obsdata.apply(lambda row: np.log10(row['Depi']), axis=1)
    grouped_byIbin = obsdata.groupby('Iobs')
    w50 = lambda x: weighted_percentile(x.values, 50, weights=obsdata.loc[x.index, "poids"].values)
    w16 = lambda x: weighted_percentile(x.values, 16, weights=obsdata.loc[x.index, "poids"].values)
    w84 = lambda x: weighted_percentile(x.values, 84, weights=obsdata.loc[x.index, "poids"].values)
    obsbin = grouped_byIbin.agg({'LogDepi':[("RAVG", w50), ("w16", w16), ("w84", w84)],
                                'EVID': [('Ndata', 'count')], 'poids2':[('w2', 'sum')]})
    obsbin.columns = obsbin.columns.get_level_values(1)
    obsbin.reset_index(inplace=True)
    obsbin.loc[:, 'StdLogR'] = (obsbin.loc[:, 'w84'] - obsbin.loc[:, 'w16'])/2
    obsbin.loc[:, 'EVID'] = obsdata.EVID.values[0]
    obsbin.loc[:, 'I0'] = I0
    obsbin.loc[:, 'QI0'] = QI0
    obsbin.loc[:, 'Depi'] = obsbin.apply(lambda row: 10**row['RAVG'], axis=1)
    obsbin.loc[:, 'Hypo'] = obsbin.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    obsbin.loc[:, 'StdI_data'] = obsbin.apply(lambda row: np.abs(-3.5)*row['StdLogR'], axis=1)
    obsbin.loc[:, 'StdI_min'] = obsbin.apply(lambda row: Stdobs['C']/np.sqrt(row['Ndata']), axis=1)
    obsbin.loc[:, 'StdI'] = obsbin.apply(lambda row: np.max([row['StdI_data'], row['StdI_min']]), axis=1)
    obsbin.loc[:, 'W'] = np.sqrt(obsbin.loc[:, 'Ndata'])*obsbin.loc[:, 'w2']*obsbin.loc[:, 'RAVG']
    obsbin = obsbin[obsbin.W==obsbin.W.max()]
    obsbin = obsbin[['EVID', 'Iobs', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'I0', 'QI0', 'Ndata']]
    obsbin.columns = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    obsbin.loc[0, :] = [obsdata.EVID.values[0], I0, 0, depth, -1, QI0, I0, QI0, 0]
    return obsbin

def RF84(obsdata, depth, Ic, I0, QI0):
    obsdata = obsdata[obsdata.Iobs>=Ic]
    obsdata['Depi'].replace(0, 0.5, inplace=True)
    #print(obsdata.EVID.values[0])
    obsdata.loc[:, 'poids'] = obsdata.apply(lambda row: 1/Stdobs[row['QIobs']]**2, axis=1)
    obsdata.loc[:, 'poids2'] = obsdata.apply(lambda row: row['poids']**2, axis=1)
    obsdata.loc[:, 'LogDepi'] = obsdata.apply(lambda row: np.log10(row['Depi']), axis=1)
    grouped_byIbin = obsdata.groupby('Iobs')
    w84 = lambda x: weighted_percentile(x.values, 84, weights=obsdata.loc[x.index, "poids"].values)
    w50 = lambda x: weighted_percentile(x.values, 50, weights=obsdata.loc[x.index, "poids"].values)
    w98 = lambda x: weighted_percentile(x.values, 98, weights=obsdata.loc[x.index, "poids"].values)
    obsbin = grouped_byIbin.agg({'LogDepi':[("RAVG", w84), ("w50", w50), ("w98", w98)],
                                'EVID': [('Ndata', 'count')], 'poids2':[('w2', 'sum')]})
    obsbin.columns = obsbin.columns.get_level_values(1)
    obsbin.reset_index(inplace=True)
    obsbin.loc[:, 'StdLogR'] = (obsbin.loc[:, 'w98'] - obsbin.loc[:, 'w50'])/2
    obsbin.loc[:, 'EVID'] = obsdata.EVID.values[0]
    obsbin.loc[:, 'I0'] = I0
    obsbin.loc[:, 'QI0'] = QI0
    obsbin.loc[:, 'Depi'] = obsbin.apply(lambda row: 10**row['RAVG'], axis=1)
    obsbin.loc[:, 'Hypo'] = obsbin.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    obsbin.loc[:, 'StdI_data'] = obsbin.apply(lambda row: np.abs(-3.5)*row['StdLogR'], axis=1)
    obsbin.loc[:, 'StdI_min'] = obsbin.apply(lambda row: Stdobs['C']/np.sqrt(row['Ndata']), axis=1)
    obsbin.loc[:, 'StdI'] = obsbin.apply(lambda row: np.max([row['StdI_data'], row['StdI_min']]), axis=1)
    obsbin.loc[:, 'W'] = np.sqrt(obsbin.loc[:, 'Ndata'])*obsbin.loc[:, 'w2']*obsbin.loc[:, 'RAVG']
    obsbin = obsbin[obsbin.W==obsbin.W.max()]
    obsbin = obsbin[['EVID', 'Iobs', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'I0', 'QI0', 'Ndata']]
    obsbin.columns = ['EVID', 'I', 'Depi', 'Hypo', 'StdLogR', 'StdI', 'Io', 'Io_std', 'Ndata']
    obsbin.loc[0, :] = [obsdata.EVID.values[0], I0, 0, depth, -1, QI0, I0, QI0, 0]
    return obsbin
    
#def Binning_Obs(obsdata, depth, Ic, I0, QI0, evid):
#    # Old version for binning RAVG
#        Stdobs = {'A':0.5,'B':0.577,'C':0.710,'D':1.0,'I':1.5,'K':2.0}
#        colonnes_binn = ['EVID','Hypo','I','Io','QIo','StdI','StdLogR','Ndata']
#        Depi = []
#        for epi in obsdata['Depi'].values:
#            Depi.append(epi)
#        Depi = np.array(Depi)
#        
#        Hypo = libr.Distance_c(depth, Depi)
#        IOBS = []
#        for iob in obsdata['Iobs'].values:
#            IOBS.append(iob)
#        Iobs = np.array(IOBS)
#
#        QIOBS = []
#        for qiob in obsdata['QIobs'].values:
#            QIOBS.append(Stdobs[qiob])
#        QIobs=np.array(QIOBS)
#
#        depth = float(depth)
#        evid = int(evid)
#        Ic = float(Ic)
#        SortieBinn = libr.RAVG_c(Iobs,Hypo,QIobs,I0,QI0,depth,evid, Ic,30)
#        ObsBinn = pd.DataFrame(data=SortieBinn, columns=colonnes_binn)
#        ObsBinn = ObsBinn[ObsBinn['EVID'] != 0]
#        return ObsBinn

#obsdata = pd.read_csv('../Testpy_dataset/obs_test.txt', sep='\t')
#I0 = 8
#QI0 = 0.5
#Ic = 6
#depth = 0

# test vitesse RAVG : 0.026 s Binning_Obs : 0.012 s