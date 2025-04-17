# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:06:50 2015

@author: Ludmila Provost
"""
import numpy as np
import pandas as pd
import library as a
try:
    import mpl_toolkits.basemap.pyproj as pyproj
except:
    import pyproj
import WLSIC

# Standard deviation of epicentral intensity based on quality factors
Std ={'A':0.5,'B':0.5,'C':0.5,'E':0.750}
# Root mean square of the inverse of the weight associated to IDP based on quality factors.
# Aim to represent equivalent standard deviations (based on weights)
# For example, the weight given to IDPs with associated quality factor A is 4. 
#Equivalent standard deviation is equal to sqrt(1/4)=0.5
Stdobs ={'A':0.5,'B':0.577,'C':0.710}

class PDF():
    def __init__(self, depth_min_ref, depth_max_ref):
        self.PdfNClassH = 50
        self.PdfNClassM = 60
        self.PdfNClassIo = 100
        
        self.PdfMinLogH = np.log10(depth_min_ref)
        self.PdfMaxLogH = np.log10(depth_max_ref)
        self.PdfMinMm = 2
        self.PdfMaxMm = 8
        self.PdfMinIo = 2
        self.PdfMaxIo = 12

        self.PdfGridSizeLogH = (self.PdfMaxLogH - self.PdfMinLogH) / float(self.PdfNClassH)
        self.EchelleLogHRef = np.logspace(self.PdfMinLogH, self.PdfMaxLogH, self.PdfNClassH)
        self.PdfGridSizeM = (self.PdfMaxMm - self.PdfMinMm) / float(self.PdfNClassM)
        self.EchelleMRef = np.linspace(self.PdfMinMm, self.PdfMaxMm, self.PdfNClassM)
        self.PdfGridSizeIo = (self.PdfMaxIo - self.PdfMinIo) / float(self.PdfNClassIo)
        self.EchelleIoRef = np.linspace(self.PdfMinIo, self.PdfMaxIo, self.PdfNClassIo)
        
        self.PDF_HM = np.zeros((self.PdfNClassH, self.PdfNClassM))
        self.PDF_HIo = np.zeros((self.PdfNClassH, self.PdfNClassIo))
        self.PDF_HMIo = np.zeros((self.PdfNClassH, self.PdfNClassM, self.PdfNClassIo))
    
    def initPDFlaw(self):
        self.PDF_HMLaw = np.zeros((self.PdfNClassH, self.PdfNClassM))
    
    def update_PDF(self, Triplets, Poids_branche):
        for Mt,Ht,It,ptriplet in zip(Triplets['Magnitude'],Triplets['Profondeur'],Triplets['Io'],Triplets['Weights']):
            indexH = round((np.log10(Ht)-self.PdfMinLogH)/self.PdfGridSizeLogH,0)        
            indexH = max([0.,indexH])
            diffH = abs(self.EchelleLogHRef-Ht)
            indexH = np.where(diffH==diffH.min())[0][0]
            
            indexM = round((Mt-self.PdfMinMm)/self.PdfGridSizeM,0)
            indexM = max([0.,indexM])
            diffM = abs(self.EchelleMRef-Mt)
            indexM = np.where(diffM==diffM.min())[0][0]

            indexIo = round((It-self.PdfMinIo)/self.PdfGridSizeIo,0)
            indexIo = max([0.,indexIo])
            diffIo = abs(self.EchelleIoRef-It)
            indexIo = np.where(diffIo==diffIo.min())[0][0]
    
            self.PDF_HM[indexH][indexM] = self.PDF_HM[indexH][indexM] + ptriplet
            self.PDF_HMLaw[indexH][indexM] = self.PDF_HMLaw[indexH][indexM] + ptriplet/Poids_branche
            self.PDF_HIo[indexH][indexIo] = self.PDF_HIo[indexH][indexIo] + ptriplet
            self.PDF_HMIo[indexH][indexM][indexIo] = self.PDF_HMIo[indexH][indexM][indexIo] + ptriplet
        
    def normalize_PDF(self):
        self.PDF_HM = self.PDF_HM/self.PDF_HM.sum()
        self.PDF_HIo = self.PDF_HIo/self.PDF_HIo.sum()
        self.PDF_HMIo = self.PDF_HMIo/self.PDF_HMIo.sum()
            



def read_empe2(Nomfichier):
    """
    Read the .txt files which contain the IPEs
    
    Formulation of the IPEs : I = C1 + C2.M + Beta.log10(Dhypo) + Gamma.Dhypo,
    with Dhypo the hypocentral distance
    
    :param Nomfichier: name and path of the IPEs file
    :type Nomfichier: str
    
    :return: 5 arrays containing the Beta coefficients, the C1 coefficients,
    the C2 coefficients, the weights and the Gamma coefficients
    of the IPEs listed in the input file
    """
    empe = np.genfromtxt(Nomfichier,skip_header=4,usecols=(0,1,2,3,4),dtype=[('Weights','f8'),('C1','f8'),('C2','f8'),('Beta','f8'),('Gamma','f8')])
    Beta = empe['Beta']
    c1 = empe['C1']
    c2 = empe['C2']
    poids = empe['Weights']
    gamma = empe['Gamma']
    try:
        Beta = np.array([float(Beta)])
        c1 = np.array([float(c1)])
        c2 = np.array([float(c2)])
        poids = np.array([float(poids)])
        gamma = np.array([float(gamma)])
    except TypeError:
        pass
    return Beta,c1,c2,poids,gamma


_GEOD = pyproj.Geod(ellps='WGS84')
def CalcDist(lon1,lat1,lon2,lat2):
    """
    Compute distance between two points from their latitude and longitude in WGS84 coordinates
    
    :param lon1: longitude of the first point
    :param lat1 : latitude of the first point
    :param lon2: longitude of the second point
    :param lat2 : latitude of the second point
    """
    try:
        _GEOD.inv(lon1,lat1,lon2,lat2)[2]/1000.
    except ValueError:
        return -1
    return _GEOD.inv(lon1,lat1,lon2,lat2)[2]/1000.

    
def Binning_Obs(ObsComplet, Param_Evt, binning, depth, I0):
    """
    Bin intensity data
    
    Prepare the inputs and call a function from library
    
    :param ObsComplet: table containing information about IDP of one earthquake,
    with the following columns: epicentral distance of the IDPs, intensity values of the IDPs,
    associated standard deviation of the IDPs
    :param Param_Evt: dictionnary containing the following information for one earthquake: intensity
    of completeness, epicentral intensity standard deviation, earthquake ID
    :param binning: a library.py function of binning
    :param depth: depth of the hypocenter of the earthquake
    :param I0: epicentral intensity
    :type ObsComplet: pandas.DataFrame
    :type Param_Evt: dict
    :type binning: function
    :type depth: float
    :type I0: float
    """
    colonnes_binn = ['EVID','Hypo','I','Io','QIo','StdI','StdlogR','Ndata']
    Depi = []
    for epi in ObsComplet['Depi'].values:
        Depi.append(epi)
    Depi = np.array(Depi)
    Hypo = a.Distance_c(depth,Depi)
    
    IOBS=[]
    for iob in ObsComplet['Iobs'].values:
        IOBS.append(iob)
    Iobs=np.array(IOBS)
    
    QIOBS=[]
    for qiob in ObsComplet['QIobs'].values:        
        QIOBS.append(Stdobs[qiob])        
    QIobs=np.array(QIOBS)
    
    Ic = Param_Evt['Ic']
    
    SortieBinn = binning(Iobs,Hypo,QIobs,I0,Param_Evt['QI0'],depth,int(Param_Evt['EVID']),Ic,30)
    ObsBinn = pd.DataFrame(data=SortieBinn,columns=colonnes_binn)
    #print Iobs,Hypo
    ObsBinn = ObsBinn[ObsBinn['EVID']!=0]
    #print ObsBinn
    return ObsBinn

def add_I02obsbin(evt, depth, StdI_0):
    """
    Add the epicentral intensity as a supplementary isoseist with an associated
    epicentral distance equal to 0.

    Parameters
    ----------
    evt : class
        Class object with the metadata associated to one earthquake,  the associated IDP
        and isoseists.
    depth : float
        depth of the earthquake.
    StdI_0 : float
        Standard deviation associated to the epicentral intensity.

    Raises
    ------
    ValueError
        Raise a value error when the epicentral intensity was not added.

    Returns
    -------
    ObsBin : pandas.DataFrame
        DataFrame containing the isoseist of one earthquake, with epicentral intensity
        added as a supplementary isoseist with an associated
        epicentral distance equal to 0.

    """
    ObsBin = evt.ObsBinn.copy()
    #print(ObsBin)
    len01 = len(ObsBin)
    ObsBin = pd.concat([ObsBin, pd.DataFrame({'EVID' : [evt.evid],
                            'Depi' : [0],
                            'Hypo': [depth],
                            'I': [evt.Io_ini],
                            'StdI': [StdI_0],
                            'Io': [evt.Io_ini],
                            'Io_std': [evt.QI0],
                            'StdLogR': [99],
                            'Ndata': [0]})],  ignore_index=True)
    len02 = len(ObsBin)
    if len02-len01!=1:
        raise ValueError('Pb avec ajout I0 obsbin')
    return ObsBin

def SearchBestStartDepth(evt, beta, c1, c2, gamma,
                         depth_min, depth_max, nbre_prof_test):
    """
    Search the start depth and magnitude for the depth and magnitude inversion
    
    This function uses the IDPs of one earthquake and one IPE to search starting values
    for the depth and magnitude inversions. Use a WRMS grid search: the couple with the minimal WRMS will be chosen as starting values
    
    :param evt: Class object with the metadata associated to one earthquake,  the associated IDP
    and isoseists.
    :param method_bin: name of the binning strategy
    :param beta: beta coefficient of the IPE
    :param c1: c1 coefficent of the IPE
    :param c2: c2 coefficient of the IPE
    :param gamma: gamma coefficient of the IPE
    :param depth_min: allowed minimal depth for the WRMS search and the depth inversion
    :param depth_max: allowed maximal depth for the WRMS search and the depth inversion
    :param nbre_prof_test: number of tested depths for the WRMS search

    
    :type evt: class
    :type methode_bin: str
    :type beta: float
    :type c1: float
    :type c2: float
    :type gamma: float
    :type depth_min: float
    :type depth_max: float
    :type nbre_prof_test: int
    
    
    :return: starting values of depth and magnitude
    """
    nObsMin = 1
    
    prof_testees = np.linspace(depth_min ,depth_max, nbre_prof_test)
    wrms_array = 1000*np.ones(len(prof_testees))
    magnitude = np.zeros(len(prof_testees))
    EvtOk = False
    for ii, depth in enumerate(prof_testees):
        StdI_0 = max([Std['A'], evt.QI0])
        #evt.Binning_Obs(depth, evt.Ic, method_bin=method_bin)
        evt.updatedepth_obsin(depth)
        ObsBin = add_I02obsbin(evt, depth, StdI_0)

        if ObsBin.shape[0]>nObsMin:

            Mag, StdM = WLSIC.WLSIC_oneEvt(ObsBin, depth, 1, beta, gamma, c1, c2).do_wls_M()

            Dhypo = ObsBin['Hypo'].values
            Iobs = ObsBin['I'].values
            Wd = 1./ObsBin['StdI'].values
            Ipred = c1 + c2*Mag + beta*np.log10(Dhypo)+gamma*Dhypo
            wrms = np.sum(Wd[:-1]*(Iobs[:-1]-Ipred[:-1])**2)/np.sum(Wd[:-1])
            wrms_array[ii] = np.sqrt(wrms)
            magnitude[ii] = Mag
            EvtOk = True
    if EvtOk:
        index_minWRMS = np.where(wrms_array == wrms_array.min())
        start_depth = prof_testees[index_minWRMS[0]]
        start_mag = magnitude[index_minWRMS[0]]
    if not EvtOk:
        print('Evt ' + str(int(evt.evid)) + ':')
        print('Not enough data for the search of Best start depth')
        index_depth10 = np.where(prof_testees == 10)
        start_depth = prof_testees[index_depth10[0]]
        start_mag = magnitude[index_depth10[0]]
        return False,False
    return start_depth[0],start_mag[0]




def weight_percentile(data, weight, percentile):
    """
    Compute weighted percentile of a distribution of data
    
    :param data: array of data
    :param weight: array of the weights associated to the data
    :param percentile: wished percentile (between 0 and 1)
    :type data: numpy.array
    :type weight: numpy.array
    :type percentile: float
    :return: the value of data corresponding to the weighted percentile
    """
    weight = weight/sum(weight)
    d, w = (list(t) for t in zip(*sorted(zip(data, weight))))
    p = np.array(w).cumsum()
    y = np.interp(percentile, p, d)
    return y

def ok4depthinv(ObsBin, mag, depth,
                C1, C2, Beta, gamma,
                depth_min, depth_max):
    Ipred = C1 + C2*mag + Beta*np.log10(ObsBin['Hypo'].values)+gamma*ObsBin['Hypo'].values
    dIwls = ObsBin['I'].values-Ipred
    Hypo = ObsBin['Hypo'].values
    g = (depth/Hypo)*(Beta/np.log(10)/Hypo+gamma)
    Wdbin = 1./(ObsBin['StdI'].values**2)
    Wdbin = Wdbin/Wdbin.sum()
    if (depth>=depth_max):
       test_inv_d = sum(g*Wdbin*dIwls)
       if -test_inv_d>0:
           test_inv = True
       else:
           test_inv = False
    elif (depth<=depth_min):
       test_inv_d = sum(g*Wdbin*dIwls) 
       if test_inv_d>0:
           test_inv = True
       else:
           test_inv = False
    else:
        test_inv= True
    return test_inv

def ok4I0inv(evt, ObsBin, mag, depth,
             C1, C2, Beta, gamma):
    Wdbin = 1./(ObsBin['StdI'].values**2)
    Wdbin = Wdbin/Wdbin.sum()
    Ipred = evt.Io_inv + Beta * np.log10(ObsBin['Hypo'].values/depth) + gamma*(ObsBin['Hypo'].values-depth)
    dIwls = ObsBin['I'].values-Ipred
    ndata = len(dIwls)
    g = np.ones(ndata)
    test_inv = False
    if (evt.Io_inv>=evt.Io_sup):
        test_inv_d = sum(g*Wdbin*dIwls)
        if -test_inv_d>0:
            test_inv = True
        else:
            test_inv = False
    elif (evt.Io_inv<=evt.Io_inf):
        test_inv_d = sum(g*Wdbin*dIwls) 
        if test_inv_d>0:
            test_inv = True
        else:
            test_inv = False
    else:
            test_inv = True
    return test_inv

def inversion_MHI0(evt,
                   C1, C2, Beta, gamma,
                   start_depth, start_mag,
                   depth_min, depth_max,
                   imposed_depth, I0_option):
    """
    Inverse magnitude and depth (and I0) based on IDPs and one IPE for one earthquake
    
    Depth, epicentral intensity and magnitude are sequentially inverted in this order. Two options
    allow not inverting depth and/or epicentral intensity. Those inversions are iterative and stop
    if the solutions converge or if more than 100 iterations are needed to converge.
    Depth and I0 are inverted within limits.
    
    :param NumEvt: ID of the earthquake
    :param DataObs: table with the IDPs and associated standard deviations
    :param methode_bin: name of the binning strategy
    :param binning: a library.py function of binning
    :param I_value: deprecated parameter, should disapear in future versions
    :param C1: C1 coefficient of the IPE
    :param C2: C2 coefficient of the IPE
    :param Beta: Beta coefficent of the IPE
    :param gamma: Gamma coefficient of the IPE
    :param start_depth: starting value of depth used in the inversion process
    :param start_mag: starting value of magnitude used in the inversion process
    :param Param_Evt: dictionnary containing the following information for one earthquake: intensity
    of completeness, epicentral intensity standard deviation, earthquake ID, epicentral intensity
    :param depth_min: minimal limit for the depth inversion
    :param depth_max: maximal limit for the depth inversion
    :param imposed_depth: imposed depth option. To invert depth, value of imposed_depth should be False, else
    the value of imposed_depth should be equal to the imposed depth.
    :param I0_option: inversion of I0 option. If equal to True, I0 is inverted. If equal to False, I0 is not inverted.
    
    :type NumEvt: float
    :type DataObs: pandas.DataFrame
    :type methode_bin: str
    :type binning: function
    :type I_value: float
    :type C1: float
    :type C2: float
    :type Beta: float
    :type gamma: float
    :type start_depth: float
    :type start_mag: float
    :type Param_Evt: dict
    :type depth_min: float
    :type depth_max: float
    :type imposed_depth: Boolean or float
    :type I0_option: Boolean
    
    :return: ivnerted magnitude, depth and I0 with the magnitude and depth standard deviations,
    updated Param_Evt, binned IDPs computed with inverted depth
    """
    # Initialization 
    depth = start_depth
    mag = start_mag
    I0 = evt.Io_ini

    iteration = 0
    MaxIter = 100.
    NbreMinIter = 3.
    StdH_fin = 0
    StdM_fin = 0
    StdI0_fin = 0
        

    evt.updatedepth_obsin(depth)
    StdI_0 = max([Std['A'], evt.QI0_inv])
    StdI_0 = np.sqrt(StdI_0/(0.1*Std['A']))
    ObsBin = add_I02obsbin(evt, depth, StdI_0)        

    while iteration<MaxIter:
        print(iteration)
        iteration+=1
        # Depth inversion
        if not imposed_depth:         
            test_inv = False
            test_inv = ok4depthinv(ObsBin, mag, depth,
                                   C1, C2, Beta, gamma,
                                   depth_min, depth_max)
            if test_inv:
                resH = WLSIC.WLSIC_oneEvt(ObsBin, depth, mag, Beta, gamma, C1, C2).do_wlsic_depth(depth_min, depth_max)
                depth = resH[0][0]
                #StdI_0 = max([Std['A'], Param_Evt['QI0_inv']])
                #resH = WLSIC.WLSIC_oneEvt(ObsBin, depth, mag, Beta, gamma, C1, C2).do_wlsic_depthM_std()
                StdH_fin = np.sqrt(np.diag(resH[1][0]))
            else:
               StdH_fin = 0.
               depth = depth
            evt.updatedepth_obsin(depth)
            ObsBin = add_I02obsbin(evt, depth, StdI_0) 

        else:
            depth = imposed_depth
            StdH_fin = 1
        # I0 inversion
        if I0_option :         
            test_inv = ok4I0inv(evt, ObsBin, mag, depth,
                                C1, C2, Beta, gamma)
            if test_inv:
                resI0 = WLSIC.WLSIC_Kov_oneEvt(ObsBin, depth, Beta, gamma, I0).do_wlsic_I0(evt.Io_inf, evt.Io_sup)
                I0 = resI0[0][0]
                StdI0_fin = np.sqrt(resI0[1][0])
                evt.QI0_inv = StdI0_fin
            else:
                I0 = I0
                evt.Io_inv = I0
                StdI0_fin = 0.5
                evt.QI0_inv = StdI0_fin
            
        else:
            I0 = evt.Io_ini
            evt.Io_inv = I0
            StdI0_fin = 0.5
            evt.QI0_inv = StdI0_fin

        # Magnitude inversion
        StdI_0 = max([Std['A'], evt.QI0_inv])
        StdI_0 = np.sqrt(StdI_0/(0.1*Std['A']))
        ObsBin.loc[ObsBin.Ndata==0, 'StdI'] = StdI_0
        ObsBin.loc[ObsBin.Ndata==0, 'I'] = evt.Io_inv
        try:    
            resM = WLSIC.WLSIC_oneEvt(ObsBin, depth, mag, Beta, gamma, C1, C2).do_wls_M()
            mag  = resM[0][0]
            StdM_fin = np.sqrt(np.diag(resM[1][0]))
        except:
                print('Singular')
                Singular=True
        # Convergence test
        if iteration == 1:
            I0conv = []
            DEPTHconv = []
            MagConv = []
        
        if len(I0conv)<NbreMinIter+1:
            I0conv.append(I0)
            DEPTHconv.append(depth)
            MagConv.append(mag)
        else:
            I0conv[0] = float(I0conv[1])
            I0conv[1] = float(I0conv[2])
            I0conv[2] = float(I0conv[3])
            I0conv[3] = I0
            MagConv[0] = float(MagConv[1])
            MagConv[1] = float(MagConv[2])
            MagConv[2] = float(MagConv[3])
            MagConv[3] = mag
            DEPTHconv[0] = float(DEPTHconv[1])
            DEPTHconv[1] = float(DEPTHconv[2])
            DEPTHconv[2] = float(DEPTHconv[3])
            DEPTHconv[3] = depth
            
        if len(I0conv) == NbreMinIter+1:                                                  
            ConvrateI0 = sum(abs(np.diff(I0conv)))
            ConvrateDepth = sum(abs(np.diff(DEPTHconv)))
            ConvrateMag = sum(abs(np.diff(MagConv)))
            ConvrateI0 = ConvrateI0/NbreMinIter    
            ConvrateDepth = ConvrateDepth/NbreMinIter
            ConvrateMag = ConvrateMag/NbreMinIter
            if (ConvrateDepth<=0.05)and(ConvrateI0<=0.01)and(ConvrateMag<=0.01):                             
                break
        return mag, depth, I0, StdM_fin, StdH_fin, evt
