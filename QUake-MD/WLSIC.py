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

#try:
#    import Fonctions_bin as a
#except:
#    import Fonctions_bin_p37 as a
#from Modules_Getbeta import CalcDist, read_obsfile, read_evtfile, read_critfile


class WLSIC_Kov_oneEvt():
    """
    Set of functions that inverse depth and epicentral intensity
    from macroseismic data for a given earthquake.
    The mathematical formulation is the Koveslighety equation:
        
        I = I0 + BETA.log10(Hypo/H) + GAMMA.(Hypo-H)
    
    where I is the intensity value, I0 the epicentral intensity, BETA the
    geometric attenuation coefficient, Hypo the hypocentral distance, H the
    hypocentral depth and GAMMA the intrisic attenuation coefficient.
    The endog data are the epicentral distance and the exog data the associated
    epicentral distance.
    """
    def __init__(self, ObsBinn, depth, Beta, Gamma, I0):
        """
        :param Obsbinn: dataframe with the binned intensity data for one earthquake.
                        This dataframe should at least have have the following
                        columns : 'I', 'Depi', 'StdI', which are respectively
                        the binned intensity, the associated epicentral distance
                        and the associated standard deviation.
        :param depth: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :param Beta: the geometric attenuation coefficient of the Koveslighety equation
        :param Gamma: the intresic attenuation coefficient of the Koveslighety equation
        :param I0: the epicentral intensity of the considered earthquake
        :type Obsbinn: pandas.DataFrame
        :type depth: float
        :type Beta: float
        :type Gamma: float
        :type I0: float 
        """
        self.beta = Beta
        self.gamma = Gamma
        self.Obsbin = ObsBinn
        self.depth = depth
        self.I0 = I0
        
    def EMIPE_H(self, Depi, H):
        """
        Function used to inverse depth
        :param Depi: epicentral distances associated to the binned intensity data
        :param depth: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :type Depi: numpy.array
        :type depth: float
        """
        I = self.I0 + self.beta*np.log10(np.sqrt(Depi**2+H**2)/H) + self.gamma*(np.sqrt(Depi**2+H**2)-H)
        return I

    def EMIPE_JACdH(self, Depi, H):
        """
        The jacobian fuunction associated to EMIPE_H
        :param Depi: epicentral distances associated to the binned intensity data
        :param depth: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :type Depi: numpy.array
        :type depth: float
        """
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        g = self.beta*(tmpValue**2-1.)/(H*np.log(10)) +  self.gamma*(tmpValue-1)
        return g.reshape(len(Depi),1)
    
    def EMIPE_I0(self, Depi, I0):
        """
        Function used to inverse epicentral intensity.
        :param Depi: epicentral distances associated to the binned intensity data
        :param I0: epicentral intensity
        :type Depi: numpy.array
        :type I0: float
        """
        I = I0 + self.beta*np.log10(np.sqrt(Depi**2+self.depth**2)/self.depth)+ self.gamma*(np.sqrt(Depi**2+self.depth**2)-self.depth)
        return I

    def EMIPE_JACdI0(self, Depi, I0):
        """
        The jacobian fuunction associated to EMIPE_I0
        :param Depi: epicentral distances associated to the binned intensity data
        :param I0: epicentral intensity
        :type Depi: numpy.array
        :type I0: float
        """
        g = np.ones(len(Depi))
        return g.reshape(len(Depi), 1)
    
    def do_wlsic_depth(self, depth_inf, depth_sup):
        """
        Function used to launch depth inversion within limits.
        
        :param depth_inf: lower depth limit of inversion
        :param depth_sup: upper depth limit of inversion
        :type depth_inf: float
        :type depth_sup: float
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted depth. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)). 
        """
        Ibin = self.Obsbin['I'].values
        Depi = self.Obsbin['Depi'].values
        resH = curve_fit(self.EMIPE_H, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdH, bounds=(depth_inf, depth_sup),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH
    
    def do_wlsic_I0(self, I0_inf, I0_sup):
        """
        Function used to launch epicentral intensity inversion within limits
        :param depth_inf: lower epicentral intensity limit of inversion
        :param depth_sup: upper epicentral intensity limit of inversion
        :type I0_inf: float
        :type I0_sup: float
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted epicentral intensity. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.Obsbin['I'].values
        Depi = self.Obsbin['Depi'].values
        resI0 = curve_fit(self.EMIPE_I0, Depi, Ibin, p0=self.I0,
                                  jac= self.EMIPE_JACdI0, bounds=(I0_inf, I0_sup),
                                  sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                  xtol=1e-2, loss='soft_l1')
        return resI0
    
    
class WLS_Kov():
    """
    Set of functions that inverse the coefficients of the Koveslighety equation:
    
        I = I0 + BETA.log10(Hypo/H) + GAMMA.(Hypo-H)
    
    where I is the intensity value, I0 the epicentral intensity, BETA the
    geometric attenuation coefficient, Hypo the hypocentral distance, H the
    hypocentral depth and GAMMA the intrinsic attenuation coefficient.
    The endog data are the epicentral distance and the exog data the associated
    epicentral distance.
    """
    def __init__(self, ObsBin_plus, Beta, Gamma):
        """
        :param ObsBin_plus: dataframe with the binned intensity data for all calibration earthquakes.
                        This dataframe should at least have have the following
                        columns : 'I', 'Depi', 'Depth', 'Io',  'StdI 'and 'eqStd'
                        which are respectively the binned intensity, the associated epicentral distance,
                        the associated depth, the associated epicentral intensity,
                        the associated standard deviation and the associated
                        inverse square root of the weights used in the inversion.
        :param Beta: the geometric attenuation coefficient of the Koveslighety equation
        :param Gamma: the intresic attenuation coefficient of the Koveslighety equation
        :type ObsBin_plus: pandas.DataFrame
        :type Beta: float
        :type Gamma: float
        """
        #Variable ObsBin_plus doit contenir les Obsbin de tous les evts de la calibration, avec une colonne depth en plus
        self.ObsBin_plus = ObsBin_plus
        self.beta = Beta
        self.gamma = Gamma
        
    def EMIPE_beta(self, X, beta):
        """
        Function used to inverse the geometric attenuation coefficient
        :param X: matrix that contains epicentral distance, depth and epicentral
                  intensity associated to the binned intensity
        :param beta: geometric attenuation coefficient
        :type X: numpy.array
        :type beta: float
        """
        Depi, depths, I0s = X
        try:
            logterm = np.sqrt(Depi**2+depths**2)/depths
        except AttributeError:
            tmp = (Depi**2+depths**2)**0.5/depths
            logterm = np.array([])
            for tt in tmp:
                logterm = np.append(logterm, tt)
            
        I = I0s + beta*np.log10(logterm)
        return I
    
    def EMIPE_gamma(self, X, gamma):
        """
        Function used to inverse the intrisic attenuation coefficient
        :param X: matrix that contains epicentral distance, depth and epicentral
                  intensity associated to the binned intensity
        :param gamma: intrisic attenuation coefficient
        :type X: numpy.array
        :type gamma: float
        """
        Depi, depths, I0s = X
        try:
            logterm = np.sqrt(Depi**2+depths**2)/depths
        except AttributeError:
            tmp = (Depi**2+depths**2)**0.5/depths
            logterm = np.array([])
            for tt in tmp:
                logterm = np.append(logterm, tt)
            
        I = I0s + self.beta*np.log10(logterm) + gamma*(np.sqrt(Depi**2+depths**2)-depths)
        return I
    
    def EMIPE_beta_gamma(self, X, beta, gamma):
        """
        Function used to inverse the attenuation coefficients
        :param X: matrix that contains epicentral distance, depth and epicentral
                  intensity associated to the binned intensity
        :param beta: geometric attenuation coefficient
        :param gamma: intrisic attenuation coefficient
        :type X: numpy.array
        :type gamma: float
        :type beta: float
        """
        Depi, depths, I0s = X
        I = I0s + beta*np.log10(np.sqrt(Depi**2+depths**2)/depths) + gamma*(np.sqrt(Depi**2+depths**2)-depths)
        return I
        
    def do_wls_beta(self):
        """
        Function used to launch geometric attenuation coefficient inversion
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                 the inverted beta. pcov is a 2-D array and 
                 the estimated covariance of popt. The diagonals provide the
                 variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
                 In the case of a weighting scheme different as 'Ponderation dI',
                 please use the do_wls_beta_std() function to retrieve
                 the covariance matrix based on the intensity data standard
                 deviation just after inverting beta with the present function.
        """
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [np.array(Depi), np.array(depths), np.array(I0s)]
        resBeta = curve_fit(self.EMIPE_beta, X, Ibin, p0=self.beta,
                                  sigma=self.ObsBin_plus['eqStd'].values, absolute_sigma=True,
                                  xtol=1e-3)
        return resBeta
    
    def do_wls_beta_std(self):
        """
        Function used to compute the covariance matrix associated to the 
        geomteric attenuation coefficient based on the standard deviations associated
        to the intensity data.
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                 the inverted beta. pcov is a 2-D array and 
                 the estimated covariance of popt. The diagonals provide the
                 variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [np.array(Depi), np.array(depths), np.array(I0s)]
        resBeta = curve_fit(self.EMIPE_beta, X, Ibin, p0=self.beta, bounds=(self.beta-0.0001, self.beta+0.0001),
                                  sigma=self.ObsBin_plus['StdI'].values, absolute_sigma=True,
                                  xtol=1e-3)
        return resBeta
    
    def do_wls_gamma_std(self):
        """
        Function used to compute the covariance matrix associated to the 
        intrinsic attenuation coefficient based on the standard deviations associated
        to the intensity data.
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                 the inverted gamma. pcov is a 2-D array and 
                 the estimated covariance of popt. The diagonals provide the
                 variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [np.array(Depi), np.array(depths), np.array(I0s)]
        #
        resGamma = curve_fit(self.EMIPE_gamma, X, Ibin, p0=self.gamma, bounds=(self.gamma-1e+6, self.gamma+1e-6),
                                  sigma=self.ObsBin_plus['StdI'].values, absolute_sigma=True,
                                  xtol=1e-4)
        if resGamma[0] > 0:
            resGamma[0] = 0
        return resGamma
    
    def do_wls_gamma(self):
        """
        Function used to launch intrisic attenuation coefficient inversion
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                 the inverted gamma. pcov is a 2-D array and 
                 the estimated covariance of popt. The diagonals provide the
                 variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
                 In the case of a weighting scheme different as 'Ponderation dI',
                 please use the do_wls_gamma_std() function to retrieve
                 the covariance matrix based on the intensity data standard
                 deviation just after inverting beta with the present function.
        """
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [np.array(Depi), np.array(depths), np.array(I0s)]
        #
        resGamma = curve_fit(self.EMIPE_gamma, X, Ibin, p0=self.gamma, bounds=(-np.inf, 1e-6),
                                  sigma=self.ObsBin_plus['eqStd'].values, absolute_sigma=True,
                                  xtol=1e-4)
        if resGamma[0] > 0:
            resGamma[0] = 0
        return resGamma
    
    def do_wls_betagamma(self):
        """
        Function used to launch geometric and intrisic attenuation coefficients inversion
        
        :return: [popt, pcov] list with popt a (2, 1) shape array containing
                 the inverted beta (popt[0]) and gamma (popt[1]). pcov is a 2-D array and 
                 the estimated covariance of popt. The diagonals provide the
                 variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
                 In the case of a weighting scheme different as 'Ponderation dI',
                 please use the do_wls_betagamma_std() function to retrieve
                 the covariance matrix based on the intensity data standard
                 deviation just after inverting beta with the present function.
        """
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [Depi, depths, I0s]
        resBetaGamma = curve_fit(self.EMIPE_beta_gamma, X, Ibin, p0=[self.beta, self.gamma],
                                  bounds=([-np.inf, -np.inf], [np.inf, 0]),
                                  sigma=self.ObsBin_plus['eqStd'].values, absolute_sigma=True,
                                  xtol=1e-3)
        return resBetaGamma
    
    def do_wls_betagamma_std(self):
        """
        Function used to compute the covariance matrix associated to the 
        geomteric and intrinsic attenuation coefficient based on the standard deviations associated
        to the intensity data.
        :return: [popt, pcov] list with popt a (2, 1) shape array containing
                 the inverted beta (popt[0]) and gamma (popt[1]).
                 pcov is a 2-D array and the estimated covariance of popt.
                 The diagonals provide the variance of the parameter estimate. 
                 To compute one standard deviation errors on the parameters use
                 perr = np.sqrt(np.diag(pcov)).
        """
        
        Ibin = self.ObsBin_plus['I'].values
        Depi = self.ObsBin_plus['Depi'].values
        depths = self.ObsBin_plus['Depth'].values
        I0s = self.ObsBin_plus['Io'].values
        X = [Depi, depths, I0s]
        resBetaGamma = curve_fit(self.EMIPE_beta_gamma, X, Ibin, p0=[self.beta, self.gamma],
                                  bounds=([self.beta-0.0001, self.gamma-1e-6], [self.beta+0.0001, self.gamma+1e-6]),
                                  sigma=self.ObsBin_plus['StdI'].values, absolute_sigma=True,
                                  xtol=1e-3)
        return resBetaGamma


class WLSIC_oneEvt():
    """
    Set of functions that inverse depth and magnitude
    from macroseismic data for a given earthquake.
    The mathematical formulation is :
        
        I = C1 + C2.Mag + BETA.log10(Hypo) + GAMMA.Hypo
    
    where I is the intensity value, C1 and C2 the magnitude coefficients, M the
    magnitude, BETA the geometric attenuation coefficient, Hypo the hypocentral
    distance and GAMMA the intrisic attenuation coefficient.
    The endog data are the epicentral distance and the exog data the associated
    epicentral distance.
    """
    def __init__(self, ObsBinn, depth, mag, Beta, Gamma, C1, C2):
        """
        :param Obsbinn: dataframe with the binned intensity data for one earthquake.
                        This dataframe should at least have have the following
                        columns : 'I', 'Depi', 'StdI', which are respectively
                        the binned intensity, the associated epicentral distance
                        and the associated standard deviation.
        :param depth: hypocenter's depth of the considered earthquake. In the case
                      of depth inversion, this value is the initial depth value
        :param mag: magnitude of the considered earthquake. In the case
                      of magnitude inversion, this value is the initial magnitude value
        :param Beta: the geometric attenuation coefficient of the Koveslighety equation
        :param Gamma: the intresic attenuation coefficient of the Koveslighety equation
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :type Obsbinn: pandas.DataFrame
        :type depth: float
        :type mag: float
        :type Beta: float
        :type Gamma: float
        :type C1: float
        :type C2: float 
        """
        self.beta = Beta
        self.gamma = Gamma
        self.Obsbin = ObsBinn
        self.depth = depth
        self.C1 = C1
        self.C2 = C2
        self.mag = mag
        
    def EMIPE_H(self, Depi, H):
        """
        Function used to inverse depth
        :param Depi: epicentral distances associated to the binned intensity data
        :param H: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :type Depi: numpy.array
        :type H: float
        """
        I = self.C1 + self.C2*self.mag + self.beta*np.log10(np.sqrt(Depi**2+H**2))+self.gamma*np.sqrt(Depi**2+H**2)
        return I
    
    def EMIPE_M(self, Hypo, M):
        """
        Function used to inverse magnitude
        :param Depi: epicentral distances associated to the binned intensity data
        :param M: magnitude of the considered earthquake.
        :type Depi: numpy.array
        :type M: float
        """
        Hypo = Hypo.astype(np.float64)
        I = self.C1 + self.C2*M+ self.beta*np.log10(Hypo)+self.gamma*Hypo
        return I

    def EMIPE_JACdH(self, Depi, H):
        """
        Jacobian function used to inverse depth
        :param Depi: epicentral distances associated to the binned intensity data
        :param H: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :type Depi: numpy.array
        :type H: float
        """
        Hypo = np.sqrt(Depi**2+H**2)
        tmpValue = H/Hypo
        g = (tmpValue)*((self.beta/(np.log(10)*Hypo))+self.gamma)
        return g.reshape(len(Depi),1)

    def EMIPE_JACdM(self, Depi, H):
        """
        Jacobian function to inverse magnitude
        :param Depi: epicentral distances associated to the binned intensity data
        :param M: magnitude of the considered earthquake.
        :type Depi: numpy.array
        :type M: float
        """        
        g = self.C2*np.ones(len(Depi))
        return g.reshape(len(Depi),1)

    def do_wlsic_depth(self, depth_inf, depth_sup):
        """
        Function used to launch depth inversion within limits.
        
        :param depth_inf: lower depth limit of inversion
        :param depth_sup: upper depth limit of inversion
        :type depth_inf: float
        :type depth_sup: float
        
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted depth. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.Obsbin['I'].values
        Depi = self.Obsbin['Depi'].values
        resH = curve_fit(self.EMIPE_H, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdH, bounds=(depth_inf, depth_sup),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH
    
    def do_wlsic_depthM_std(self):
        """
        Function used to compute the covariance matrix associated to the 
        inverted depth based on the standard deviations associated
        to the intensity data.
        
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted depth. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.Obsbin['I'].values
        Depi = self.Obsbin['Depi'].values
        resH = curve_fit(self.EMIPE_H, Depi, Ibin, p0=self.depth,
                                 jac= self.EMIPE_JACdH, bounds=(self.depth-0.01, self.depth+0.01),
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 ftol=5e-2)
        return resH

    def do_wls_M(self):
        """
        Function used to launch magnitude inversion.
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted magnitude. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)).

        """
        Ibin = self.Obsbin['I'].values
        Hypo = self.Obsbin['Hypo'].values
        #np.sqrt(Depi)
        resM = curve_fit(self.EMIPE_M, Hypo, Ibin, p0=self.mag,
                                 jac= self.EMIPE_JACdM, 
                                 sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                                 xtol=1e-3)
        return resM
    
    def do_wls_M_std(self):
        """
        Function used to compute the covariance matrix associated to the 
        inverted magnitude based on the standard deviations associated
        to the intensity data.
        
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted magnitude. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)).
        """
        Ibin = self.Obsbin['I'].values
        Depi = self.Obsbin['Depi'].values
        resM = curve_fit(self.EMIPE_M, Depi, Ibin, p0=self.mag, bounds=(self.mag-0.001, self.mag+0.001),
                         jac= self.EMIPE_JACdM, 
                         sigma=self.Obsbin['StdI'].values, absolute_sigma=True,
                         ftol=1e-3)
        return resM


class WLS():
    """
    Set of functions that inverse the coefficient of the following equation
    from macroseismic data of a calibration dataset.
    The mathematical formulation is :
        
        I = C1 + C2.Mag + BETA.log10(Hypo) + GAMMA.Hypo
    
    where I is the intensity value, C1 and C2 the magnitude coefficients, M the
    magnitude, BETA the geometric attenuation coefficient, Hypo the hypocentral
    distance and GAMMA the intrisic attenuation coefficient.
    """
    def __init__(self, ObsBin_plus, C1, C2, Beta, Gamma):
        """
        :param ObsBin_plus: dataframe with the binned intensity data for all calibration earthquakes.
                        This dataframe should at least have have the following
                        columns : 'I', 'Depi', 'Depth', 'Mag', 'X', 'eqStd' and 'eqStdM'
                        which are respectively the binned intensity, the associated epicentral distance,
                        the associated depth, the associated magnitude, the X parameter used in C1/C2 inversion,
                        the associated inverse square root of the weights used in the inversion
                        of Gamma and the associated square root of the weights used in the inversion
                        of C1/C2, C1/C2/Beta, C1/C2/Beta/Gamma.
        :param Beta: the geometric attenuation coefficient
        :param Gamma: the intresic attenuation coefficient
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :type ObsBin_plus: pandas.DataFrame
        :type Beta: float
        :type Gamma: float
        :type C1: float
        :type C2: float
        """
        self.Obsbin_plus = ObsBin_plus
        self.C1 = C1
        self.C2 = C2
        self.Cb = -C1/C2
        self.Ca = 1/C2
        self.beta = Beta
        self.gamma = Gamma
        
    def EMIPE_CaCb(self, X, Ca, Cb):
        """
        Function used to calibrate C1 and C2 with magnitude as exog data
        through Ca and Cb coefficients.
        :param Ca: first magnitude equivalent coefficient. Equal to 1/Ca
        :param Cb: second magnitude equivalent coefficient. Equal to -Cb/Ca
        :param X: equivalent data, equal to the weighted mean of I - BETA.log10(Hypo) + GAMMA.Hypo
                  by earthquake. Weights are equal to the inverse of the intensity
                  standard deviation squared.
        :type X: numpy.array
        :type Ca: float
        :type Cb: float
        """
        MAG = Ca*X + Cb
        return MAG
    
    def EMIPE_C1C2BetaGamma(self, X, C1, C2, Beta, Gamma):
        """
        Function used to inverse the attenuation and magnitude coefficients
        :param X: matrix that contains magnitude, epicentral distance and depth
                  intensity associated to the binned intensity
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :param Beta: geometric attenuation coefficient
        :param Gamma: intrisic attenuation coefficient
        :type X: numpy.array
        :type Gamma: float
        :type Beta: float
        :type C1: float
        :type C2: float
        """
        mags, depi, depths = X
        hypos = np.sqrt(depi**2 + depths**2)
        I = C1 + C2*mags + Beta*np.log10(hypos) + Gamma*hypos
        return I
    
    def EMIPE_C1C2Beta(self, X, C1, C2, Beta):
        """
        Function used to inverse the magnitude coefficients and the geometric attenuation
        coefficient
        :param X: matrix that contains magnitude, epicentral distance and depth
                  intensity associated to the binned intensity
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :param Beta: geometric attenuation coefficient
        :type X: numpy.array
        :type Beta: float
        :type C1: float
        :type C2: float
        """
        mags, depi, depths = X
        hypos = np.sqrt(depi**2 + depths**2)
        I = C1 + C2*mags + Beta*np.log10(hypos)
        return I
    
    
    
    def EMIPE_gamma(self, X, Gamma):
        """
        Function used to inverse the intrisic attenuation coefficient
        :param X: matrix that contains magnitude, epicentral distance and depth
                  intensity associated to the binned intensity
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :param Gamma: intrisic attenuation coefficient
        :type X: numpy.array
        :type Gamma: float
        :type C1: float
        :type C2: float
        """
        mags, depi, depths = X
        hypos = np.sqrt(depi**2 + depths**2)
        I = self.C1 + self.C2*mags + self.beta*np.log10(hypos) + Gamma*hypos
        return I
    
    def do_wls_C1C2(self):
        """
        Function used to launch the inversion of the magnitude coefficient C1 and C2
        The DataFrame with the intensity data will need a "X" column, equal to
        the weighted mean of I - BETA.log10(Hypo) + GAMMA.Hypo by earthquake.
        Weights are equal to the inverse of the intensity standard deviation squared
        
        :return: (2, 1) shape list with in position [0] the C1 coefficient and in
                 position [1] the C2 coefficient
        """
        resCaCb = curve_fit(self.EMIPE_CaCb, self.Obsbin_plus.X.values, self.Obsbin_plus.Mag.values, 
                            p0=[self.Ca, self.Cb], sigma=self.Obsbin_plus['eqStdM'].values, absolute_sigma=True,
                            xtol=1e-3)
        C2 = 1/resCaCb[0][0]
        C1 = -resCaCb[0][1]/resCaCb[0][0]
        # a verifier pour le calcul matrice de covariance!
        #StdC2 = resCaCb[1][1]
        #StdC1 = resCaCb[1][0]/resCaCb[1][1]
        return [C1, C2]
    
    def do_wls_Gamma(self):
        """
        Function used to launch the inversion of the intrisic attenuation coefficient Gamma
        
        :return: [popt, pcov] list with popt a (1, 1) shape array containing
                the inverted gamma. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)). Values in covariance
                matrix are not accurate when using a weighting scheme different
                as 'Ponderation dI'. 
        """
        X = [self.Obsbin_plus.Mag.values,
             self.Obsbin_plus.Depi.values,
             self.Obsbin_plus.Depth.values]
        Ibin = self.Obsbin_plus['I'].values
        Gamma = curve_fit(self.EMIPE_gamma, X, Ibin, 
                                  p0=[self.gamma],
                                  sigma=self.Obsbin_plus['eqStd'].values,
                                  absolute_sigma=True,
                                  bounds=(-np.inf, 0),
                                  xtol=1e-4)
        
        return Gamma
    
    def do_wls_C1C2BetaGamma(self):
        """
        Function used to launch the inversion of all coefficients
        
        return: [popt, pcov] list with popt a (4, 1) shape array containing
                the inverted coefficient with popt[0] the C1 coefficient,
                popt[1] the C2 coefficient, popt[2] the beta coefficient and
                the popt[3] the gamma coefficient. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)). pcov values are not accurate with
                the eqStdM as sigma. A function has to be developped to compute
                accurate pcov
        """
        X = [self.Obsbin_plus.Mag.values,
             self.Obsbin_plus.Depi.values,
             self.Obsbin_plus.Depth.values]
        Ibin = self.Obsbin_plus['I'].values
        C1C2BetaGamma = curve_fit(self.EMIPE_C1C2BetaGamma, X, Ibin, 
                                  p0=[self.C1, self.C2, self.beta, self.gamma],
                                  bounds=([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, 9e-5]),
                                  sigma=self.Obsbin_plus['eqStdM'].values,
                                  absolute_sigma=True,
                                  xtol=1e-3)# 1e-6
        if C1C2BetaGamma[0][3] > 0:
            C1C2BetaGamma[0][3] = 0
        return C1C2BetaGamma
    
    def do_wls_C1C2Beta(self, sigma='none'):
        """
        Function used to launch the inversion of all coefficients, except the intrinsic
        attenuation coefficient gamma.
        
         return: [popt, pcov] list with popt a (3, 1) shape array containing
                the inverted coefficient with popt[0] the C1 coefficient,
                popt[1] the C2 coefficient and popt[2] the beta coefficient. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)). pcov values are not accurate with
                the eqStdM as sigma. A function has to be developped to compute
                accurate pcov
        """
        
            
        X = [self.Obsbin_plus.Mag.values,
             self.Obsbin_plus.Depi.values,
             self.Obsbin_plus.Depth.values]
        Ibin = self.Obsbin_plus['I'].values
        if sigma == 'none':
            sigma = self.Obsbin_plus['eqStdM'].values
        C1C2Beta = curve_fit(self.EMIPE_C1C2Beta, X, Ibin, 
                             p0=[self.C1, self.C2, self.beta],
                             bounds=([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
                             sigma=sigma,
                             absolute_sigma=True,
                             xtol=1e-3)
        return C1C2Beta

    def EMIPE_JACdC1C2BetaH(self, X, C1, C2, beta, H):
        """
        Jacobian function used to inverse depth
        :param Depi: epicentral distances associated to the binned intensity data
        :param H: hypocenter's depth of the considered earthquake.
                      Should be greater than 0
        :type Depi: numpy.array
        :type H: float
        """
        Depi, Mag = X[:2]
        aH = X[2:][0]
        Hypo = np.sqrt(Depi**2+(aH*H).sum(axis=0)**2)
        tmpValue = (aH*H).sum(axis=0)/Hypo
        gC1 = np.ones(len(Depi))
        gC2 = Mag
        gbeta = np.log10(Hypo)
        #gH = (tmpValue)*((self.beta/(np.log(10)*Hypo))+self.gamma)
        GH = np.array([])
        for ahh, hh in zip(aH, H):
            Hypo = np.sqrt(Depi**2+ahh*hh**2)
            tmpValue = ahh*hh/Hypo
            gH = (tmpValue)*((beta/(np.log(10)*Hypo)))
            gH = np.nan_to_num(gH)
            #GH = np.tile(gH, (len(Depi),1))
            try:
                GH = np.vstack((GH, gH))
            except ValueError:
                GH = np.concatenate((GH, gH))
        #g = np.array([gC1, gC2, gbeta, GH])
        g = np.vstack((gC1, gC2))
        g = np.vstack((g, gbeta))
        g = np.vstack((g, GH))
        return g.reshape(len(Depi),3+len(Depi))
    
    def EMIPE_C1C2BetaH(self, X, C1, C2, Beta, H1, H2, H3, H4, H5, H6, H7, H8,
                        H9, H10, H11, H12, H13, H14, H15, H16, H17, H18, H19, H20,
                        H21, H22, H23, H24, H25, H26, H27, H28, H29, H30, H31):
        """
        Function used to inverse the magnitude coefficients and the geometric attenuation
        coefficient
        :param X: matrix that contains magnitude, epicentral distance and depth
                  intensity associated to the binned intensity
        :param C1: first magnitude coefficient
        :param C2: second magnitude coefficient
        :param Beta: geometric attenuation coefficient
        :type X: numpy.array
        :type Beta: float
        :type C1: float
        :type C2: float
        """
        mags, depi = X[:2]
        ah = X[2:][0]
        H = np.vstack(([H1]*len(depi), [H2]*len(depi), [H3]*len(depi), [H4]*len(depi),
                       [H5]*len(depi), [H6]*len(depi), [H7]*len(depi), [H8]*len(depi),
                       [H9]*len(depi), [H10]*len(depi), [H11]*len(depi), [H12]*len(depi),
                       [H13]*len(depi), [H14]*len(depi), [H15]*len(depi), [H16]*len(depi),
                       [H17]*len(depi), [H18]*len(depi), [H19]*len(depi), [H20]*len(depi),
                       [H21]*len(depi), [H22]*len(depi), [H23]*len(depi), [H24]*len(depi),
                       [H25]*len(depi), [H26]*len(depi), [H27]*len(depi), [H28]*len(depi),
                       [H29]*len(depi), [H30]*len(depi), [H31]*len(depi)
                       ))
#        print(H.shape)
#        print(ah.shape)
        #print((H*ah).sum(axis=0))
        #ah --> (n, len(Depi)) array avec n le nombre de EQ. Chaque ligne contient
        #des 1 et des 0 et correspond a un EQ. 1 est attribue aux indices de
        # obsbin_plus.EVID ==evid concerne.
        #H --> (n, len(Depi)) array avec n le nombre de EQ, chaque ligne contient obsbin_plus.Depth
        hypos = np.sqrt(depi**2 + (H*ah).sum(axis=0)**2)
        I = C1 + C2*mags + Beta*np.log10(hypos)
        return I

    def do_wls_C1C2BetaH(self, sigma='none'):
        """
        Function used to launch the inversion of all coefficients, except the intrinsic
        attenuation coefficient gamma.
        
         return: [popt, pcov] list with popt a (3, 1) shape array containing
                the inverted coefficient with popt[0] the C1 coefficient,
                popt[1] the C2 coefficient and popt[2] the beta coefficient. pcov is a 2-D array and 
                the estimated covariance of popt. The diagonals provide the
                variance of the parameter estimate. 
                To compute one standard deviation errors on the parameters use
                perr = np.sqrt(np.diag(pcov)). pcov values are not accurate with
                the eqStdM as sigma. A function has to be developped to compute
                accurate pcov
        """
        #print(self.Obsbin_plus.columns)
        aH = np.array([])
        liste_evid = np.unique(self.Obsbin_plus.EVID.values)
        depths = np.zeros(31)
        Hmin = np.zeros(31)
        Hmax = np.zeros(31)
        for compt, evid in enumerate(liste_evid):
            ind = (self.Obsbin_plus.EVID == evid)
            depth = self.Obsbin_plus[ind]['Depth'].values[0]
            hmin = self.Obsbin_plus[ind]['Hmin'].values[0]
            hmax = self.Obsbin_plus[ind]['Hmax'].values[0]
            depths[compt] = depth
            Hmin[compt] = hmin
            Hmax[compt] = hmax
            zeros = np.zeros(len(self.Obsbin_plus.EVID))
            zeros[ind] = 1
            try:
                aH = np.vstack((aH, zeros))
            except ValueError:
                 aH = np.concatenate((aH, zeros))
    
        X = [self.Obsbin_plus.Mag.values,
             self.Obsbin_plus.Depi.values,
             aH]
        Ibin = self.Obsbin_plus['I'].values
        if sigma == 'none':
            sigma = self.Obsbin_plus['eqStdM'].values
        p0 = np.append(np.array([self.C1, self.C2, self.beta]),
                       depths)
        bounds_inf = np.append(np.array([-np.inf, -np.inf, -np.inf]),
                       Hmin)
        bounds_sup = np.append(np.array([np.inf, np.inf, np.inf]),
                       Hmax)
        C1C2BetaH = curve_fit(self.EMIPE_C1C2BetaH, X, Ibin, 
                             p0=p0,
                             bounds=(bounds_inf, bounds_sup),
                             sigma=sigma,
                             absolute_sigma=True,
                             xtol=1e-3)
        return C1C2BetaH

    def do_wls_C1C2BetaH_std(self, sigma):
        aH = np.array([])
        liste_evid = np.unique(self.Obsbin_plus.EVID.values)
        depths = np.zeros(31)
        Hmin = np.zeros(31)
        Hmax = np.zeros(31)
        for compt, evid in enumerate(liste_evid):
            ind = (self.Obsbin_plus.EVID == evid)
            depth = self.Obsbin_plus[ind]['Depth'].values[0]
            hmin = self.Obsbin_plus[ind]['Hmin'].values[0]
            hmax = self.Obsbin_plus[ind]['Hmax'].values[0]
            depths[compt] = depth
            Hmin[compt] = hmin
            Hmax[compt] = hmax
            zeros = np.zeros(len(self.Obsbin_plus.EVID))
            zeros[ind] = 1
            try:
                aH = np.vstack((aH, zeros))
            except ValueError:
                 aH = np.concatenate((aH, zeros))
    
        X = [self.Obsbin_plus.Mag.values,
             self.Obsbin_plus.Depi.values,
             aH]
        Ibin = self.Obsbin_plus['I'].values
        
        p0 = np.append(np.array([self.C1, self.C2, self.beta]),
                       depths)
        bounds_inf = np.append(np.array([-np.inf, -np.inf, -np.inf]),
                       Hmin)
        bounds_sup = np.append(np.array([np.inf, np.inf, np.inf]),
                       Hmax)
        C1C2BetaH = curve_fit(self.EMIPE_C1C2BetaH, X, Ibin, 
                             p0=p0,
                             bounds=(bounds_inf, bounds_sup),
                             sigma=sigma,
                             absolute_sigma=True,
                             xtol=1e-3)
        return C1C2BetaH[1]
    
    def compute_2Dsigma(self, eta, col='eqStdM'):
        sigma = np.diag(self.Obsbin_plus[col].values)
        sigma[sigma==0] = eta
        return sigma
    
    def do_odr_C1C2Beta(self):
        # Pour prendre en compte l'incertitude su Mw et Depi
        # Pas d'incertitude sur I????
        pass
    
    def do_odr_C1C2BetaGamma(self):
        pass