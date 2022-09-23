# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:54:31 2021

@author: PROVOST-LUD
"""
import pytest
import WLSIC
import pandas as pd
import numpy as np

def test_Kov_do_wlsic_I0_gammaeq0():
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_coeff.txt')
    
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        I0 = obsbin.Io.values[0]
        resI0 = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth, Beta, 0, I0).do_wlsic_I0(I0-1, I0+1)
        assert np.round(resI0[0][0], 3) == np.round(I0, 3)
        


def test_Kov_do_wlsic_depth_gammaeq0():
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_coeff.txt')
    
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        depth_ini = 10
        obsbin.loc[:, 'Hypo'] = np.sqrt(obsbin['Depi'].values**2 + depth_ini**2)
        Hmin = 1
        Hmax = 20
        I0 = obsbin.Io.values[0]
        resH = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth_ini, Beta, 0, I0).do_wlsic_depth(Hmin, Hmax)
        assert np.round(resH[0][0], 3) == np.round(depth, 3)
        
def test_Kov_do_wlsic_I0_gamma():
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_coeff.txt')
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    Gamma = coeff_data['Gamma'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        I0 = obsbin.Io.values[0]
        I0_ini = I0-0.5
        resI0 = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth, Beta, Gamma, I0_ini).do_wlsic_I0(I0-1, I0+1)
        assert resI0[0][0] == pytest.approx(I0, 0.001)
        
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset01_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset01_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset01_coeff.txt')
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    Gamma = coeff_data['Gamma'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        I0 = obsbin.Io.values[0]
        I0_ini = I0-0.5
        resI0 = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth, Beta, Gamma, I0_ini).do_wlsic_I0(I0-1, I0+1)
        assert resI0[0][0] == pytest.approx(I0, 0.001)
        
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset02_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset02_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset02_coeff.txt')
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    Gamma = coeff_data['Gamma'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        I0 = obsbin.Io.values[0]
        I0_ini = I0-0.5
        resI0 = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth, Beta, Gamma, I0_ini).do_wlsic_I0(I0-1, I0+1)
        assert resI0[0][0] == pytest.approx(I0, 0.001)
        
def test_Kov_do_wlsic_depth_gamma():
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset04_coeff.txt')
    
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    Gamma = coeff_data['Gamma'].values[0]
    for evid in liste_evt:
        obsbin = obs_data[obs_data.EVID==evid]
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        depth_ini = 10
        obsbin.loc[:, 'Hypo'] = np.sqrt(obsbin['Depi'].values**2 + depth_ini**2)
        I0 = obsbin.Io.values[0]
        Hmin = 1
        Hmax = 20
        resH = WLSIC.WLSIC_Kov_oneEvt(obsbin, depth_ini, Beta, Gamma, I0).do_wlsic_depth(Hmin, Hmax)
        assert np.round(resH[0][0], 1) == np.round(depth, 1)

def test_Kov_do_wls_beta():        
    obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_obs.txt')
    evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_evt.txt')
    coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset03_coeff.txt')
    
    obs_data.loc[:, 'eqStd'] = obs_data.loc[:, 'StdI']
    
    liste_evt = evt_data.EVID.values
    Beta = coeff_data['Beta'].values[0]
    beta_ini = -3
    for evid in liste_evt:
        depth = evt_data[evt_data.EVID==evid].H.values[0]
        obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
            
    resBeta = WLSIC.WLS_Kov(obs_data, beta_ini, 0).do_wls_beta()
    assert resBeta[0][0] == pytest.approx(Beta, 0.01)
    
    
def test_Kov_do_wls_gamma():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:    
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
        obs_data.loc[:, 'eqStd'] = obs_data.loc[:, 'StdI']
        
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        gamma_ini = 0
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
                
        resGamma = WLSIC.WLS_Kov(obs_data, Beta, gamma_ini).do_wls_gamma()
        assert np.round(resGamma[0][0], 5) == np.round(Gamma, 5)
    
def test_Kov_do_wls_beta_gamma():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
        obs_data.loc[:, 'eqStd'] = obs_data.loc[:, 'StdI']
        
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        beta_ini = -3
        gamma_ini = 0
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
                
        resBetaGamma = WLSIC.WLS_Kov(obs_data, beta_ini, gamma_ini).do_wls_betagamma()
        assert resBetaGamma[0][0] == pytest.approx(Beta, 0.01)
        assert np.round(resBetaGamma[0][1], 5) == np.round(Gamma, 5)
        
def test_do_wlsic_depth():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
            
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        C1 = coeff_data['C1'].values[0]
        C2 = coeff_data['C2'].values[0]
        for evid in liste_evt:
            obsbin = obs_data[obs_data.EVID==evid]
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            mag = evt_data[evt_data.EVID==evid].Mag.values[0]
            depth_ini = 10
            Hmin = 1
            Hmax = 20
            resH = WLSIC.WLSIC_oneEvt(obsbin, depth_ini, mag, Beta, Gamma, C1, C2).do_wlsic_depth(Hmin, Hmax)
            assert resH[0][0] == pytest.approx(depth, 0.01)
            
def test_do_wlsC1C2():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
        obs_data.loc[:, 'eqStd'] = obs_data.loc[:, 'StdI']
        
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        C1 = coeff_data['C1'].values[0]
        C2 = coeff_data['C2'].values[0]
        c1_ini = 1
        c2_ini = 1
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            mag = evt_data[evt_data.EVID==evid].Mag.values[0]
            temp = obs_data[obs_data.EVID==evid]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
            obs_data.loc[obs_data.EVID==evid, 'Mag'] = mag
            obs_data.loc[obs_data.EVID==evid, 'eqStdM'] = 1/len(liste_evt)
            hypo_tmp = np.sqrt(temp.Depi.values**2 + depth**2)
            X_tmp = temp.I.values - Beta*np.log10(hypo_tmp) - Gamma*hypo_tmp
            obs_data.loc[obs_data.EVID==evid, 'X'] = np.average(X_tmp, weights=1/temp.StdI.values**2)
        data4inversion = obs_data.groupby('EVID').mean()
        resC1C2 = WLSIC.WLS(data4inversion, c1_ini, c2_ini, Beta, Gamma).do_wls_C1C2()
        assert resC1C2[0] ==  pytest.approx(C1, 0.001)
        assert resC1C2[1] ==  pytest.approx(C2, 0.001)

def test_do_wls_C1C2BetaGamma():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
            
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        C1 = coeff_data['C1'].values[0]
        C2 = coeff_data['C2'].values[0]
        c1_ini = 1
        c2_ini = 1
        beta_ini = -3
        gamma_ini = -0.1
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            mag = evt_data[evt_data.EVID==evid].Mag.values[0]
            temp = obs_data[obs_data.EVID==evid]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
            obs_data.loc[obs_data.EVID==evid, 'Mag'] = mag
            obs_data.loc[obs_data.EVID==evid, 'eqStd'] = temp.StdI.values
            obs_data.loc[obs_data.EVID==evid, 'eqStdM'] = temp.StdI.values
            
        resC1C2BetaGamma = WLSIC.WLS(obs_data, c1_ini, c2_ini, beta_ini, gamma_ini).do_wls_C1C2BetaGamma()
        assert resC1C2BetaGamma[0][0] == pytest.approx(C1, 0.001)
        assert resC1C2BetaGamma[0][1] == pytest.approx(C2, 0.001)
        assert resC1C2BetaGamma[0][2] == pytest.approx(Beta, 0.001)
        assert np.round(resC1C2BetaGamma[0][3], 5) == np.round(Gamma, 5)
        
def test_do_wls_Gamma():
    liste_test_id = [ '01', '02', '03', '04']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
            
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        Gamma = coeff_data['Gamma'].values[0]
        C1 = coeff_data['C1'].values[0]
        C2 = coeff_data['C2'].values[0]
        gamma_ini = -0.1
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            mag = evt_data[evt_data.EVID==evid].Mag.values[0]
            temp = obs_data[obs_data.EVID==evid]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
            obs_data.loc[obs_data.EVID==evid, 'Mag'] = mag
            obs_data.loc[obs_data.EVID==evid, 'eqStd'] = temp.StdI.values
            
        resGamma = WLSIC.WLS(obs_data, C1, C2, Beta, gamma_ini).do_wls_Gamma()

        assert np.round(resGamma[0][0], 5) == np.round(Gamma, 5)
        
def test_do_wls_C1C2Beta():
    liste_test_id = [ '03']
    for test_id in liste_test_id:
        obs_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_obs.txt')
        evt_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_evt.txt')
        coeff_data = pd.read_csv('../Testpy_dataset/pytest_dataset' + test_id + '_coeff.txt')
            
        liste_evt = evt_data.EVID.values
        Beta = coeff_data['Beta'].values[0]
        C1 = coeff_data['C1'].values[0]
        C2 = coeff_data['C2'].values[0]
        c1_ini = 1
        c2_ini = 1
        beta_ini = -3
        for evid in liste_evt:
            depth = evt_data[evt_data.EVID==evid].H.values[0]
            mag = evt_data[evt_data.EVID==evid].Mag.values[0]
            temp = obs_data[obs_data.EVID==evid]
            obs_data.loc[obs_data.EVID==evid, 'Depth'] = depth
            obs_data.loc[obs_data.EVID==evid, 'Mag'] = mag
            obs_data.loc[obs_data.EVID==evid, 'eqStd'] = temp.StdI.values
            obs_data.loc[obs_data.EVID==evid, 'eqStdM'] = temp.StdI.values
            
        resC1C2Beta = WLSIC.WLS(obs_data, c1_ini, c2_ini, beta_ini, 0).do_wls_C1C2Beta()
        assert resC1C2Beta[0][0] == pytest.approx(C1, 0.001)
        assert resC1C2Beta[0][1] == pytest.approx(C2, 0.001)
        assert resC1C2Beta[0][2] == pytest.approx(Beta, 0.001)