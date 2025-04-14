# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:44:20 2019

@author: baize-funck-ame
"""


import numpy as np
import pandas as pd
import library_bin as libr
from tkinter import messagebox as tkm
try:
    from mpl_toolkits.basemap import pyproj
except KeyError:
    import os
    os.environ['PROJ_LIB'] = '/opt/Anaconda3/pkgs/proj4-5.2.0-h470a237_1/share/proj'
    from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogLocator
import copy

font0 = FontProperties()

class PlotEvt():
    def __init__(self, Ffp):
        self.FfP = Ffp

    def build(self, evid):
        # Standard deviation of epicentral intensity based on quality factors
        Std ={'A':0.5,'B':0.5,'C':0.5,'E':0.750}
        
        self.evid = evid
        try:
            EvtFile = self.FfP.EvtFile
            ObsFile = self.FfP.ObsFile
        except AttributeError:
            tkm.showerror("Error", "Unexpected error, please reload the input data files by pushing Reset")
        if hasattr(self.FfP,'Parameter file'):
            ParameterFile = self.FfP.ParameterFile
        
        # if self.evid < 0 or not isinstance(self.evid,int):
        #     tkm.showerror("Error", "Please enter an integer")
        #     return
        
        if np.size(EvtFile[EvtFile['EVID']==self.evid].values) == 0 :
            print("L'evid " + str(self.evid) + " n'existe pas")
            tkm.showerror("Error","Evid " + str(self.evid) + " doesn't exist")
            return
        
        self.day = EvtFile[EvtFile['EVID']==self.evid]['Day'].values[0]
        self.month = EvtFile[EvtFile['EVID']==self.evid]['Month'].values[0]
        self.year = EvtFile[EvtFile['EVID']==self.evid]['Year'].values[0]
        
        self.QI0name = EvtFile[EvtFile['EVID']==self.evid]['QI0'].values[0]
        if not isinstance(self.QI0name,str):
            tkm.showerror("Error","invalid QI0")
        if not ((self.QI0name >= 'A' and self.QI0name <= 'C') or self.QI0name == 'K' or self.QI0name == 'E' or self.QI0name == 'I'):
            tkm.showerror("","invalid QI0")
            return;
        self.QI0 = Std[(self.QI0name)]
        print("QI0name:")
        print(self.QI0name)
        self.QPos = EvtFile[EvtFile['EVID']==self.evid]['QPos'].values[0]
        if isinstance(self.QPos, bytes):
            self.QPos = self.QPos.decode('utf8')
        if not isinstance(self.QPos, str):
            tkm.showerror("Error","invalid QPos")
            print(self.QPos)
            return;
        if not ((self.QPos >= 'A' and self.QPos <= 'D') or self.QPos == 'I' or self.QPos == 'K' or self.QPos == 'E'):
            tkm.showerror("Error","invalid QPos")
            print(self.QPos)
            return;
        print("QPos")
        print(self.QPos)
        self.Io_ini = EvtFile[EvtFile['EVID']==self.evid]['I0'].values[0]
        if self.Io_ini < 0 or self.Io_ini > 12:
            tkm.showerror("Error","invalid I0")
            return;
        print("I0 catalogue")
        print(self.Io_ini )
        #voir pur maj Iinf Isup LImitforsampling
        self.Io_inf = self.Io_ini - 2*self.QI0
        self.Io_sup = self.Io_ini + 2*self.QI0
        self.I0 = EvtFile[EvtFile['EVID']==self.evid]['I0'].values[0]
        """
        # Under development
        if hasattr(self.FfP, 'Parameterfile'):
            print("Use Parameterfile")
            self.Ic = ParameterFile[ParameterFile['EVID']==self.evid]['Ic'].values[0]
            if not (isinstance(self.Ic, float)):
                tkm.showerror("Error","Invalid Ic")
                return
            self.Nobs = ParameterFile[ParameterFile['EVID']==self.evid]['NObs'].values[0]
            self.Nfelt = ParameterFile[ParameterFile['EVID']==self.evid]['NFelt'].values[0]
        else:
            self.Nfelt = float(np.size(ObsFile[(ObsFile['EVID']==self.evid) & (ObsFile['Iobs'] == -1)].values)/6)
        """
        
        self.Nobs = float(np.size(ObsFile[(ObsFile['EVID']==self.evid) & (ObsFile['Iobs'] > 0)].values)/6)
        
        if not (isinstance(self.Nobs, float)) :
            tkm.showerror("Error", "Invalid Nobs ")
            return
        self.Lat_evt = float(EvtFile[EvtFile['EVID']==self.evid]['Lat'].values[0])
        self.Lon_evt = float(EvtFile[EvtFile['EVID']==self.evid]['Lon'].values[0])
        if not (isinstance(self.Lat_evt, float) and isinstance(self.Lon_evt,float)):
            tkm.showerror("Error", "invalid lat or lon")
            
        if not 'Depi' in ObsFile.columns:
            ObsFile.loc[ObsFile['EVID']==self.evid, 'Depi'] = ObsFile.loc[ObsFile['EVID']==self.evid].apply(lambda row:CalcDist(row['Lon'],row['Lat'],self.Lon_evt,self.Lat_evt),axis=1)
        ObsFile.loc[ObsFile['EVID']==self.evid,'I0'] = self.Io_ini
        
        date = EvtFile[EvtFile['EVID']==self.evid]['Year'].values[0]
        ObsFile.loc[ObsFile['EVID']==self.evid,'Year'] = date
        
        self.Obsevid = ObsFile[ObsFile['EVID']==self.evid]
    
    def Binning_Obs(self, depth, Ic, method_bin='ROBS'):
        #print(self.Obsevid.head())
        if method_bin == 'RAVG':
            self.ObsBinn = libr.RAVG(self.Obsevid, depth, Ic, self.I0, self.QI0)
        elif method_bin == 'ROBS':
            self.ObsBinn = libr.ROBS(self.Obsevid, depth, Ic, self.I0, self.QI0)
        elif method_bin == 'RP50':
            self.ObsBinn = libr.RP50(self.Obsevid, depth, Ic, self.I0, self.QI0)
        elif method_bin == 'RP84':
            self.ObsBinn = libr.RP84(self.Obsevid, depth, Ic, self.I0, self.QI0)
        elif method_bin == 'RF50':
            self.ObsBinn = libr.RF50(self.Obsevid, depth, Ic, self.I0, self.QI0)
        elif method_bin == 'RF84':
            self.ObsBinn = libr.RF84(self.Obsevid, depth, Ic, self.I0, self.QI0)
            
    def updatedepth_obsin(self, depth):      
        self.ObsBinn.loc[:, 'Hypo'] = self.ObsBinn.apply(lambda row: np.sqrt(row['Depi']**2+depth**2), axis=1)
    
    def add_invMHI0_variables(self, QI0_inv, Io_inv, Ic):
        self.QI0_inv = QI0_inv
        self.Io_inv = Io_inv
        self.Ic = Ic
    
    def deepCopy(self):
        evtcopy = PlotEvt(self.FfP)
        print(u"Recréation d'evt success")
        evtcopy.build(self.evid)
        print("build success")
        return evtcopy
#        return copy.deepcopy(self)
        
    def getFig(self):
        try:
            return self.fig
        except:
            try:
                self.PlotEvt('RAVG', 0,  option_coupure=False, option_felt=False) 
            except:
                tkm.showerror("Error","Problem with creation of the map")
            return self.fig
    
    def PlotEvt(self, methode_bin, depth, option_coupure=False, option_felt=False):
        
        if not type(option_felt) is bool:
            self.Obsevid.loc[self.Obsevid['Iobs']==-1, 'Iobs'] = option_felt
          
        # print('COUCOU 2, NOBS')
        # print(self.Nobs)
        self.fig = plt.figure(figsize=[12, 8]) 
        gs = GridSpec(1, 2, width_ratios=[1.5,1])
        axe0 = plt.subplot(gs[0,0])
        axe1 = plt.subplot(gs[0,1])
        
        jahr = self.Obsevid['Year'].values[0]
        
        latmin = self.Obsevid['Lat'].min(axis=0)
        latmax = self.Obsevid['Lat'].max(axis=0)
        lonmin = self.Obsevid['Lon'].min(axis=0)
        lonmax = self.Obsevid['Lon'].max(axis=0)
        
        latevt = self.Lat_evt
        lonevt = self.Lon_evt

        latmin = min([latmin, latevt])
        latmax = max([latmax, latevt])
        lonmin = min([lonmin, lonevt])
        lonmax = max([lonmax, lonevt])
        
        difflat = (latmax + 0.5) - (latmin - 0.5)
        difflon = (lonmax + 0.5) - (lonmin - 0.5)

        
        carte = Basemap(resolution='l',projection='merc',
                                llcrnrlat=latmin-0.5,urcrnrlat=latmax+0.5,
                                llcrnrlon=lonmin-0.5,urcrnrlon=lonmax+0.5,
                                lat_0=latevt,lon_0=lonevt,ax=axe0)
        pas = round(difflon/3)
        if pas <= 0:
            pas=0.5
        carte.drawmeridians(np.arange(round(lonmin-0.5),round(lonmax+0.5),pas), labels=[False,True,True,False])
        pas = round(difflat/3)
        if pas <= 0:
            pas=0.5
        carte.drawparallels(np.arange(round(latmin-0.5),round(latmax+0.5),pas), labels=[True,False,False,True])
        carte.drawcountries(linewidth=0.8, zorder=3)
        carte.drawcoastlines(linewidth=0.8)
        
        longscale = lonmin-0.5 + (lonmax - lonmin+1)/3
        latscale = latmin-0.5 + (latmax - latmax + 1)/5
        larg = CalcDist(lonmin-0.5,latmin, lonmax+0.5,latmin)
        if larg >= 250:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 200, barstyle='fancy', zorder=10)
        elif larg >= 200:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 100, barstyle='fancy', zorder=10)
        elif larg >= 150:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 75, barstyle='fancy', zorder=10)
        elif larg >= 100:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 50, barstyle='fancy', zorder=10)
        elif larg >= 50:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 25, barstyle='fancy', zorder=10)
        elif larg >= 20:
            carte.drawmapscale(longscale,latscale, lonmin, latmin, 10, barstyle='fancy', zorder=10)
        else:
            carte.drawmapscale(longscale, latscale, lonmin, latmin, 3, barstyle='fancy', zorder=10)
        
        carte.shadedrelief(zorder=1)

   
        Lonevt,Latevt = carte(lonevt, latevt)
#        print(lonevt,latevt)
        self.Obsevid = self.Obsevid.sort_values('Iobs')
        
        for index,row in self.Obsevid.iterrows():
            Lon,Lat = carte(row['Lon'],row['Lat'])

             #Italian style
            if (row['Iobs']==11):
                couleur = 'k'  
            elif (row['Iobs']>=9)and(row['Iobs']<=10.5):
                couleur = 'DarkViolet' 
            elif (row['Iobs']>=9)and(row['Iobs']<=10):
                couleur = 'Crimson'
            elif (row['Iobs']>=8)and(row['Iobs']<=8.5):
                couleur = 'Red'
            elif (row['Iobs']>=7)and(row['Iobs']<=7.5):
                couleur = 'Orange'
            elif (row['Iobs']>=6)and(row['Iobs']<=6.5):
                couleur = 'Yellow'
            elif (row['Iobs']>=5)and(row['Iobs']<=5.5):
                couleur = 'Lime'
            elif (row['Iobs']>=4)and(row['Iobs']<=4.5):
                couleur = 'Green'
            elif (row['Iobs']>=3)and(row['Iobs']<=3.5):
                couleur = 'LightSkyBlue'    
            elif (row['Iobs']>=2)and(row['Iobs']<=2.5):
                couleur = 'DeepSkyBlue'
            elif (row['Iobs']>=1)and(row['Iobs']<=1.5):
                couleur = 'Snow'
            else:
                continue

            carte.scatter(Lon, Lat,
                       facecolors=couleur,
                       marker='o',
                       s=30, edgecolors='DimGray', zorder=5)

        carte.plot(Lonevt,Latevt,color='k',marker='$\star$',
                   markersize=15, markerfacecolor='none', markeredgewidth=1.5, zorder=6)

        font1 = font0.copy()
        font1.set_weight('heavy')
       
        #Legende de la carte / italian style:
        legI11 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'k')
        legI10 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'DarkViolet')
        legI9 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Crimson')
        legI8 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Red')
        legI7 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Orange')
        legI6 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Yellow')
        legI5 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Lime')
        legI4 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Green')
        legI3 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'LightSkyBlue')
        legI2 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'DeepSkyBlue')
        legI1 = mlines.Line2D([],[],marker='o',markeredgecolor='DimGray',linestyle='',color= 'Snow')
        legIo = mlines.Line2D([],[],marker='$\star$',linestyle='',color= 'k',
                              markersize=10, markerfacecolor='none', markeredgewidth=1)

        # Italian style
        self.leg0 = axe0.legend((legIo,legI1,legI2,legI3,legI4,legI5,legI6,legI7,legI8,legI9,legI10, legI11),
                    (u"Epicenter",u" I / I-II",u"II / II-III",u"III / III-IV",u"IV / IV-V",u"V / V-VI",
                     u"VI / VI-VII",u"VII / VII-VIII",u"VIII / VIII-IX",u"IX / IX-X",u"X / X-XI",u"XI"),
                     numpoints=1,scatterpoints = 1, bbox_to_anchor=(0, -0.5, 1, 0.5),
                     mode='expand', ncol=2, prop={'size':10}, borderpad=1.5)
        
    
        Hypo = np.sqrt(self.Obsevid['Depi']**2+depth**2)
        
#        print('Nbre felt et distance epi moyenne:')
#        print(len(self.Obsevid[self.Obsevid['Iobs']==-1]))
#        print(self.Obsevid[self.Obsevid['Iobs']==-1]['Depi'].mean())
        axe1.semilogx(Hypo,self.Obsevid['Iobs'].values,'o', alpha=1,
                              markeredgecolor='#1f77b4', markerfacecolor='None', ms=4, label='IDP')
                
        self.Binning_Obs(0, 0)
        axe1.errorbar(self.ObsBinn['Hypo'].values, self.ObsBinn['I'].values, fmt='d',
                                yerr=self.ObsBinn['StdI'].values, mec='k', markersize=8, label='RAVG isoseismal')

        
        maxhypo = max(Hypo)
        
        axe1.set_xlabel('Epicentral distance [km]', size=12)
        axe1.set_ylabel('Macroseismic intensity', size=12)
       
        if maxhypo < 100:
            maxhypo = 100
        axe1.set_xlim([0.1, maxhypo+100])
        start, stop = axe1.get_xlim()

        axe1.xaxis.set_major_formatter(ScalarFormatter())
        locmaj = LogLocator(base=10.0, subs=(0.1, 0.2, 0.5, 1, 2, 5, 10 ))
        axe1.xaxis.set_major_locator(locmaj)
        axe1.xaxis.set_minor_formatter(NullFormatter())
        axe1.set_ylim([1, 12])
        axe1.legend(loc=3)
        
        self.fig.subplots_adjust(left=None, bottom=0.25, right=None, top=None, wspace=None, hspace=None)
        axe1.set_anchor('N')
        axe0.set_anchor('N')
        SuperTitre = 'Num evt: ' + str(int(self.evid)) + ', Year: ' + str(int(jahr)) + ', Io=' + str(self.Io_ini) + ', QI0=' + str(self.QI0name) + ', QPos=' + str(self.QPos) + ', Nobs=' + str(self.Nobs) 
        self.fig.suptitle(SuperTitre,fontsize=12)
        
        
    def Binning_Obs_old(self, depth, Ic):
        Stdobs = {'A':0.5,'B':0.577,'C':0.710,'D':1.0,'I':1.5,'K':2.0}
        colonnes_binn = ['EVID','Hypo','I','Io','QIo','StdI','StdLogR','Ndata']
        Depi = []
        for epi in self.Obsevid['Depi'].values:
            Depi.append(epi)
        Depi = np.array(Depi)
        Hypo = libr.Distance_c(depth, Depi)
#        print("Hypo:")
#        print(Hypo)
        IOBS = []
        for iob in self.Obsevid['Iobs'].values:
            IOBS.append(iob)
        Iobs = np.array(IOBS)
#        print("Iobs:")
#        print(Iobs)
        QIOBS = []
        for qiob in self.Obsevid['QIobs'].values:
            QIOBS.append(Stdobs[qiob])
        QIobs=np.array(QIOBS)
#        print("QIobs:")
#        print(QIobs)
        
        if hasattr(self,'Ic') and Ic == -1:
            Ic = self.Ic
        
#        print("Ic:")
#        print(Ic)
        
        I0 = float(self.Io_ini)
        QI0 = float(self.QI0)
        depth = float(depth)
        evid = int(self.evid)
        Ic = float(Ic)
        SortieBinn = libr.RAVG_c(Iobs,Hypo,QIobs,I0,QI0,depth,evid, Ic,30)
        self.ObsBinn = pd.DataFrame(data=SortieBinn, columns=colonnes_binn)
        self.ObsBinn = self.ObsBinn[self.ObsBinn['EVID'] != 0]
        #print(self.ObsBinn)


class FilesForPlot():
    def __init__(self, Evtname, Obsname, Parametername=""):
        try:
            self.EvtFile = pd.read_csv(Evtname, sep=';')
#            self.EvtFile = np.genfromtxt(Evtname, skip_header=1, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8),
#                                                          dtype=[('EVID', 'f8'), ('I0', 'f8'), ('QI0', 'U8'), ('Lon', 'f8'), 
#                                                                       ('Lat', 'f8'), ('QPos', 'S8'), ('Day', 'f8'), ('Month', 'f8'),
#                                                                       ('Year', 'f8')])
        except:
            tkm.showerror('Error',"Problem with Event File")
        #self.EvtFile = pd.DataFrame(self.EvtFile)
         
        try:
            self.ObsFile = pd.read_csv(Obsname, sep=';')
        except:
            tkm.showerror("Error","Problem with Observation File")
        
        if not Parametername == "":
            try:
                self.CritiqueFile = pd.read_csv(Parametername, delimiter='\t')
            except:
                tkm.showerror("Error","Problem with Ic/Depth File")
        
        
def searchByDate(EvtName, day, month, year):
#    EvtFile = np.genfromtxt(EvtName,skip_header=1,usecols=(0,1,2,3,4,5,6,7,8),dtype=[('EVID','f8'),('I0','f8'),('QI0','U8'),('Lon','f8'),('Lat','f8'),('QPos','S8'),('Day','f8'),('Month','f8'),('Year','f8')])
#    EvtFile = pd.DataFrame(EvtFile)
    EvtFile = pd.read_csv(EvtName, sep=';')
        
    # tests date correcte 
    if day < 0:
        print("Invalid date")
        tkm.showerror("Error", "Invalid date")
        return
    if month <= 0 or month > 12:
        if not month == 0:
            tkm.showerror("Error", "Invalid date")
            print("Invalid date")
            return
    elif month == 2 :
        if year % 4 != 0:
            if day > 28:
                tkm.showerror("Error", "Invalid date")
                print("Invalid date")
                return
        elif year % 100 !=0:
            if day > 29:
                tkm.showerror("Error", "Invalid date")
                print("Invalid date")
                return
        elif year % 400 != 0:
            if day > 28:
                tkm.showerror("Error", "Invalid date")
                print("Invalid date")
                return
        elif day > 29:
            tkm.showerror("Error", "Invalid date")
            print("Invalid date")
            return
    elif month <= 7:
        if month % 2 == 1:
            if day > 31:
                tkm.showerror("Error", "Invalid date")
                print("Invalid date")
                return
        if day > 30:
            tkm.showerror("Error", "Invalid date")
            print("Invalid date")
            return
    elif month % 2 == 1:
        if day > 30:
            tkm.showerror("Error", "Invalid date")
            print("Invalid date")
            return
    elif day > 31:
        tkm.showerror("Error", "Invalid date")
        print("Invalid date")
        return
                
    # filtres pour si on ne connaît que l'année, ou que le mois et l'année avec affichage du nom
    if month <= 0:
        setEvid = EvtFile[EvtFile['Year']==year]['EVID'].values
    elif day <= 0:
        setEvid = EvtFile[(EvtFile['Year']==year) & (EvtFile['Month'] == month)]['EVID'].values
    else:
        setEvid = EvtFile[(EvtFile['Year']==year) & (EvtFile['Month'] == month) & (EvtFile['Day']== day)]['EVID'].values
    #print(setEvid)
    return setEvid


def searchByLoc(EvtName, lon, lat, rad):
    EvtFile = pd.read_csv(EvtName, sep=';')
#    EvtFile = np.genfromtxt(EvtName,skip_header=1,usecols=(0,1,2,3,4,5,6,7,8),dtype=[('EVID','f8'),('I0','f8'),('QI0','U8'),('Lon','f8'),('Lat','f8'),('QPos','S8'),('Day','f8'),('Month','f8'),('Year','f8')])
#    EvtFile = pd.DataFrame(EvtFile)
        
    #tests sur la latitude et la longitude
    print(CalcDist(lon,lat,EvtFile['Lon'],EvtFile['Lat']))
    
    setEvid = EvtFile[CalcDist(lon,lat,EvtFile['Lon'],EvtFile['Lat']) <= rad]
    print(setEvid)
    return setEvid




    
######## Fonctions (shoot et equi) pour dessiner un cercle sur une map, viennent de l'ORB
def shoot(lon, lat, azimuth, maxdist=None):        
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
     
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        print("Only N-S courses are meaningful, starting at a pole!")
     
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
     
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
    
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
             d / 4. - cz) * sy * d + tu
     
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
    
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
     
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
     
    return (glon2, glat2, baz)

    
def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
 
    #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    plt.plot(X,Y, **kwargs)
    m.plot(X,Y, **kwargs)
    
_GEOD = pyproj.Geod(ellps='WGS84')
def CalcDist(lon1,lat1,lon2,lat2):
    return _GEOD.inv(lon1, lat1, lon2, lat2)[2]/1000.
