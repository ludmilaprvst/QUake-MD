#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:19:58 2019

@author: baize-funck-ame
"""

import tkinter as tk
from tkinter import messagebox as tkm
import tkinter.filedialog as fd
import PlotEvtObject as peo
import QUakeMD as qumd
import queue as que
import threading as th

class QUakeMdInterface(tk.Toplevel):
    """
    Manage the user interface of QUake-MD
    """
    def __init__(self,parent):
        tk.Toplevel.__init__(self,parent)
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        self.createView()
    
    def createView(self):
        self.canvas = tk.Canvas(self, width=1100, height=750, scrollregion=(0,0,1200,750), bg='white')
        # Create the widgets
        self.labelFileEvt = tk.Label(self, text="Evt file", bg='white')
        self.labelFileObs = tk.Label(self, text="Obs file", bg='white')
        #self.labelFileParameter = tk.Label(self, text="Ic/Dc file", bg='white')
        
        self.variableEvt = tk.StringVar()
        self.variableObs = tk.StringVar()
        #self.variableParameter  = tk.StringVar()
        self.entryEvt = tk.Entry(self, textvariable=self.variableEvt, width=30)
        self.entryObs = tk.Entry(self, textvariable=self.variableObs, width=30)
        #self.entryParameter  = tk.Entry(self, textvariable=self.variableParameter , width=30)
        
        self.buttonFileEvt = tk.Button(self, text=u"Browse", command=self.onButtonEvtClick)
        self.buttonFileObs = tk.Button(self, text=u"Browse", command=self.onButtonObsClick)
        #self.buttonFileParameter  = tk.Button(self, text=u"Browse", command=self.onButtonParameter Click)
        self.buttonLaunch = tk.Button(self, text=u"Launch", command=self.onButtonLaunchClick)
        self.buttonReset = tk.Button(self, text=u"Reset", command=self.onButtonResetClick, state='disabled')
        
        self.labelQuestion = tk.Label(self, text="Do you want to test one event or all ?", bg='white')
        self.labelEvid = tk.Label(self, text="Enter the ID of the event.", bg='white')
        self.variableEvid = tk.IntVar()
        self.entryEvid = tk.Entry(self, textvariable=self.variableEvid, width=10, state='readonly')
        self.buttonNotAll = tk.Button(self, text=u"One", command=self.onButtonNotAllClick, state='disabled')
        self.buttonAll = tk.Button(self, text=u"All", command=self.onButtonAllClick, width=5, state='disabled')
        
        self.labelDepthMin = tk.Label(self, text="Do you want to fix a minimal depth ?", bg='white')
        self.labelDepthMax = tk.Label(self, text="Do you want to fix a maximal depth ?", bg='white')
        self.variableDepthMin = tk.IntVar()
        self.variableDepthMax = tk.IntVar()
        self.entryDepthMin = tk.Entry(self, textvariable=self.variableDepthMin, width=10, state='readonly')
        self.entryDepthMax = tk.Entry(self, textvariable=self.variableDepthMax, width=10, state='readonly')
        self.buttonDepthMinYes = tk.Button(self, text="Yes", command=self.onButtonDepthMinYesClick, state='disabled')
        self.buttonDepthMaxYes = tk.Button(self, text="Yes", command=self.onButtonDepthMaxYesClick, state='disabled')
        self.buttonDepthMinNo = tk.Button(self, text="No", command=self.onButtonDepthMinNoClick, state='disabled')
        self.buttonDepthMaxNo = tk.Button(self, text="No", command=self.onButtonDepthMaxNoClick, state='disabled')
        
        self.labelIc = tk.Label(self, text="Do you want to fix an intensity of completeness ?", bg='white')
        self.variableIc = tk.IntVar()
        self.entryIc = tk.Entry(self, textvariable=self.variableIc, width=10, state='readonly')
        self.buttonIcYes = tk.Button(self, text="Yes", command=self.onButtonIcYesClick, state='disabled')
        self.buttonIcNo = tk.Button(self, text="No", command=self.onButtonIcNoClick, state='disabled')
        
        self.labelDirectory = tk.Label(self, text="In which directory do you want to save your data ?", bg='white')
        self.variableDirectory = tk.StringVar()
        self.entryDirectory = tk.Entry(self, textvariable=self.variableDirectory, width=20, state='readonly')
        self.buttonDirectory = tk.Button(self, text="Browse", command=self.onButtonDirectoryClick, state='disabled')
        
        self.labelEqFile = tk.Label(self, text="Select the files with the equation", bg='white')
        self.labelCoeff = tk.Label(self, text="Select the rating of each file", bg='white')
        self.labelIbin = tk.Label(self, text="Intensity bin. strategy", bg='white')
        
        self.listeVariableEq = []
        self.listeVariableCoeff = []
        self.listeEntryEq = []
        self.listeEntryCoeff = []
        self.listeBrowse = []
        self.liste_binning = []
        
        self.buttonAdd = tk.Button(self, text="Add", command=self.addLineEq, state='disabled')
        self.buttonDelete = tk.Button(self, text="Delete", command=self.deleteLineEq, state='disabled')
        
        self.buttonStart = tk.Button(self, text="Start", command=self.onButtonStartClick, state='disabled', width=10, height=2)
        
        self.xDefilB = tk.Scrollbar(self, orient='horizontal', command=self.canvas.xview)
        self.yDefilB = tk.Scrollbar(self, orient='vertical', command=self.canvas.yview)
        
        self.canvas['xscrollcommand'] = self.xDefilB.set
        self.canvas['yscrollcommand'] = self.yDefilB.set
        # Place Widgets on the GUI
        self.grid()
        self.canvas.grid(column=0, row=0, sticky='nsew')
        self.xDefilB.grid(column=0, row=1, sticky='ew')
        self.yDefilB.grid(column=1, row=0, sticky='ns')     
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.canvas.create_window(80, 20, window=self.labelFileEvt, anchor='ne')
        self.canvas.create_window(80, 60, window=self.labelFileObs, anchor='ne')
        #self.canvas.create_window(80, 100, window=self.labelFileParameter , anchor='ne')
        
        self.canvas.create_window(90, 20, window=self.entryEvt, anchor='nw')
        self.canvas.create_window(90, 60, window=self.entryObs, anchor='nw')
        #self.canvas.create_window(90, 100, window=self.entryParameter , anchor='nw')
        
        self.canvas.create_window(360, 16, window=self.buttonFileEvt, anchor='n')
        self.canvas.create_window(360, 56, window=self.buttonFileObs, anchor='n')
        #self.canvas.create_window(360, 96, window=self.buttonFileParameter , anchor='n')
        self.canvas.create_window(360, 100, window=self.buttonLaunch, anchor='n')
        self.canvas.create_window(320, 100, window=self.buttonReset, anchor='ne')
        
        self.canvas.create_window(210, 190, window=self.labelQuestion, anchor='n')
        self.canvas.create_window(120, 210, window=self.labelEvid, anchor='n')
        self.canvas.create_window(125, 240, window=self.entryEvid, anchor='ne')
        self.canvas.create_window(135, 236, window=self.buttonNotAll, anchor='nw')
        self.canvas.create_window(310, 236, window=self.buttonAll, anchor='n')
        
        self.canvas.create_window(210, 290, window=self.labelDepthMin, anchor='n')
        self.canvas.create_window(125, 320, window=self.entryDepthMin, anchor='ne')
        self.canvas.create_window(135, 316, window=self.buttonDepthMinYes, anchor='nw')
        self.canvas.create_window(310, 316, window=self.buttonDepthMinNo, anchor='n')
        
        self.canvas.create_window(210, 370, window=self.labelDepthMax, anchor='n')
        self.canvas.create_window(125, 400, window=self.entryDepthMax, anchor='ne')
        self.canvas.create_window(135, 396, window=self.buttonDepthMaxYes, anchor='nw')
        self.canvas.create_window(310, 396, window=self.buttonDepthMaxNo, anchor='n')
        
        self.canvas.create_window(210, 450, window=self.labelIc, anchor='n')
        self.canvas.create_window(125, 480, window=self.entryIc, anchor='ne')
        self.canvas.create_window(135, 476, window=self.buttonIcYes, anchor='nw')
        self.canvas.create_window(310, 476, window=self.buttonIcNo, anchor='n')
        
        self.canvas.create_window(210, 530, window=self.labelDirectory, anchor='n')
        self.canvas.create_window(230, 560, window=self.entryDirectory, anchor='ne')
        self.canvas.create_window(240, 556, window=self.buttonDirectory, anchor='nw')
        
        self.canvas.create_window(550, 20, window=self.labelEqFile, anchor='n')
        self.canvas.create_window(920, 20, window=self.labelCoeff, anchor='n')
        self.canvas.create_window(750, 20, window=self.labelIbin, anchor='n')
        
        self.canvas.create_window(900, 56, window=self.buttonAdd, anchor='n')
        self.canvas.create_window(820, 56, window=self.buttonDelete, anchor='n')
        self.canvas.create_window(900, 106, window=self.buttonStart, anchor='n')
        self.addLineEq()
        self.listeEntryEq[0].config(state='readonly')
        self.listeEntryCoeff[0].config(state='readonly')
        self.listeBrowse[0].config(state='disabled')
    
    def onButtonEvtClick(self):
        filename = fd.askopenfilename(title="Choose the Event file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.variableEvt.set(filename)
    
    def onButtonObsClick(self):
        filename = fd.askopenfilename(title="Choose the Observation file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.variableObs.set(filename)
        
#    def onButtonCritiqueClick(self):
#        filename = fd.askopenfilename(title="Choose the Ic/Dc file", filetypes=[('txt files','.txt'),('all files','.*')])
#        self.variableCritique.set(filename)
        
    def onButtonLaunchClick(self):
        if self.variableEvt.get() == "" or self.variableObs.get() == "":
            tkm.showerror("Need files","We need at least Event and Observation files")
            return
        try:
            #self.files = peo.FilesForPlot(self.variableEvt.get(), self.variableObs.get(), self.variableCritique.get()) 
            self.files = peo.FilesForPlot(self.variableEvt.get(), self.variableObs.get())
        except:
            return
        
        self.entryEvt.config(state='disabled')
        self.entryObs.config(state='disabled')
        #self.entryParameter.config(state='disabled')
        self.buttonFileEvt.config(state='disabled')
        self.buttonFileObs.config(state='disabled')
        #self.buttonFileParamater.config(state='disabled')
        
        self.buttonLaunch.config(state='disabled')
        self.buttonReset.config(state='normal')
        
        self.entryEvid.config(state='normal')
        self.buttonNotAll.config(state='normal')
        self.buttonAll.config(state='normal')
        
        self.entryDepthMin.config(state='normal')
        self.buttonDepthMinYes.config(state='normal')
        self.buttonDepthMinNo.config(state='normal')
        
        self.entryDepthMax.config(state='normal')
        self.buttonDepthMaxYes.config(state='normal')
        self.buttonDepthMaxNo.config(state='normal')
        
        self.entryIc.config(state='normal')
        self.buttonIcYes.config(state='normal')
        self.buttonIcNo.config(state='normal')
        try:
            self.files.ParameterFile
        except:
            self.buttonIcNo.config(state='disabled')
        
        self.entryDirectory.config(state='normal')
        self.buttonDirectory.config(state='normal')
        
        self.listeEntryEq[0].config(state='normal')
        self.listeEntryCoeff[0].config(state='normal')
        self.listeBrowse[0].config(state='normal')
        
        self.buttonAdd.config(state='normal')
        self.buttonStart.config(state='normal')
        
    def onButtonResetClick(self):      
        self.entryEvt.config(state='normal')
        self.entryObs.config(state='normal')
        #self.entryParameter.config(state='normal')
        
        self.buttonFileEvt.config(state='normal')
        self.buttonFileObs.config(state='normal')
        #self.buttonFileParameter.config(state='normal')
        
        self.buttonLaunch.config(state='normal')
        
        self.variableEvt.set("")
        self.variableObs.set("")
        #self.variableParameter.set("")
        self.variableEvid.set(0)
        self.variableDepthMin.set(0)
        self.variableDepthMax.set(0)
        self.variableIc.set(0)
        self.variableDirectory.set("")
        while len(self.listeVariableEq) > 1:
            self.deleteLineEq()
        self.listeEntryEq[0].config(state='readonly')
        self.listeEntryCoeff[0].config(state='readonly')
        self.listeBrowse[0].config(state='disabled')
           
        self.buttonReset.config(state='disabled')
        
        self.entryEvid.config(state='disabled')
        self.buttonNotAll.config(state='disabled')
        self.buttonAll.config(state='disabled')
        
        self.entryDepthMin.config(state='readonly')
        self.buttonDepthMinYes.config(state='disabled')
        self.buttonDepthMinNo.config(state='disabled')
        
        self.entryDepthMax.config(state='readonly')
        self.buttonDepthMaxYes.config(state='disabled')
        self.buttonDepthMaxNo.config(state='disabled')
        
        self.entryIc.config(state='readonly')
        self.buttonIcYes.config(state='disabled')
        self.buttonIcNo.config(state='disabled')
        
        self.entryDirectory.config(state='disabled')
        self.buttonDirectory.config(state='disabled')
        
        self.buttonAdd.config(state='disabled')
        self.buttonDelete.config(state='disabled')
        self.buttonStart.config(state='disabled')
        
    def onButtonNotAllClick(self):
        self.entryEvid.config(state='normal')
        
    def onButtonAllClick(self):
        self.variableEvid.set(0)
        self.entryEvid.config(state='readonly')
        
    def onButtonDepthMinYesClick(self):
        self.entryDepthMin.config(state='normal')
        
    def onButtonDepthMinNoClick(self):
        self.variableDepthMin.set(0)
        self.entryDepthMin.config(state='readonly')
        
    def onButtonDepthMaxYesClick(self):
        self.entryDepthMax.config(state='normal')
        
    def onButtonDepthMaxNoClick(self):
        self.variableDepthMax.set(0)
        self.entryDepthMax.config(state='readonly')
        
    def onButtonIcYesClick(self):
        self.entryIc.config(state='normal')
        
    def onButtonIcNoClick(self):
        self.variableIc.set(0)
        self.entryIc.config(state='readonly')
        
    def onButtonDirectoryClick(self):
        directory = fd.askdirectory()
        self.variableDirectory.set(directory)
        
    def addLineEq(self):
        liste_Ibin = ['RAVG', 'ROBS', 'RP50', 'RP84', 'RF50', 'RF84']
        if not (len(self.listeVariableEq) == len(self.listeVariableCoeff) and 
            len(self.listeVariableCoeff) == len(self.listeEntryCoeff) and
            len(self.listeEntryCoeff) == len(self.listeEntryEq) and
            len(self.listeEntryEq) == len(self.listeBrowse)):
            tkm.showerror("Error", "Problem with counting lines")
            return
        index = len(self.listeVariableEq)
        variableEq = tk.StringVar()
        variableCoeff = tk.DoubleVar()
        variableIbinchoice = tk.StringVar()
        entryEq = tk.Entry(self, textvariable=variableEq, width=22)
        entryCoeff = tk.Entry(self, textvariable=variableCoeff, width=10)
        comboIbin = tk.ttk.Combobox(self, values=liste_Ibin, textvariable=variableIbinchoice,
                                    state='readonly',
                                    width=10)
        browse = tk.Button(self, text="Browse")
        
        self.listeVariableEq.append(variableEq)
        self.listeVariableCoeff.append(variableCoeff)
        self.listeEntryEq.append(entryEq)
        self.listeEntryCoeff.append(entryCoeff)
        self.listeBrowse.append(browse)
        self.liste_binning.append(comboIbin)
        
        def onBrowseButtonClickAux(evt, i=index):
            return self.onBrowseButtonClick(evt, i)
        browse.bind('<Button-1>', onBrowseButtonClickAux)
        
        h = len(self.listeVariableEq) * 40 + 20
        """
        self.canvas.create_window(550, 20, window=self.labelEqFile, anchor='n')
        self.canvas.create_window(920, 20, window=self.labelCoeff, anchor='n')
        self.canvas.create_window(750, 20, window=self.labelIbin, anchor='n')
        """
        self.canvas.create_window(550, h, window=entryEq, anchor='n')
        self.canvas.create_window(650, h - 4, window=browse, anchor='n')
        self.canvas.create_window(750, h - 4, window=comboIbin, anchor='n')
        self.canvas.create_window(900, h, window=entryCoeff, anchor='n')
        self.canvas.create_window(900, h + 36, window=self.buttonAdd, anchor='n')
        self.canvas.create_window(820, h + 36, window=self.buttonDelete, anchor='n')
        self.canvas.create_window(900, h + 86, window=self.buttonStart, anchor='n')
        if len(self.listeVariableEq) > 1:
            self.buttonDelete.config(state='normal')
        self.replaceTheScrollbar()
          
    def onBrowseButtonClick(self, evt, index):
        filename = fd.askopenfilename(title="Choose the Equation file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.listeVariableEq[index].set(filename)
    
    def deleteLineEq(self):
        if not (len(self.listeVariableEq) == len(self.listeVariableCoeff) and 
            len(self.listeVariableCoeff) == len(self.listeEntryCoeff) and
            len(self.listeEntryCoeff) == len(self.listeEntryEq) and
            len(self.listeEntryEq) == len(self.listeBrowse)):
            tkm.showerror("Error", "Problem with counting lines")
            return
        lastindex = len(self.listeVariableEq) - 1
        
        self.listeEntryEq[lastindex].destroy()
        self.listeEntryCoeff[lastindex].destroy()
        self.listeBrowse[lastindex].destroy()
        self.liste_binning[lastindex].destroy()
        
        self.listeVariableEq = self.listeVariableEq[:lastindex]
        self.liste_binning = self.liste_binning[:lastindex]
        self.listeVariableCoeff = self.listeVariableCoeff[:lastindex]
        self.listeEntryEq = self.listeEntryEq[:lastindex]
        self.listeEntryCoeff = self.listeEntryCoeff[:lastindex]
        self.listeBrowse = self.listeBrowse[:lastindex]
        
        h = len(self.listeVariableEq) * 40 + 56
        self.canvas.create_window(900, h, window=self.buttonAdd, anchor='n')
        self.canvas.create_window(820, h, window=self.buttonDelete, anchor='n')
        self.canvas.create_window(900, h + 50, window=self.buttonStart, anchor='n')
        if len(self.listeVariableEq) <= 1:
            self.buttonDelete.config(state='disabled')
        
    def onButtonStartClick(self):
        if (self.variableIc.get() < 0) or (self.variableDepthMin.get() < 0) or (self.variableDepthMax.get() < 0):
             tkm.showerror("Error", "Problem with Ic, Depth Min/max values")
             return
         
        sumVC = 0
        for i in range(len(self.listeVariableCoeff)):
            sumVC += self.listeVariableCoeff[i].get()
         
        if not sumVC == 1:
            tkm.showerror("Error", "The sum of rates must be 1")
            return
            
        config = self.entryIc.config()
        print('State of entryIc: ' + config['state'][4])
        if config['state'][4] == 'normal':
            ic = self.variableIc.get()
        else:
            ic = False
        
        config = self.entryDepthMin.config()
        if config['state'][4] == 'normal':
            depthmin = self.variableDepthMin.get()
        else:
            depthmin = False
        
        config = self.entryDepthMax.config()
        if config['state'][4] == 'normal':
            depthmax = self.variableDepthMax.get()
        else:
            depthmax = False
        
        print("Ic : " + str(ic))
        print("Depth Min : " + str(depthmin))
        print("Depth Max : " + str(depthmax))
        
        queue = que.Queue()
        config = self.entryEvid.config()
        if config['state'][4] == 'normal':
            try:
                self.evt
                if self.evt.evid == self.variableEvid.get():
                    queue.put(self.evt)
                else: 
                    evt = peo.PlotEvt(self.files)
                    evt.build(self.variableEvid.get())
                    queue.put(evt)
            except:
                evt = peo.PlotEvt(self.files)
                evt.build(self.variableEvid.get())
                queue.put(evt)
        elif config['state'][4] == 'readonly':
            evt = peo.PlotEvt(self.files)
            liste_evt = self.files.EvtFile['EVID'].values
            for evt_id in liste_evt:
                evt = peo.PlotEvt(self.files)
                evt.build(evt_id)
                queue.put(evt)
        
        listeEq = []
        for i in range(len(self.listeVariableEq)):
            listeEq.append(self.listeVariableEq[i].get())
            
        print(listeEq)
        
        listeIbin = []
        for i in range(len(self.liste_binning)):
            listeIbin.append(self.liste_binning[i].get())
        print(listeIbin)
        listeCoeff = []
        for i in range(len(self.listeVariableCoeff)):
            listeCoeff.append(self.listeVariableCoeff[i].get())
        
        print(listeCoeff)
        
        
        folder = self.variableDirectory.get()
        if folder == "":
            tkm.showerror("Error", "Need a folder")
        # Ajouter listeIbin
        thread = QUakeThread(queue, self.variableEvt.get(), self.variableObs.get(), 
                     listeEq, listeCoeff, listeIbin, ic, 
                     folder, depthmin,
                     depthmax)
        
        thread.start()
        
#        qumd.QUakeMD(queue, self.variableEvt.get(), self.variableObs.get(), 
#                     self.variableCritique.get(), listeEq, listeCoeff, Ic=ic, 
#                     output_folder=folder, depth_min_ref=depthmin,
#                     depth_max_ref=depthmax)
        

        pass
        
    def replaceTheScrollbar(self):
        self.xDefilB.destroy()
        self.yDefilB.destroy()
        
        self.xDefilB = tk.Scrollbar(self, orient='horizontal', command=self.canvas.xview)
        self.yDefilB = tk.Scrollbar(self, orient='vertical', command=self.canvas.yview)
        
        self.canvas['xscrollcommand'] = self.xDefilB.set
        self.canvas['yscrollcommand'] = self.yDefilB.set
        
        self.xDefilB.grid(column=0, row=1, sticky='ew')
        self.yDefilB.grid(column=1, row=0, sticky='ns')
        
    def build(self, EvtName, ObsName,  files, event, evid, Parametername=""):
        #self.createView()
        self.variableEvt.set(EvtName)
        self.variableObs.set(ObsName)
        #self.variableCritique.set(CritiqueName)
        
        self.files = files
        self.evt = event

        self.entryEvt.config(state='disabled')       
        self.entryObs.config(state='disabled')
        #self.entryCritique.config(state='disabled')
        self.buttonFileEvt.config(state='disabled')
        self.buttonFileObs.config(state='disabled')
        #self.buttonFileCritique.config(state='disabled')
        
        self.buttonLaunch.config(state='disabled')
        self.buttonReset.config(state='normal')
        
        self.entryEvid.config(state='normal')
        self.buttonNotAll.config(state='normal')
        self.buttonAll.config(state='normal')
        
        self.entryDepthMin.config(state='normal')
        self.buttonDepthMinYes.config(state='normal')
        self.buttonDepthMinNo.config(state='normal')
        
        self.entryDepthMax.config(state='normal')
        self.buttonDepthMaxYes.config(state='normal')
        self.buttonDepthMaxNo.config(state='normal')
        
        self.entryIc.config(state='normal')
        self.buttonIcYes.config(state='normal')
        self.buttonIcNo.config(state='normal')
#        try:
#            self.files.CritiqueFile
#        except:
#            self.buttonIcNo.config(state='disabled')
            
        self.entryDirectory.config(state='normal')
        self.buttonDirectory.config(state='normal')
        
        self.listeEntryEq[0].config(state='normal')
        self.liste_binning[0].config(state='normal')
        self.listeEntryCoeff[0].config(state='normal')
        self.listeBrowse[0].config(state='normal')
        self.buttonAdd.config(state='normal')
        self.buttonStart.config(state='normal')
        
        self.variableEvid.set(evid)
        
class QUakeThread(th.Thread):
    def __init__(self, queue, Evtname, Obsname,  listeEq, listeCoeff, listeIbin, 
                 ic, folder, depthmin, depthmax, Parametername=""):
        th.Thread.__init__(self)
        self.queue = queue
        self.Evtname = Evtname
        self.Obsname = Obsname
        self.Parametername = Parametername
        self.listeEq = listeEq
        self.listeCoeff = listeCoeff
        self.listeIbin = listeIbin
        self.ic = ic
        self.folder = folder
        self.depthmin = depthmin
        self.depthmax = depthmax
        
    def run(self):
        qumd.QUakeMD(self.queue, self.Evtname, self.Obsname, self.listeEq, self.listeCoeff, self.listeIbin,
                     Ic=self.ic, 
                     output_folder=self.folder, depth_min_ref=self.depthmin,
                     depth_max_ref=self.depthmax, Parametername=self.Parametername)
        tkm.showinfo("QUake-MD","The calculus of the magnitude and the depth is over")
    
if __name__ == "__main__":
    app = QUakeMdInterface(None)
    app.title('QUake-MD')
    app.mainloop()
