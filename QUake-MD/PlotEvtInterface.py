# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 10:33:43 2019

@author: baize-funck-ame
"""

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox as tkm
import tkinter.filedialog as fd
import PlotEvtObject as peo
import QUakeMDInterface as qumdi
#import PlotEvtInterfaceGraph as peig
import threading as th

class plotEvtInterface(tk.Toplevel):
    def __init__(self, parent):
        tk.Toplevel.__init__(self, parent)
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        self.createView()
        self.createController()
    
    def createView(self):
        self.canvas = tk.Canvas(self, width=1500, height=750, scrollregion=(0,0,1500,750), bg='white')
        
        self.labelFileEvt = tk.Label(self, text="Evt file", bg='white')
        self.labelFileObs = tk.Label(self, text="Obs file", bg='white')
        self.labelFileCritique = tk.Label(self, text="Ic/Dc file", bg='white')
        
        self.variableEvt = tk.StringVar()
        self.variableObs = tk.StringVar()
        self.variableCritique = tk.StringVar()
        self.entryEvt = tk.Entry(self, textvariable=self.variableEvt, width=30)
        self.entryObs = tk.Entry(self, textvariable=self.variableObs, width=30)
        self.entryCritique = tk.Entry(self, textvariable=self.variableCritique, width=30)
        
        self.buttonFileEvt = tk.Button(self, text=u"Browse", command=self.onButtonEvtClick)
        self.buttonFileObs = tk.Button(self, text=u"Browse", command=self.onButtonObsClick)
        self.buttonFileParameter = tk.Button(self, text=u"Browse", command=self.onButtonParameterClick)
        self.buttonLaunch = tk.Button(self, text=u"Launch", command=self.onButtonLaunchClick)
        self.buttonReset = tk.Button(self, text=u"Reset", command=self.onButtonResetClick, state='disabled')
        
        #Par evid
        labelEvid = tk.Label(self, text="Which evid do you want to visualize ?", bg='white')
        self.evid =  tk.IntVar()
        self.entryEvid = tk.Entry(self, textvariable=self.evid, width=13, state='readonly')
        self.buttonEvid = tk.Button(self, text="OK", command=self.onButtonEvidClick, state='disabled')
        
        #Par date
        labelDate = tk.Label(self, text="When did the earthquake happen ?", bg='white')
        self.day = tk.IntVar()
        self.month = tk.IntVar()
        self.year = tk.IntVar()
        self.entryDay = tk.Entry(self, textvariable=self.day, width=5, state='readonly')
        self.entryMonth = tk.Entry(self, textvariable=self.month, width=5, state='readonly')
        self.entryYear = tk.Entry(self, textvariable=self.year, width=5, state='readonly')
        self.buttonDate = tk.Button(self, text="OK", command=self.onButtonDateClick, state='disabled')
        
        
#        labelLoc = tk.Label(self, text="Where did the earthquake happen ?", bg='white')
#        self.variableLon = tk.DoubleVar()
#        self.variableLat = tk.DoubleVar()
#        self.variableRad = tk.DoubleVar()
#        self.entryLon = tk.Entry(self, textvariable=self.variableLon, width=5)
#        self.entryLat = tk.Entry(self, textvariable=self.variableLat, width=5)
#        self.entryRad = tk.Entry(self, textvariable=self.variableRad, width=5)
#        self.buttonLoc = tk.Button(self, command=self.onButtonLocClick, text="OK")
        
        self.liste = fd.Listbox(self, width=15, state='disabled')
        self.buttonChoice = tk.Button(self, text="OK", command=self.onButtonChoiceClick, state='disabled')
        
        self.buttonQUakeMD = tk.Button(self, text="Start QUake-MD", command=self.onButtonStartClick, state='disabled')
        
        self.variableNameFig = tk.StringVar()
        self.entryNameFig = tk.Entry(self, textvariable=self.variableNameFig, state='disabled')
        self.buttonSave = tk.Button(self, text="Save", command=self.onButtonSaveClick, state='disabled')
        self.figPlt = 0
        self.leg0 = 0
        self.plot = 0
        
        self.xDefilB = tk.Scrollbar(self, orient='horizontal', command=self.canvas.xview)
        self.yDefilB = tk.Scrollbar(self, orient='vertical', command=self.canvas.yview)
        
        self.canvas['xscrollcommand'] = self.xDefilB.set
        self.canvas['yscrollcommand'] = self.yDefilB.set
        
        # Place Widgets    
        self.grid()
        self.canvas.grid(column=0, row=0, sticky='nsew')
        self.xDefilB.grid(column=0, row=1, sticky='ew')
        self.yDefilB.grid(column=1, row=0, sticky='ns')     
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.canvas.create_window(80, 20, window=self.labelFileEvt, anchor='ne')
        self.canvas.create_window(80, 60, window=self.labelFileObs, anchor='ne')
        #self.canvas.create_window(80, 100, window=self.labelFileParameter, anchor='ne')
        
        self.canvas.create_window(90, 20, window=self.entryEvt, anchor='nw')
        self.canvas.create_window(90, 60, window=self.entryObs, anchor='nw')
       #self.canvas.create_window(90, 100, window=self.entryParameter, anchor='nw')
        
        self.canvas.create_window(360, 16, window=self.buttonFileEvt, anchor='n')
        self.canvas.create_window(360, 56, window=self.buttonFileObs, anchor='n')
        #self.canvas.create_window(360, 96, window=self.buttonFileParameter, anchor='n')
        self.canvas.create_window(360, 96, window=self.buttonLaunch, anchor='n')
        self.canvas.create_window(320, 96, window=self.buttonReset, anchor='ne')
        
        self.canvas.create_window(5, 180, window=labelEvid, anchor='nw')
        self.canvas.create_window(230, 180, window=self.entryEvid, anchor='nw')
        self.canvas.create_window(360, 176, window=self.buttonEvid, anchor='n')
        
        self.canvas.create_window(5, 220, window=labelDate, anchor='nw')
        self.canvas.create_window(210, 220, window=self.entryDay, anchor='nw')
        self.canvas.create_window(250, 220, window=self.entryMonth, anchor='nw')
        self.canvas.create_window(290, 220, window=self.entryYear, anchor='nw')
        self.canvas.create_window(360, 216, window=self.buttonDate, anchor='n')
        
#        self.canvas.create_window(5, 260, window=labelLoc, anchor='nw')
#        self.canvas.create_window(210, 260, window=self.entryLon, anchor='nw')
#        self.canvas.create_window(250, 260, window=self.entryLat, anchor='nw')
#        self.canvas.create_window(290, 260, window=self.entryRad, anchor='nw')
#        self.canvas.create_window(360, 256, window=self.buttonLoc, anchor='n')
#        Add 50 to the ordinate below if you add the localisation
        
        self.canvas.create_window(150, 260, window=self.liste, anchor='nw')
        self.canvas.create_window(300, 330, window=self.buttonChoice)
        
        self.canvas.create_window(200, 470, window=self.buttonQUakeMD, anchor='n')
        
        
        
        self.canvas.create_window(650, 625, window=self.entryNameFig, anchor='nw')
        self.canvas.create_window(800, 621, window=self.buttonSave, anchor='nw')

    def onButtonEvtClick(self):
        filename = fd.askopenfilename(title="Choose the Event file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.variableEvt.set(filename)
    
    def onButtonObsClick(self):
        filename = fd.askopenfilename(title="Choose the Observation file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.variableObs.set(filename)
        
    def onButtonCritiqueClick(self):
        filename = fd.askopenfilename(title="Choose the Ic/Dc file", filetypes=[('txt files','.txt'),('all files','.*')])
        self.variableCritique.set(filename)
        
    def onButtonLaunchClick(self):
        if self.variableEvt.get() == "" or self.variableObs.get() == "":
            tkm.showerror("Need files","We need at least Event and Observation files")
            return
        try:
            self.files = peo.FilesForPlot(self.variableEvt.get(), self.variableObs.get(), self.variableCritique.get())  
        except:
            return
        self.entryEvt.config(state='disabled')
        self.entryObs.config(state='disabled')
        self.entryCritique.config(state='disabled')
        self.buttonFileEvt.config(state='disabled')
        self.buttonFileObs.config(state='disabled')
        self.buttonFileCritique.config(state='disabled')
        
        self.buttonLaunch.config(state='disabled')
        self.buttonReset.config(state='normal')
        
        self.entryEvid.config(state='normal')
        self.buttonEvid.config(state='normal')    
        self.entryDay.config(state='normal')
        self.entryMonth.config(state='normal')
        self.entryYear.config(state='normal')
        self.buttonDate.config(state='normal')
        
    def onButtonResetClick(self):      
        self.entryEvt.config(state='normal')
        self.entryObs.config(state='normal')
        self.entryCritique.config(state='normal')
        
        self.buttonFileEvt.config(state='normal')
        self.buttonFileObs.config(state='normal')
        self.buttonFileCritique.config(state='normal')
        
        self.buttonLaunch.config(state='normal')
        
        self.variableEvt.set("")
        self.variableObs.set("")
        self.variableCritique.set("")
        self.evid.set(0)
        self.day.set(0)
        self.month.set(0)
        self.year.set(0)
        
        if self.liste.size() > 0:
            self.liste.delete(0, self.liste.size()-1)
            
        self.buttonReset.config(state='disabled')
        
        self.entryEvid.config(state='disabled')
        self.buttonEvid.config(state='disabled')
        
        self.entryDay.config(state='disabled')
        self.entryMonth.config(state='disabled')
        self.entryYear.config(state='disabled')
        self.buttonDate.config(state='disabled')
        
        self.liste.config(state='disabled')
        self.buttonChoice.config(state='disabled')
        
        self.buttonQUakeMD.config(state='disabled')
        
    def onButtonEvidClick(self):
        try:
            self.evid.get()
        except:
            tkm.showerror('Error', "Please, enter an integer")
        self.UseThread()
                
    def onButtonDateClick(self):
        self.liste.config(state='normal')
        self.buttonChoice.config(state='normal')
        
        try:
            day = self.day.get()
        except:
            tkm.showerror('Error', "Please, an integer for the day")
        try:
            month = self.month.get()
        except:
            tkm.showerror('Error', "Please, an integer for the month")
        try:
            year = self.year.get()
        except:
            tkm.showerror('Error', "Please, an integer for the year")
        
        setEvid = peo.searchByDate(self.variableEvt.get(), day, month, year)
        
        if self.liste.size() > 0:
            self.liste.delete(0, self.liste.size()-1)
            
        for i in range(len(setEvid)):
            self.liste.insert(i, setEvid[i])
        
    def onButtonLocClick(self):
        self.liste.config(state='normal')
        self.buttonChoice.config(state='normal')
        
        try:
            variableLon = self.variableLon.get()
        except:
            tkm.showerror('Error', "Please, a numeric value for the longitude")
        try:
            variableLat = self.variableLat.get()
        except:
            tkm.showerror('Error', "Please, an numeric value for the latitude")
        try:
            variableRad = self.variableRad.get()
        except:
            tkm.showerror('Error', "Please, an numeric value for the radium")
        
        setEvid = peo.searchByLoc(self.variableEvt.get(), variableLon, variableLat, variableRad)
        
        if self.liste.size() > 0:
            self.liste.delete(0, self.liste.size()-1)
            
        for i in range(len(setEvid)):
            self.liste.insert(i, setEvid[i])
            
    def onButtonChoiceClick(self):    
        choiceInt = self.liste.curselection()
        if not len(choiceInt) == 1:
            tkm.showinfo("Information","You have to select one event")
            print("You have to select one event")
            return
        
        self.evid.set(self.liste.get(choiceInt[0]))
        print(self.evid)
        self.UseThread()    
        
    def onButtonSaveClick(self):
        if self.variableNameFig.get() == "":
            self.variableNameFig.set("Evt_" + str(self.evid.get()))
            
        #faire un petit thread
        self.figPlt.savefig(self.variableNameFig.get() + '.png', dpi=150, bbox_inches='tight', bbox_extra_artist=[self.leg0])
        
    def onButtonStartClick(self):
        quake = qumdi.QUakeMdInterface(None)
        quake.build(self.variableEvt.get(), self.variableObs.get(), self.variableCritique.get(), self.files, self.plot.deepCopy(), self.evid.get())
        quake.title('QUake-MD')
        
    
    def UseThread(self):
        useTh = WorkThread(self)
        useTh.start()
        #self.progressbar.start()
        #useTh.join()
        
    def onPressEvidEnter(self, event):
        self.onButtonEvidClick()
    
    def onPressYearEnter(self, event):
        self.onButtonDateClick()
        
    def onPressListEnter(self, event):
        self.onButtonChoiceClick()
            
    def createController(self):
        self.entryEvid.bind("<Return>", self.onPressEvidEnter)
        self.entryYear.bind("<Return>", self.onPressYearEnter)
        self.liste.bind("<Return>", self.onPressListEnter)

    
class WorkThread(th.Thread):
    def __init__(self, pInterface):
        th.Thread.__init__(self)
        self.pInterface = pInterface
    def run(self):
        self.pInterface.entryNameFig.config(state='disabled')
        self.pInterface.buttonSave.config(state='disabled')
        self.pInterface.buttonQUakeMD.config(state='disabled')
        
        self.pInterface.labelprogressbar = tk.Label(self.pInterface, text='Plotting the data:', bg='white')
        self.pInterface.progressbar = ttk.Progressbar(self.pInterface, orient='horizontal', mode='indeterminate')
        self.pInterface.canvas.create_window(200, 600, window=self.pInterface.labelprogressbar, anchor='n')
        self.pInterface.canvas.create_window(200, 620, window=self.pInterface.progressbar, anchor='n')
        self.pInterface.progressbar.start()
        
        self.pInterface.plot = peo.PlotEvt(self.pInterface.files)
        self.pInterface.plot.build(self.pInterface.evid.get())
        
        self.pInterface.figPlt = self.pInterface.plot.getFig()
        self.pInterface.leg0 = self.pInterface.plot.leg0
        
        fig = FigureCanvasTkAgg(self.pInterface.figPlt, self.pInterface)
        
        self.pInterface.figure = fig.get_tk_widget()
        
        self.pInterface.canvas.create_window(450, 1, window=self.pInterface.figure, anchor='nw')
            
        self.pInterface.progressbar.stop()
        self.pInterface.labelprogressbar.destroy()
        self.pInterface.progressbar.destroy()
        
        self.pInterface.entryNameFig.config(state='normal')
        self.pInterface.buttonSave.config(state='normal')
        self.pInterface.buttonQUakeMD.config(state='normal')
        self.replaceTheScrollbar()
        
    def replaceTheScrollbar(self):
        
        self.pInterface.xDefilB.destroy()
        self.pInterface.yDefilB.destroy()
        
        self.pInterface.xDefilB = tk.Scrollbar(self.pInterface, orient='horizontal', command=self.pInterface.canvas.xview)
        self.pInterface.yDefilB = tk.Scrollbar(self.pInterface, orient='vertical', command=self.pInterface.canvas.yview)
        
        self.pInterface.canvas['xscrollcommand'] = self.pInterface.xDefilB.set
        self.pInterface.canvas['yscrollcommand'] = self.pInterface.yDefilB.set
        
        self.pInterface.xDefilB.grid(column=0, row=1, sticky='ew')
        self.pInterface.yDefilB.grid(column=1, row=0, sticky='ns') 
        


if __name__ == "__main__":
    app = plotEvtInterface(None)
    app.title('Data Visualisation')
    app.mainloop()