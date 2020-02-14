#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:09:49 2019

@author: Baize-Funck Amelie

Script creating the start interface
"""

import tkinter as tk
import PlotEvtInterface as pei
import QUakeMDInterface as qmdi

class applicationInterface(tk.Tk):
    def __init__(self,parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        canvas = tk.Canvas(self, width=500, height=200, scrollregion=(0,0,1200,750), bg='white')
        canvas.pack()
        
        label = tk.Label(self, text="Use directly QUake-MD or visualize your data before.", bg='white')
        buttonVisualisation = tk.Button(self, text='Data visualization', command=self.onButtonVisualisationClick)
        buttonQUakeMD = tk.Button(self, text='QUake-MD', command=self.onButtonQUakeMDClick)
        
        canvas.create_window(240, 40, window=label)
        canvas.create_window(120, 100, window=buttonVisualisation)
        canvas.create_window(360, 100, window=buttonQUakeMD)
        
    def onButtonVisualisationClick(self):
        visual = pei.plotEvtInterface(None)
        visual.title('Data visualization')
        
    def onButtonQUakeMDClick(self):
        quake = qmdi.QUakeMdInterface(None)
        quake.title('QUake-MD')
        
        
if __name__ == "__main__":
    app = applicationInterface(None)
    app.title('Visualisation and Quake-MD')
    app.mainloop()
