# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:11:43 2019

@author: AdminLocal
"""

#à garder tant que l'algo n'échelle est pas fait
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class plotEvtInterfaceGraph(tk.Tk):
    def __init__(self,parent, fig):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        canvas = FigureCanvasTkAgg(fig, self)
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        
        self.title("Graph")