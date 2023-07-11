.. QUake-MD Interface UserManual documentation master file, created by
   sphinx-quickstart on Mon Jul 10 12:47:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QUake-MD Interface's documentation!
=========================================================

**QUake-MD GUI** creates a graphical interface between a user and the **QUake-MD** tool.

The QUake-MD tool estimates a weighted space of magnitude (M), depth (H) and epicentral intensity 
(I0) based on historical intensity data points (IDP), their associated uncertainties and empirical 
intensity prediction equations (IPE). The final space of solution aims to be representative of the IDPs 
and their quality and of IPE epistemic uncertainties. 

QUake-MD GUI offers in addition to a QUake-MD interface a data visualization tool. The IDP 
visualization is an important step for M/H/I0 estimates from macroseismic data. It is strongly 
recommended to visualize the macroseismic data before M/H/I0 estimates.

The user manual describes first how to run QUake-MD GUI, then how to use the data visualization 
part of QUake-MD GUI and finally how to use the QUake-MD part.

.. note::

   The QUake-MD methodology is published in Provost and Scotti 2020:
   Ludmila Provost, Oona Scotti; QUake‐MD: Open‐Source Code to Quantify Uncertainties in Magnitude–Depth Estimates of Earthquakes from Macroseismic Intensities.
   Seismological Research Letters 2020;; 91 (5): 2520–2530. doi: https://doi.org/10.1785/0220200064

.. toctree::
   :maxdepth: 3
   :caption: Contents:
   
   installation
   launch



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
