.. QUake-MD_Interface_UserManual documentation master file, created by
   sphinx-quickstart on Mon Jul 10 12:47:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

QUake-MD Interface's documentation
==================================

Launch QUake-MD Interface
-------------------------

To run QUake-MD GUI, open an ipython terminal/Konsole, move to the QUake-MD GUI repertory and enter::

	run Quake-MD/AppliInterface.py


The following start window will open:

.. image:: accueil_QUakeMD.PNG
    :width: 300px
    :align: center
    :height: 150px
    :alt: alternate text

To access to the data visualization tool, click the *Data Visualization* button. 
To access to the QUake-MD tool, click the *QUake-MD* button.

.. note::

   The QUake-MD window can open from the Data Visualization window. 

Data Visualization
------------------

On the *Data Visualization* window, the user can select a macroseismic database, 
and then select one event from this database by event ID or date of the earthquake occurrence.

Once the event selected, the software will draw a macroseismic map with the IDP and a plot of the 
epicentral projection of the IDP. Additional information will be also projected on the window, like the 
epicentral intensity (I0) quality and the epicenter localization quality.

One example of *Data Visualization* :

.. image:: fenetre_DataVisualization.PNG
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text

Input data format
^^^^^^^^^^^^^^^^^

Different data are needed to visualize the macroseismic data associated to an historical earthquake. 
First data are the epicentral parameter, i.e. longitude and latitude of the epicenter, the epicentral 
intensity, quality associated to epicentral localization and intensity. IDPs, i.e. intensity values 
associated to a locality are also needed, meaning an intensity value, quality of this intensity value and 
the longitude and latitude of the associated locality. An earthquake ID is also necessary to link all the 
data. 

These data have to be stored in two separate files: one Event file and one Observation file. These 
two files and their format are described in the two following paragraphs.


Event file
""""""""""

The Event file contains data which describe the epicenter parameter of historical earthquakes, i.e. the epicenter localisation, time of earthquake occurrence and value of the intensity at the epicenter. 
The first line contain the names of each column. Each following line corresponds to one earthquake.

The mandatory columns are: 
 * The ID of the event (name of the column: EVID),
 * The macroseismic epicentral intensity (name of the column: I0) ,.
 * The quality of the previous element, graded A (good quality), B (middle quality), C (poor quality) or E (very poor quality) (name of the column: QI0).
 * The longitude of the epicenter in WGS84 (name of the column: Lon)
 * The latitude of the epicenter in WGS84 (name of the column: Lat)
 * The quality of the position, graded A (very good location quality), B (good location quality), C (middle location quality), D (poor location quality), E (very poor location quality) or I (location very uncertain, set arbitrary based on one IDP) (name of the column: QPos)
 * The day of the earthquake (name of the column: Day)
 * The month of the earthquake (name of the column: Month)
 * The year of the earthquake (name of the column: Year)
 * A name, surrounded by quotation marks and used for legibility (name of the column: Name). 

Each column is separated by the ; sign. One example of input Event file can be found in the **Example** folder of 
**QUake-MD** repository (Evt.example.txt).


Observation file
""""""""""""""""

The Observation file contains all observations by localities associated to the earthquakes stored in the 
Event file. The first line contain the names of each column. Each following line describes an observation at one locality so many lines can be 
associated to one earthquake.
 
The mandatory columns are: 
 * The ID of the event (name of the column: EVID),
 * The value of intensity at this IDP. It’s equalled to -1 when the earthquake was felt but no value could be attributed (lack of information) (name of the column: IObs),
 * The quality of the intensity, graded A(good quality), B (middle quality) or C (poor quality) (name of the column: QIobs),
 * The longitude of the locality in WGS84 (name of the column: Lon),
 * The latitude of the locality in WGS84 (name of the column: Lat).

Each column is separated by the ; sign. One example of input Event file can be found in the **Example** folder of 
**QUake-MD** repository (Obs.example.txt).

Using Data Visualization interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To use this interface, you have to select the Event and the Observation file names:

 .. image:: Browse_inputdatafiles.PNG
    :width: 280px
    :align: center
    :height: 100px
    :alt: alternate text

To begin, you have to select the 2 files, entering the path in the text zone directly or by clicking on the Browse button and choosing the corresponding file.
Then, just click on the Launch button. If a problem is detected while the files are opened, an error window will appear.

Once files are launched, you can select your event. 
You can directly typing its ID in the white field in front of the *Which event do you want to visualize?* sentence. Then click on the corresponding Ok button or press Enter.
You can also search it by date. You can enter a complete date or a couple month/year or just the year. For that, make sure other fields are completed with a 0. 
Then, click on the Ok button or press Enter. A list with the corresponding earthquakes ID will appearbelow in the dedicated field. Click on the one you want to visualize, then click on the Ok 
button in front of the earthquakes list or press Enter.

 .. image:: DV_selectEQ.PNG
    :width: 250px
    :align: center
    :height: 400px
    :alt: alternate text

Now, wait some seconds, calculus take a while. The map and the graph will appear on the window, or an error will occur and you will be adverted.

 .. image:: DV_map_figure.PNG
    :width: 500px
    :align: center
    :height: 500px
    :alt: alternate text

The output is a figure with left the macroseismic map for the selected event. Each color corresponds in an intensity level. The legend of the macroseismic map is below the map. Right is an epicentral 
projection of the IDP. The title of the figure indicates the event ID, the year of occurrence, the I0 value, the associated quality, the epicenter location quality, the number of IDP Nobs and the number
Nfelt of Felt testimonies (no intensity level could be attributed to the locality and the earthquake was felt).
 
You can save the figure with the **Save** button.
To launch the **QUake-MD** interface, you can click on the **Start QUake-MD** button.

QUake-MD interface
------------------

**QUake_MD** interface can be launched from the start window or from the **Data Visualization** interface 
window.

With the **QUake-MD** interface, the user selects the intensity data and Intensity Prediction Equations (IPEs) to compute probability density functions (PDF) in magnitude and hypocentral depth from intensity data with
the **QUake-MD** method.
The user can compute PDF for one or several earthquakes.

Input data format
^^^^^^^^^^^^^^^^^

**QUake_MD** will need some inputs files. The intensity data needed are stored in two files: an Event file and an Observation file.
The Observation file has the same format as the input Observation file of **Data Visualization**. The Event file need some more field than the one in
**Data Visualization**. Finally, IPEs files are needed.


Observation file
""""""""""""""""
See **Data Visualization** Observation input file.


Event file
""""""""""

See **Data Visualization** Event input file.

IPE file
""""""""

The IPE files contain the IPEs which will be used in the inversion process of Quake-MD, with the following 
mathematical formulation::

 	I = C1 + C2 M + β log10(Dhypo) + γ Dhypo

With C1 and C2 the magnitude coefficient, M the magnitude, β the geometrical attenuation coefficient, 
log the decimal logarithm operator, Dhypo the hypocentral distance and γ the intrinsic attenuation 
coefficient.

The IPEs files have the following format: one head line with information about the origin of the IPE 
(useful just for the user, Quake-MD GUI will not use this line), one blank line, one line with the 
column’s names, another blank line and 5 columns with the IPEs:
 * The first column is the weight of the IPE
 * The second column is the C1 coefficient
 * The third column is the C2 coefficient
 * The fourth column is the β coefficient
 * The fifth column is the γ coefficient

Each column is separated by tabulation. The sum of the weight column must be equal to one.


Using QUake-MD interface
^^^^^^^^^^^^^^^^^^^^^^^^

Loading the intensity data files
""""""""""""""""""""""""""""""""

To use this interface, you have to select the Event and the Observation file names:

 .. image:: Browse_inputdatafiles.PNG
    :width: 280px
    :align: center
    :height: 100px
    :alt: alternate text

To begin, you have to select the 2 files, entering the path in the text zone directly or by clicking on the Browse button and choosing the corresponding file.
Then, just click on the Launch button. If a problem is detected while the files are opened, an error window will appear.

If QUake_MD interface is launched from the Data Visualization interface window, the intensity files are already loaded. 

Selecting earthquakes and their inversion parameters
""""""""""""""""""""""""""""""""""""""""""""""""""""

Two options are available. **QUake-MD** can run one one earthquake or several earthquakes.
To run one earthquake, press the *One* button and enter the earthquake ID in the dedicated field.

 .. image:: QUakeMD_select_EQ.PNG
    :width: 280px
    :align: center
    :height: 100px
    :alt: alternate text

To run on several earthquakes, press the *All* button. **QUake-MD** compute then PDF for all earthquakes in the Event file.

**QUake-MD** compute PDF in magnitude and depth from intensity data. Depth is inverted within an interval delemited by depth bounds. By default, the depth bounds
are equal to 1 km and 25 km. However the user can change those bounds using the dedicated fields on the interface. Only integers are allowed.
The intensity data smaller than the intensity of completeness are not taken into account for the PDF computation. By default, the intensity of completeness is equal to 3.
However the user can change the intensity of completeness using the dedicated field on the interface.

Selecting the output directory
""""""""""""""""""""""""""""""
An output folder has to be selected to store the results of QUake-MD. 

 .. image:: QUake-MD_outputfolder.PNG
    :width: 280px
    :align: center
    :height: 100px
    :alt: alternate text


Selecting the Intensity Prediction Equation
"""""""""""""""""""""""""""""""""""""""""""
The IPEs files has to be selected through the Browse button and associated to a rating.It is possible to add several IPE files with the Add button. in this case, the sum of the rating must 
be 1. Otherwise, a warning message will appear when starting the QUake-MD computation.
For the moment, the sum of the probabilities of each IPE within the file is not tested by the QUake-MD algorithm. Please check it before running QUake-MD. If the sum of the IPE's probabilities
within a IPE file is not equal to 1, the barycenter value in the output file will be not consistent with the direct median output of each IPE (output file *All_IPEs_classical_results.txt*, see
the dedicated chapter about the outputs). 

 .. image:: QUakeMD_selectIPEs.PNG
    :width: 380px
    :align: center
    :height: 200px
    :alt: alternate text

A binning intensity method should be associated to the IPE file.The different binning intensity methods available are intensity class binning method called RAVG, ROBS, RP50, RP84, RF50 and RF84.
The width of the intensity class of RAVG is equal to 1 and use overlapping windows. For the other indicators, the width of the intensity class is equal to 0.25 (only one value of intensity wihtin an intensity class).
The value of intensity and epicentral distance of each isoseist computed with the RAVG method is equal to the weighted mean of the intensity data points (IDP) within the considered intensity class.
For the other methods, the value of intensity of each isoseist is equal to the intensity value of the considered intensity class.
For the ROBS method, the epicentral distance of each isoseist is equal to the epicentral distance weighted mean of the IDP within the considered intensity class.
For the RP50 method, the epicentral distance of each isoseist is equal to the epicentral distance weighted median of the IDP within the considered intensity class.
For the RP84 method, the epicentral distance of each isoseist is equal to the epicentral distance weighted 84th percentile of the IDP within the considered intensity class.
The RF50 and RF84 method are similar to the RP50 and RP84 methods but only the most reliable and farthest isoseist is kept. Epicentral intensity associated to an epicentral distance equal to 0 is the second isoseist for the RF methods.
The most reliable and farthest isoseist is chosen on the base of three criteria: the number of IDP within the intensity class, the epicentral distance of the isoseist and the weight of the IDP within the intensity class.
The square root of the number of IDP, the sum of the IDP weights and the epicentral distance associated to each isoseist are computed and the highest procduct of those three results gives the most reliable and farthest isoseist.

Run Quake-MD
""""""""""""
 Once all previous steps completed, the **QUake-MD** algorithm can be launch with the *Start* Button (below the IPE selection area).


Output files description
^^^^^^^^^^^^^^^^^^^^^^^^

In the ouput folder, different type of output can be found:
 * A log file,
 * A summary file with the output barycenters and percentiles of the magnitude/depth PDF associated to the input parameters,
 * An event output folder with figures of the PDF and the fit of the IPE prediction to the observed intensities and files with the PDF and the gaussian outputs of the magnitude and depth inversions.

The log file
""""""""""""

The log file which report each step of the calculation and the input parameters. The name of the log file is the time of the calculation start: *yearmonth-day_hour:minute:second.txt*.

 .. image:: outputfile_log.PNG
    :width: 380px
    :align: center
    :height: 200px
    :alt: alternate text


The summary file
""""""""""""""""

The other *.txt* file is the *file_temp_.txt* and contains a summary of the results, i.e. the event ID, the epicentral intensity of the input macroseismic catalogue,
its associated quality, the intensity of completeness used, the barycenter of the magnitude solutions estimated by Quake-MD and its 
associated weighted 16th and 84th percentiles, the barycenter of the depth solutions estimated by **QUake-MD** and its associated weighted 16th and 84th percentiles and the barycenter of the epicentral 
intensity solutions estimated by **QUake-MD** and its associated weighted 16th and 84th percentiles.

 .. image:: outputfile_temp.PNG
    :width: 600px
    :align: center
    :height: 120px
    :alt: alternate text

The earthquake output folder
""""""""""""""""""""""""""""

For each earthquake, an output folder is created. The name of the folder is the ID of the earthquake.
Different files can be found in this folder:
 * an IDP binning file, which contains the intensity bins computed for the inversion of magnitude and depth,
 * an All_IPEs_classical_results.txt file, with the output median and standard deviation defining the output gaussian solution for depth and magnitude,
 * PDF files, in depth/magnitude, depth/epicentral intensity and magnitude/epicentral intensity/depth for the final PDF and in depth/magnitude for each IPE file.

**Files**

The *IDP_binning_.txt* is a table separated by ',' with 5 columns:
 * EVID: ID of the earthquake,
 * Depi: Epicentral intensity of the isoseist,
 * I: intensity value of the isoseist,
 * StdI: Uncertainty associated to the intensity value,
 * StdLogR: Uncertainty associated to the epicentral distance,
 * Ndata: Number of IDP used to compute the isoseist.

The name of the *IDP_binning_.txt* contains the name of the method used to compute the isoseists.

 .. image:: outputfile_Ibin.PNG
    :width: 450px
    :align: center
    :height: 250px
    :alt: alternate text


The *All_IPEs_classical_results.txt* file is a table separated by ',' with 11 columns:
 * NumEvt: ID of the earthquake,
 * Bin_method: the intensity binning method used,
 * C1: C1 coefficient,
 * C2: C2 coefficient,
 * Beta: Beta coefficient,
 * Gamma: Gamma coefficient,
 * Mag: median magnitude  of the gaussian output from the least square inversion,
 * StdM: associated standard deviation of the gaussian output from the least square inversion,
 * H: median depth of the gaussian output from the least square inversion,
 * StdH: associated standard deviation of the gaussian output from the least square inversion,
 * Io: median epicentral intensity of the gaussian output from the least square inversion,

This table contains the direct outputs of the inversion of magnitude and depth for each IPE used, i.e. the median and its associated standard deviation used to describe the output gaussian.

 .. image:: outputfile_classicalresults.PNG
    :width: 420px
    :align: center
    :height: 220px
    :alt: alternate text

The *HM.txt*, *HIo.txt* and *HMIo.txt* files contain the different computed PDF. The first line indicates the earthquake ID, the year of occurence and the epicentral intensity from the catalogue.
The three following lines indicates respectively the value of the PDF epicentral intensity barycenter, the value of the PDF magnitude barycenter and the value of the PDF depth (km) barycenter.
The fifth line is the column names of the following PDF. The last lines are the PDF. The content of the first two columns are indicated by the name of the fifth line. The last column is always 
the value of probability associated to the previous columns.

 .. image:: outputfile_HM.PNG
    :width: 550px
    :align: center
    :height: 350px
    :alt: alternate text

**Figures**
Four type of figures are saved. 

The first one is an *fit_intensity* figure, which compare predicted intensity with observed intensity and display also the direct magnitude and depth outputs
for one IPE file. The upper figure represent the intensity prediction equations (IPEs, curves color‐coded by inverted depth) fit to the binned IDPs (diamonds). Gray points correspond to the IDP, the pink band corresponds to the I0
uncertainty range used to filter the space of M, H solutions, and the red dots correspond to the I0 value at the end of the inversion process for each IPE (for IPEs predicting I0
outside of the accepted boundaries, no red dots are shown).
The lower figure represents associated M, H central solutions of the inversions (dots color‐coded by inverted depth).

 .. image:: 650009_fit_intensity_Law_0_RAVG.jpeg
    :width: 380px
    :align: center
    :height: 500px
    :alt: alternate text

The other figures are PDF figures. Only the final PDF is plotted in the output files. The first one is the magnitude/depth PDF:

 .. image:: HM.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: alternate text 

The second one is the corresponding epicentral intensity/depth PDF:

 .. image:: HIo.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: alternate text 

And the third one is the magnitude/epicentral intensity/depth PDF:

 .. image:: HMIo.png
    :width: 380px
    :align: center
    :height: 220px
    :alt: alternate text


.. toctree::
   :maxdepth: 3
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
