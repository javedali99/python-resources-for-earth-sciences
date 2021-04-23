
<!--# Python Resources for Earth Sciences-->
![featured](https://user-images.githubusercontent.com/15319503/115857198-63b72d80-a3fb-11eb-82bc-a37b77e7521c.png)


<a href="https://twitter.com/javedali99"><img src="https://img.shields.io/badge/twitter-%231DA1F2.svg?&style=for-the-badge&logo=twitter&logoColor=white" alt="Twitter@javedali"></a>


<div>

This repository contains information about python libraries broadly relevant to Earth Sciences (Hydrology, Meteorology, Geospatial, Climatology etc.). Packages are broadly grouped according to their function; however, many have functionality that spans multiple categories. 

If you have any comments or suggestions for additions or improvements for this repository, submit an issue , or make some changes and submit a pull request. If you can’t contribute on GitHub, [send me an email](mailto:javedali28@gmail.com). If you have an issue with one of the libraries discussed below, please contact the maintainer of that library.

---

## Content

* [Geospatial Analysis and Mapping](#geospatial-analysis-and-mapping)
* [Hydrology](#hydrology)
  + [Data Collection](#data-collection)
  + [Hydrological Modelling](#hydrological-modelling)
  + [Groundwater Modelling](#groundwater-modelling)
  + [Time Series Analysis](#time-series-analysis)
  + [Optimization, Uncertainty, Statistics](#optimization-uncertainty-statistics)
* [Meteorology](#meteorology)
* [Climatology](#climatology)
* [Seismology](#seismology)


---

# Geospatial Analysis and Mapping

- [Geopandas](https://geopandas.org/index.html): GeoPandas is an open source project to make working with geospatial data in python easier. GeoPandas extends the datatypes used by pandas to allow spatial operations on geometric types.

- [whitebox](https://github.com/giswqs/whitebox-python): The whitebox Python package is built on WhiteboxTools, an advanced geospatial data analysis platform. WhiteboxTools can be used to perform common geographical information systems (GIS) analysis operations, such as cost-distance analysis, distance buffering, and raster reclassification.

- [PySal](https://pysal.org/): A python spatial analysis library for open source and crossed platform Geospatial Data Science 

- [Shapely](https://shapely.readthedocs.io/en/latest/): Shapely is a Python package for manipulation and analysis of planar geometric objects. It is based on the widely deployed `GEOS` (the engine of PostGIS) and `JTS` libraries.

- [Rasterio](https://rasterio.readthedocs.io/en/latest/): Rasterio is a `GDAL` and Numpy-based Python library designed to make your work with geospatial raster data more productive, more fun — more Zen. It is a highly useful module for raster processing which you can use for reading and writing several different raster formats in Python. Python automatically registers all known `GDAL` drivers for reading supported formats when importing the module.

- [Georasters](https://pypi.org/project/georasters/): The GeoRasters package is a python module that provides a fast and flexible tool to work with GIS raster files. It provides the GeoRaster class, which makes working with rasters quite transparent and easy

- [Fiona](https://github.com/Toblerity/Fiona): It reads and writes geographic data files and thereby helps Python programmers integrate geographic information systems with other computer systems.

- [geemap](https://geemap.org/): A Python package for interactive mapping with Google Earth Engine, `ipyleaflet`, and `ipywidgets`.

- [Awesome Earth Engine](https://github.com/giswqs/Awesome-GEE): A curated list of Google Earth Engine resources includng many python libraries 

- [geonotebook](https://github.com/OpenGeoscience/geonotebook): Jupyter notebook extension for geospatial visualization and analysis developed by NASA

- [Verde](https://github.com/fatiando/verde): It is a Python library for processing spatial data (bathymetry, geophysics surveys, etc) and interpolating it on regular grids (i.e., gridding).

- [PcRaster](http://pcraster.geo.uu.nl/): It is a collection of software targeted at the development and deployment of spatio-temporal environmental models. 

- [PyGeoprocessing](https://pypi.org/project/pygeoprocessing/): A Python/Cython based library that provides a set of commonly used raster, vector, and hydrological operations for GIS processing. 

- [Pysheds](https://github.com/mdbartos/pysheds): Simple and fast watershed delineation in python.

- [Lidar](https://github.com/giswqs/lidar): Terrain and hydrological analysis based on LiDAR-derived digital elevation models (DEM).

- [PYWR](https://github.com/pywr/pywr): Spatial allocation tool 

- [ParTerra-Python](https://gitlab.com/deltares/parterra/parterra-python): The Participatory Terrain model deploys an algorithm to fuse together data from OpenStreetMap (OSM) and any base elevation dataset to create a high-resolution digital terrain model for any area in the world.

- [salem](https://salem.readthedocs.io/en/latest/): Adds geolocalised subsetting, masking, and plotting operations to xarray's data structures via accessors

- [Regionmask](https://regionmask.readthedocs.io/en/stable/): Plotting and creation of masks of spatial regions

- [xshape](https://xshape.readthedocs.io/en/latest/): Tools for working with shapefiles, topographies, and polygons in xarray

- [Collocate](https://github.com/cistools/collocate): Collocate xarray trajectories in arbitrary physical dimensions

- [HoloViews](http://holoviews.org/): Library designed to make data analysis and visualization seamless and simple

- [GeoViews](http://geoviews.org/): Library that makes it easy to explore and visualize geographical, meteorological, and oceanographic datasets, such as those used in weather, climate, and remote sensing research

- [Datashader](https://github.com/pyviz/datashader): Graphics pipeline system for creating meaningful representations of large datasets quickly and flexibly

- [Panel](https://panel.pyviz.org/): Create custom interactive web apps and dashboards by connecting user-defined widgets to plots, images, tables, or text

- [hvPlot](https://hvplot.pyviz.org/): A high-level plotting API for the PyData ecosystem built on HoloViews

- [EarthSim](https://earthsim.pyviz.org/): Tools for working with and visualizing environmental simulations

- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/): Easy cartographic (maps) data visualization.

- [Geoviews](http://geo.holoviews.org/): Explore and visualize geographic data using HoloViews.



</br>

# Hydrology

### Data Collection

- [HKVFEWSPY](https://github.com/HKV-products-services/hkvfewspy): Connection to the DelftFEWS servers 

- [HyRiver](https://github.com/cheginit/HyRiver): HyRiver is a software stack consisting of six Python libraries that are designed to aid in watershed analysis through web services. This project includes hydrology and climatology data within the US.

- [Openradar](https://github.com/nens/openradar): Library for processing a set of dutch, german and belgian precipitation radars into calibrated composites. 

- [Ecohydrolib](https://github.com/selimnairb/EcohydroLib): Libraries and command-line scripts for performing ecohydrology data preparation workflows. 

- [Ulmo](https://github.com/ulmo-dev/ulmo/): Clean, simple and fast access to public hydrology and climatology data.

- [PyHIS](https://pypi.org/project/pyhis/): It is a python library for querying CUAHSI*-HIS** web services

- [Wetterdienst](https://github.com/earthobservations/wetterdienst): Python Toolset For Accessing Weather Data From German Weather Service

- [ERA5-tools](https://github.com/loicduffar/ERA5-tools): Python scripts to download and view ERA5 climatologic data, as well as to extract time series (hourly to monthly data on many atmospheric and land-surface parameters)

- [CAMELS-AUS](https://github.com/csiro-hydroinformatics/camels-aus-py): Python package to easily load and use the CAMELS-AUS dataset

- [BoM Water](https://github.com/csiro-hydroinformatics/pybomwater): This package has been developed to access to the BoM Water Data Sensor Observation Service (SOS). With a goal to easily and efficiently integrate data into scientific workflows.

- [Eco-Data Manage Toolkit](https://github.com/cayolopesbc/eco-data-toolkit): It is a Python toolkit to facilitate data management for hydrology/limnology applications.



### Hydrological Modelling

- [CMF](https://github.com/philippkraft/cmf): Catchment Modelling Framework, a hydrologic modelling toolbox.

- [TopoFlow](https://github.com/peckhams/topoflow): Spatial hydrologic model (D8-based, fully BMI-compliant).

- [VIC](https://github.com/UW-Hydro/VIC): The Variable Infiltration Capacity (VIC) Macroscale Hydrologic Model. 

- [Xanthos](https://github.com/JGCRI/xanthos): Xanthos is an open-source hydrologic model, written in Python, designed to quantify and analyze global water availability.

- [WRF-Hydro](https://github.com/NCAR/wrf_hydro_py): wrfhydrpy is a Python API for the WRF-Hydro modelling system. 

- [pyDEM](https://github.com/creare-com/pydem): PyDEM is a package for topographic (terrain) analysis. It takes in digital elevation model (DEM) rasters, and it outputs quantities like slope, aspect, upstream area, and topographic wetness index.

- [EXP-HYDRO](https://github.com/sopanpatil/exp-hydro): EXP-HYDRO is a catchment scale hydrological model that operates at a daily time-step.  It takes as inputs the daily values of precipitation, air temperature, and potential evapotranspiration, and simulates daily streamflow at the catchment outlet. 

- [RRMPG](https://github.com/kratzert/RRMPG): Rainfall-Runoff modelling playground. 

- [LHMP](https://github.com/hydrogo/LHMP): Lumped Hydrological Models Playground. 

- [SMARTPy](https://github.com/ThibHlln/smartpy): Python implementation of the rainfall-runoff model SMART 

- [PyStream](https://github.com/martibosch/pystream): Python implementation of the STREAM hydrological rainfall-runoff model. 

- [HydrPy](https://github.com/hydpy-dev/hydpy): A framework for the development and application of hydrological models based on Python. 

- [Catchmod](https://pypi.org/project/pycatchmod/): CATCHMOD is widely used rainfall runoff model in the United Kingdom. It was introduced by Wilby (1994).

- [wflow](https://github.com/openstreams/wflow): wflow consists of a set of Python programs that can be run on the command line and perform hydrological simulations. The models are based on the PCRaster Python framework 

- [PyTOPKAPI](https://github.com/sahg/PyTOPKAPI): PyTOPKAPI is a BSD licensed Python library implementing the TOPKAPI Hydrological model (Liu and Todini, 2002).

- [mhmpy](https://github.com/MuellerSeb/mhmpy): A Python-API for the mesoscale Hydrological Model.

- [SuperflexPy](https://github.com/dalmo1991/superflexPy): A new open source framework for building conceptual hydrological models 

- [NeuralHydrology](https://github.com/neuralhydrology/neuralhydrology): Python library to train neural networks with a strong focus on hydrological applications

- [StreamStats](https://github.com/earthlab/streamstats): Python package for interfacing with the USGS StreamStats API.

- [hidrocomp](https://github.com/clebsonpy/HidroComp): Python library for hydrological data analysis

- [PyFlo](https://github.com/benjiyamin/pyflo): It is an open-source Python library for performing hydraulic and hydrology stormwater analysis. Features include network hydraulic grade analysis and time/iteration based storage and flood routing simulations.

- [HydroFunctions](https://github.com/mroberge/hydrofunctions): A suite of convenience functions for working with hydrology data in an interactive Python session.

- [pySTEPS](https://github.com/pySTEPS/pysteps): It is an open-source and community-driven Python library for probabilistic precipitation nowcasting, i.e. short-term ensemble prediction systems.


### Groundwater Modelling

- [Flopy](https://github.com/modflowpy/flopy): The Python interface to MODFLOW. 

- [imod-python](https://imod.xyz/): Make massive MODFLOW models. 

- [Idfpy](https://github.com/tomvansteijn/idfpy): A simple module for reading and writing iMOD IDF files. IDF is a simple binary format used by the iMOD groundwater modelling software. 

- [WellApplication](https://github.com/utah-geological-survey/WellApplication): Set of tools for groundwater level and water chemistry analysis. 

- [TIMML](https://github.com/mbakker7/timml):  A Multi-Layer, Analytic Element Model. 

- [TTim](https://github.com/mbakker7/ttim): A Multi-Layer, Transient, Analytic Element Model.  

- [PyHELP](https://github.com/cgq-qgc/pyhelp): A Python library for the assessment of spatially distributed groundwater recharge and hydrological components with HELP. 

- [PyRecharge](https://github.com/abdikaiym/pyrecharge): Spatially distributed groundwater recharge and depletion modeling framework in Python

- [Anaflow](https://github.com/GeoStat-Framework/AnaFlow): A python-package containing analytical solutions for the groundwater flow equation

- [WellTestPy](https://github.com/GeoStat-Framework/welltestpy): A python-package for handling well based field campaigns.

- [HydroGeoSines](https://github.com/HydroGeoSines/HydroGeoSines): Signal In the Noise Exploration Software for Hydrogeological Datasets.

- [Pytesmo](https://github.com/TUW-GEO/pytesmo): Python Toolbox for the Evaluation of Soil Moisture Observations.

- [Phydrus](https://github.com/phydrus/phydrus): Python implementation of the HYDRUS-1D unsaturated zone model


### Time Series Analysis

- [Hydropy](https://github.com/stijnvanhoey/hydropy): Analysis of hydrological oriented time series. 

- [Pastas](https://github.com/pastas/pastas): Analysis of hydrological time series using time series models. 

- [Hydrostats](https://github.com/BYU-Hydroinformatics/Hydrostats): Tools for use in comparison studies, specifically for use in the field of hydrology. 

- [htimeseries](https://github.com/openmeteo/htimeseries)| This module provides the HTimeseries class, which is a layer on top of pandas, offering a little more functionality.

- [efts-python](https://github.com/csiro-hydroinformatics/efts-python): A python library for reading and writing Ensemble Forecast Time Series in netCDF files.


### Optimization, Uncertainty, Statistics

- [LMFIT](https://github.com/lmfit/lmfit-py): Non-Linear Least Squares Minimization, with flexible Parameter settings, based on scipy.optimize.leastsq, and with many additional classes and methods for curve fitting. 
 
- [SPOTpy](https://github.com/thouska/spotpy): A Statistical Parameter Optimization Tool for Python. 

- [PyGLUE](http://code.activestate.com/pypm/pyglue/): Generalised Likelihood Uncertainty Estimation (GLUE) Framework. 

- [Pyemu](https://github.com/jtwhite79/pyemu): A python modules for model-independent uncertainty analyses, data-worth analyses, and interfacing with PEST(++). 

- [HPGL](http://hpgl.github.io/hpgl/): High Performance Geostatistics Library. 

- [HydroErr](https://github.com/BYU-Hydroinformatics/HydroErr): Goodness of Fit metrics for use in comparison studies, specifically in the field of hydrology. 

- [Climate-indices](https://github.com/monocongo/climate_indices): Climate indices for drought monitoring, community reference implementations in Python. 

- [HydroLM](https://github.com/mullenkamp/HydroLM): The HydroLM package contains a class and functions for automating linear regressions OLS for hydrologists.

- [PySDI](https://bitbucket.org/pysdi/pysdi/src/master/): It is a set of open source scripts that compute non-parametric standardized drought indices (SDI) using raster data sets as input data.

- [PyForecast](https://github.com/usbr/PyForecast): It is a statistical modeling tool useful in predicting monthly and seasonal inflows and streamflows. The tool collects meterological and hydrologic datasets, analyzes hundreds to thousands of predictor subsets, and returns statistical regressions between predictors and streamflows.


### Miscellaneous

- [ESMPY](https://www.earthsystemcog.org/projects/esmpy/): Earth System Modeling Framework (ESMF) Python interface 

- [PyHSPF](https://github.com/djlampert/PyHSPF): Python extensions to the Hydrological Simulation Program in Fortran (HSPF).

- [SPHY](https://github.com/WilcoTerink/SPHY): Spatial Processes in HYdrology (SPHY) model 

- [xsboringen](https://github.com/tomvansteijn/xsboringen): (In Dutch) A python library for processing and plotting borehole and CPT data, developed for open data formats in the Netherlands. 

- [PyMT](https://github.com/csdms/pymt/): It is an Open Source Python package that provides the necessary tools used for the coupling of models that expose the Basic Model Interface (BMI).

- [Landlab](https://github.com/landlab/landlab): The Landlab project creates an environment in which scientists can build a numerical landscape model without having to code all of the individual components. 

- [EFlowCalc](https://github.com/ThibHlln/eflowcalc): Calculator of Streamflow Characteristics. 

- [IRIS](https://github.com/SciTools/iris): A powerful, format-agnostic, and community-driven Python library for analysing and visualising Earth science data. 

- [Hydrointerp](https://github.com/mullenkamp/hydrointerp): A Python package for interpolating hydrologic data. 
 
- [EFlowCalc](https://github.com/ThibHlln/eflowcalc): EFlowCalc is an open-source calculator of ecological streamflow characteristics in Python. 
 
- [Hydrofunctions](https://github.com/mroberge/hydrofunctions): A suite of convenience functions for working with hydrology data in an interactive Python session.
 
- [Shyft](https://gitlab.com/shyft-os): It is the open-source toolbox for the energy-market domain, funded and supported by Statkraft. 

- [Hydroshare](https://github.com/hydroshare/hydroshare): HydroShare is a collaborative website for better access to data and models in the hydrologic sciences.

- [Hydrobox](https://github.com/VForWaTer/hydrobox): Hydrological preprocessing and analysis toolbox build upon pandas and numpy

- [Wetland](https://github.com/giswqs/wetland): It is a toolset for mapping surface water and wetland hydrological dynamics using fine-resolution aerial imagery within Google Earth Engine (GEE).

- [iRONS](https://github.com/AndresPenuela/iRONS): iRONS (interactive Reservoir Operation Notebooks and Software) is a python package that enables the simulation, forecasting and optimisation of reservoir systems.

</br>

# Meteorology

- [MetPy](https://github.com/Unidata/MetPy): It is a collection of tools in Python for reading, visualizing and performing calculations with weather data. 

- [PyEto](https://github.com/woodcrafty/PyETo): It is a Python library for calculating reference crop evapotranspiration (ETo), sometimes referred to as potential evapotranspiration (PET). The library provides numerous functions for estimating missing meteorological data. 

- [Improver](https://github.com/metoppv/improver): It is a library of algorithms for meteorological post-processing and verification. 

- [MetSim](https://github.com/UW-Hydro/MetSim): It is a meteorological simulator and forcing disaggregator for hydrologic modeling and climate applications. 

- [MELODIST](https://github.com/kristianfoerster/melodist): It is an open-source toolbox written in Python for disaggregating daily meteorological time series to hourly time steps. 

- [PyCat](https://github.com/wegener-center/pyCAT): Climate Analysis Tool written in python 

- [PySteps](https://github.com/pySTEPS/pysteps): It is a community-driven initiative for developing and maintaining an easy to use, modular, free and open source Python framework for short-term ensemble prediction systems. 

- [Evaporation](https://github.com/openmeteo/evaporation): Calculation of evaporation and transpiration. 

- [rainymotion](https://github.com/hydrogo/rainymotion): Python library for radar-based precipitation nowcasting based on optical flow techniques.

- [Metview](https://github.com/ecmwf/metview-python): Python interface to Metview, a meteorological workstation and batch system for accessing, examining, manipulating and visualising meteorological data.

- [IMPROVER](https://github.com/metoppv/improver): It is a library of algorithms for meteorological post-processing.

- [JAMS](https://github.com/mcuntz/jams_python): It is a general Python package offering miscellaneous functions in different categories, such as reading different file formats, julian date routines, or meteorological functions.

- [windspharm](https://ajdawson.github.io/windspharm/index.html): Python Spherical harmonic wind analysis

- [wrf-python](https://wrf-python.readthedocs.io/en/latest/): Python A collection of diagnostic and interpolation routines for use with output of the Weather Research and Forecasting (WRF-ARW) Model

- [scikit-downscale](https://github.com/pangeo-data/scikit-downscale): Statistical downscaling and postprocessing models for climate and weather model simulations.

- [Awesome-EarthObservation-Code](https://github.com/acgeospatial/awesome-earthobservation-code): A curated list of awesome tools, tutorials, code, helpful projects, links, stuff about Earth Observation and Geospatial stuff!


</br>

# Climatology

- [climlab](https://climlab.readthedocs.io/en/latest/): Process-oriented climate modeling

- [aospy](https://aospy.readthedocs.io/en/stable/): Automated analysis and management of gridded climate data

- [Oocgcm](https://oocgcm.readthedocs.io/en/latest/): Analysis of large gridded geophysical datasets

- [Pangaea](https://pangaea.readthedocs.io/en/latest/): xarray extension for gridded land surface & weather model output

- [xgcm](https://xgcm.readthedocs.io/en/latest/): Extends the xarray data model to understand finite volume grid cells (common in General Circulation Models) and provides interpolation and difference operations for such grids

- [OpenClimateGIS](https://www.earthsystemcog.org/projects/openclimategis/): Geospatial manipulation, subsetting, computation, and translation of spatiotemporal climate data

- [climpred](https://climpred.readthedocs.io/en/stable/): xarray wrapper for analysis of ensemble forecast models for climate prediction

- [pyOWM](https://github.com/csparpa/pyowm): PyOWM is a client Python wrapper library for OpenWeatherMap (OWM) web APIs


</br>

# Seismology

- [Madagascar](http://www.ahay.org/wiki/Installation): Multi-dimensional data processing suite

- [ObsPy](https://github.com/obspy/obspy/wiki): Framework for reading, writing and processing seismic and seismological data

- [Bruges](https://github.com/agile-geoscience/bruges/tree/master/bruges): Various geophysical equations and tools

- [Segyio](https://github.com/equinor/segyio): Fast library for seismic SEGY files

- [Pyrocko](https://github.com/pyrocko/pyrocko): Seismology toolkit

- [rsudp](https://github.com/raspishake/rsudp): Continuous ObsPy-based visual display, sudden motion monitoring, and historical replay of Raspberry Shake data

- [SeismicZFP](https://github.com/equinor/seismic-zfp): Convert SEG-Y/ZGY files to compressed [SGZ files](https://github.com/equinor/seismic-zfp/blob/master/docs/file-specification.md) & retrieve arbitrary sub-volumes from these, fast





