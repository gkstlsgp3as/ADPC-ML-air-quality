**Interpolation-Based Fusion of Sentinel-5P, SRTM, and Regulatory-Grade Ground Stations Data for Producing Spatially Continuous Maps of PM2.5 Concentrations Nationwide over Thailand**
==========

**Abstract**
Atmospheric pollution has recently drawn significant attention due to its proven adverse effects on public health and the environment. This concern has been aggravated specifically in Southeast Asia due to increasing vehicular use, industrial activity, and agricultural burning practices. Consequently, elevated PM2.5 concentrations have become a matter of intervention for national authorities who have addressed the needs of monitoring air pollution by operating ground stations. However, their spatial coverage is limited and the installation and maintenance are costly. Therefore, alternative approaches are necessary at national and regional scales. In the current paper, we investigated interpolation models to fuse PM2.5 measurements from ground stations and satellite data in an attempt to produce spatially continuous maps of PM2.5 nationwide over Thailand. Four approaches are compared, namely the inverse distance weighted (IDW), ordinary kriging (OK), random forest (RF), and random forest combined with OK (RFK) leveraging on the NO2, SO2, CO, HCHO, AI, and O3 products from the Sentinel-5P satellite, regulatory-grade ground PM2.5 measurements, and topographic parameters.

![adpc_mdpi](https://github.com/gkstlsgp3as/ADPC-ML-air-quality/assets/58411517/63465d92-4cb7-4784-907d-dfad4b9c4771)


## Codes

These are the brief explanation about each file.
- Task 6: Process S5P data-final.ipynb: get Sentinel-5P data and convert it to gridded files
- s5p_pre_processor.R: make stacks of raster files and extract values for calculating correlation.
- interpolation_method.R: interpolate PM2.5 data with IDW, OK, RF, and RK method, and validate the performance.
- Notepad_Sienna.pdf: cover all the work that was done during the internship.


For details see 
[atmosphere-13-00161-v2.pdf](https://github.com/gkstlsgp3as/ADPC-ML-air-quality/files/12116355/atmosphere-13-00161-v2.pdf)
