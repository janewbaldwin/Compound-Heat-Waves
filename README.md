# Compound-Heat-Waves

Code supporting publication in journal Earth's Future, Baldwin et al 2019 "Temporally Compound Heat Waves and Global Warming: An Emerging Hazard."

The data corresponding to this code can be downloaded from the following ftp site-- ftp://nomads.gfdl.noaa.gov/users/Jane.Baldwin/compoundheatwaves/GFDL-CM2.5-FLOR/

Code is divided into two main folders /PaperFigures and /HeatWaveStats.
/PaperFigures includes scripts corresponding to every figure in the paper. Most are Jupyter Notebooks, and some notebooks contain multiple figures. Figure 1 is created using a FERRET script (.jnl). Figure S8 is created using a typical python script (.py). Most of the figures that use derived heat wave statistics in the paper use the tn_90_3114 definition (daily minimum temperature, 90th percentile threshold, temporal structure 311) as a representative example. The derived heat wave statistics for the whole suite of definitions are available on the ftp site, and can also be recreated using the /HeatWaveStats code to analyze the raw climate model output, as described below.

/HeatWaveStats includes scripts to create the derived heat wave statistics for various heat wave definitions from the raw GCM output. Due to the computational demands of analyzing the global, daily temperature data for this task, these scripts were originally run on the GFDL analysis cluster. As such, two different types of scripts are included. Scripts without file extensions are PBS scripts which were used to submit and schedule the heat waves statistics scripts on the GFDL analysis cluster. These scripts show how to run the .py scripts that actually analyze the heat wave statistics. The code is divided between /ENSEMBLE, /CONT_2XCO2_MEANSHIFTS, and /MERRA2 which produces derived heat wave statistics from the FLOR ensemble; the Control, mean shifted Control, and 2XCO2 simulations; and the MERRA2 data respectively.

