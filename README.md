# gnm
A network model of glymphatic flow
This repository contains the MATLAB code to set up and solve a hydraulic network model for glymphatic flow described by Titof et al. (https://www.biorxiv.org/content/10.1101/2021.09.23.461519v1). 

setParam.m -- sets the variables that are used in the model 
FindP.m -- a function that iteratively solves for the total driving pressure applied to the model in order to have a median flow in the pial pvs of 18.7 microns/s
branching_hexagon_model_pext.m -- the main code that creates and solves the model
frref.m -- a function called by the model code to solve the matrix. It is from the MATLAB file exchange
bhm_plot.m -- a function for visualizing the 3D network model and the results from the model
find_stats.m -- a function called by bhm_plot.m
