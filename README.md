# gnm
A network model of glymphatic flow
This repository contains the MATLAB code to set up and solve a hydraulic network model for glymphatic flow described by Titof et al. (https://www.biorxiv.org/content/10.1101/2021.09.23.461519v1). 

setParam.m -- sets the variables that are used in the model 
FindP.m -- a function that iteratively solves for the total driving pressure applied to the model in order to have a median flow in the pial pvs of 18.7 microns/s. Requires as input the variables that change. See example code below.
branching_hexagon_model_pext.m -- the main code that creates and solves the model
frref.m -- a function called by branching_hexagon_model_pext.m to solve the matrix. It is from the MATLAB file exchange: https://www.mathworks.com/matlabcentral/fileexchange/21583-fast-reduced-row-echelon-form
bhm_plot.m -- a function for visualizing the 3D network model and the results from the model
find_stats.m -- a function called by bhm_plot.m


Example code to run in the Matlab command window:
model_results_name='model_results';

% set parameters that change
sleep_or_awake='sleep';
paren_type='high_res';
pen_perm=nan;
cap_perm=1.8e-14;
cap_ar=0.07;

% set variables that don't change and create a structure to feed into the model
[param] = setParam(pen_perm,cap_perm,cap_ar,paren_type,sleep_or_awake);

% solve for the driving pressure pext that will result in the median pial velocity being 18.7 microns/s
FindP

% run the model
[Qtotal,Rtotal]=branching_hexagon_model_pext(model_results_name,pext,param)

% plot the results
bhm_plot(model_results_name,'volume_flow_rate','b','o',0)
