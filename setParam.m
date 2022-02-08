function [param] = setParam(pen_perm,cap_perm,cap_ar,paren_type,sleep_or_awake)
% SETPARAM -- Use this function to set parameters used in branching hexagon
% model. The parameters that change can be found in lines 49-71; you can
% comment and uncomment the various sections to set the parameters you'd 
% like to change.

%% Parameters that don't change
% Set number of generations
param.g=9;

% Set PVS-to-artery area ratio
param.K=[1.4 1.4 cap_ar nan];

% Set diameter and length of pial arteries (in m)
param.d_pialart=46*1e-6;
param.l_pial_art=175*1e-6*3;

% Set length of pial to penetrating arteries (the horizontal offshoots from
% the main hexagon)
param.l_pial2pen=175*1e-6;

% Set radius and length of penetrating artery segments (in m)
param.d_penart=11*1e-6;
param.l_penart=1000*1e-6;

% Set radius and length of capillaries
%r_cap=2*1e-6; % from Blinder et al, Nature Neuroscience 2013 (pg 890)
param.r_cap=3*1e-6; % from Miyawaki et al, Nature Communications 2020

% Set number of capillaries branching off of each penetrating artery
param.ncap_per_penart=11;

% Set C_efflux: hydraulic conductivity (1/R) for flow past capillaries and
% parenchyma to ground (i.e., flow through perivenous and other efflux
% routes); we don't know the value, but it's got to be negligible compared
% to capillary/parenchymal flow, so set it to something large
% Units: mL/(min*mmHg)
param.C_efflux = 1;

cf=1.25e-10; %conversion factor for resistance to go from SI to mmHg-min/mL

%% Parameters we're changing----------------------------------------------

if strcmp(paren_type,'high_res')
    % C_par low (high resistance)
    Lp=1.8e-8; %cm^2/(mmHg-s); (hydraulic conductance of parenchyma); K low
    ef_cf=0.003; % endfoot cavity fraction; what percentage of the endfoot wall is gap?
    g=20e-9; % (m) gap width
elseif strcmp(paren_type,'low_res')
    % C_par high (low resistance)
    Lp=6.6e-6; %cm^2/(mmHg-s); (hydraulic conductance of parenchyma); K high
    ef_cf=(1-.63); % endfoot cavity fraction; what percentage of the endfoot wall is gap?
    g=ef_cf*2*pi*0.5*param.d_penart/2.5; % (m) gap width
end

% Awake vs asleep
if strcmp(sleep_or_awake,'awake')
    Lp=Lp/5.5; % uncomment this to simulate wakefulness
end

% Set permeability for pial, penetrating, and capillary PVSs (set to
% nan to model as an open space)
param.kappa=[nan pen_perm cap_perm]; % intermediate porosity for pen PVSs
%    param.kappa=[nan 1.8e-14 1.8e-14]; % porous penetrating PVSs
%    param.kappa=[nan nan 1.8e-14]; % open penetrating PVSs

printR=1; %choose whether or not to display the calculated resistance values

%%-----------------------------------------------------------------------
%  Calculations derived from parameters that change
% ---------------------------------------------------------------------

% Calculate R_par: resistance to flow from the pen nodes, through the 
% parenchyma, to venous side; units: mmHg-min/mL
L_art_to_ven=93e-6; % units m; average distance from an arteriole to venule
L_art_to_art=130e-6; % units m; average distance from an arteriole to arteriole
% Model from S8 in https://www.pnas.org/content/pnas/suppl/2017/08/23/1706942114.DCSupplemental/pnas.201706942SI.pdf?targetid=nameddest%3DSTXT
Lp=Lp/(133*100^2);
R_par=log((1-2*L_art_to_ven/param.d_penart)^2)/(2*pi*Lp*param.l_penart/param.ncap_per_penart); 
R_par=R_par*cf;
if printR, disp(['R_par = ' num2str(R_par,6)]), end
%---------------------------------------------------------------------
% Calculate R_efg: resistance to flow through the endfeet
% gaps into the parenchymal space. Goes in series with R_par (mmHg-min/mL)
mu=993*7e-7; % viscosity (kg /(m-s)); rho=993 kg/m^3 and nu=7e-7 m^2/s
T=.45e-6; % (m) thickness of the endfoot wall
R_efg=12*mu*T/(g^2*ef_cf*2*pi*(param.d_penart/2)*sqrt(1+param.K(2))*(param.l_penart/param.ncap_per_penart));
R_efg=R_efg*cf;
if printR, disp(['R_efg = ' num2str(R_efg,6)]), end
param.C_paren=1/(R_efg+R_par); % mL/(min-mmHg)
if printR, disp(['C_paren = ' num2str(param.C_paren,6)]), end

%-----------------------------------------------------------------------
% Set C_cap: hydraulic conductivity (1/R) for flow through the capillary
% perivascular space; units: mL/(min*mmHg)
R_vas=0.4; % poise/micron^3; this comes from Blinder et al 2013
mu_h2o=0.0091; % poise, viscosity for water, what Blinder et al used
r_vas=2; % microns
mult=4.*(1 - 0.863*exp(-r_vas/14.3) + 27.5*exp(-r_vas/0.351)); % multiplier to correct for the granular nature of blood; in the online methods section of Blinder et al 2013
L_vas=R_vas*pi*r_vas^4/(8*mu_h2o*mult); % microns; this is the length of an equivalent single resistor representing the resistance through the entire capillary bed
if printR, disp(['L_vas = ' num2str(L_vas,6) ' microns']), end
param.l_cap=L_vas*1e-6; % Rescaling factor for consistency with the main model

%-----------------------------------------------------------------------
% Calculate the resistance/conductance for pen and precap PVSs and store
% pen PVSs
r=param.d_penart/2; % m
r2=r*sqrt(param.K(2)+1); % radius of outer wall for circular annulus
mu_CSF=993*7e-7; % CSF viscosity (kg /(m-s)); from rho=993 kg/m^3 * nu=7e-7 m^2/s
if isnan(param.kappa(2))
    R_pen=cf*mu*8.91*param.K(2)^(-2.78)/r^4*param.l_penart/param.ncap_per_penart; % (units mL / (mmHg*min*m))
    param.C_pen=1/R_pen; % (units mmHg*min/mL)
else
    R_pen=cf*mu_CSF/(pi*param.kappa(2)*(r2^2-r^2))*param.l_penart/param.ncap_per_penart; % (units mL / (mmHg*min*m))
    param.C_pen=1/R_pen; % (units mmHg*min/mL)
end
if printR, disp(['R_pen = ' num2str(R_pen,6)]), end
% Precap PVSs
r=param.r_cap; % m
r2=r*sqrt(param.K(3)+1); % radius of outer wall for circular annulus
if isnan(param.kappa(3))
    R_cap=cf*8*mu_CSF*L_vas*1e-6/pi/(r2^4-r^4-(r2^2-r^2)^2/log(r2/r)); % mmHg*min/mL
else
    R_cap=cf*mu_CSF*L_vas*1e-6/(param.kappa(3)*pi*(r2^2-r^2)); % mmHg*min/mL
end
param.C_cap = 1/R_cap; % mL/(min*mmHg)
if printR, disp(['R_cap = ' num2str(R_cap,6)]), end


end