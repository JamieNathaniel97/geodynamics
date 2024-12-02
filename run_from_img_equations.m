%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target grid size z-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% calculate sediment parameters
% set variables (values found online)
rho_particle = 2650; % density of sand, gravel & silt
rho_air  = 1.225;    % density of air
sa_phi   = 0.253;    % sand porosity fraction
gr_phi   = 0.325;    % gravel porosity fraction
si_phi   = 0.17;     % silt porosity fraction

sa_particle_kT = 8;      % sand particle thermal conductivity 
gr_particle_kT = 5;      % gravel particle thermal conductivity
si_particle_kT = 3.5;    % silt particle thermal conductivity
air_kT         = 0.025;  % air thermal conductivity

% calculate bulk density 
rho_sand   = sa_phi * rho_air + (1-sa_phi) * rho_particle;
rho_gravel = gr_phi * rho_air + (1-gr_phi) * rho_particle;
rho_silt   = si_phi * rho_air + (1-si_phi) * rho_particle;

% calculate effective thermal conductivity. Idk how to make this non
% negative
sa_kT = sa_particle_kT * (air_kT + (1 - sa_phi) * (sa_particle_kT - air_kT))/...
                         (air_kT + (1 - sa_phi) * (air_kT - sa_particle_kT));
gr_kT = gr_particle_kT * (air_kT + (1 - gr_phi) * (gr_particle_kT - air_kT))/...
                         (air_kT + (1 - sa_phi) * (air_kT - sa_particle_kT));
si_kT = si_particle_kT * (air_kT + (1 - si_phi) * (si_particle_kT - air_kT))/...
                         (air_kT + (1 - sa_phi) * (air_kT - sa_particle_kT));

% material properties for each rock unit (update based on your calibration)
matprop = [
% unit  conductivity (kT)   density (rho)   heat capacity (Cp)   heat production (Hr)
   1	      3.678            2697.6	           845	                4.172    %HE1
   2	      1	               2000	               1000	                1        %Bg
   3	      sa_kT          rho_sand              1000	                1        %sand
   4	      3.218            2703.5	           845	                5.575    %HE2
   5	      gr_kT          rho_gravel            1000	                1        %gravel
   6	      1	               2000	               1000	                1        %clay
   7	      si_kT          rho_silt              1000	                1        %silt
   8	      1	               2000	               1000	                1        %Ms
   9	      1e-6             1000	               1000	                0];      % air/water
       
% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3),Nz,Nx);
Cp     = reshape(matprop(units,4),Nz,Nx);
kT     = reshape(matprop(units,2),Nz,Nx);
Hr     = reshape(matprop(units,5),Nz,Nx);

% continue setting remaining model parameters, then call model routine

