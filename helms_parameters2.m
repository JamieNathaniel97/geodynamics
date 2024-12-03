%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;


% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target grid size z-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit (update based on your calibration)

matprop = [
% unit conductivity(kT) density(rho)  heat capacity(Cp)  heat production(Hr)
   1	      3.678       2697.6             1000	          4.172         %HE1
   2	      1	          2000	             1000	          1             %Gneiss
   3	      1	          2000	             1000	          1             %Sand
   4	      3.218       2703.5             1000	          5.575         %HE2
   5	      1	          2000	             1000	          1             %Gravel
   6	      1	          2000	             1000	          1             %Clay
   7	      1	          2000	             1000	          1             %Silt
   8	      1	          2000	             1000	          1             %Mud
   9	      1e-6        1000	             1000	          0];           % air/water
         
% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3),Nz,Nx); %density
Cp     = reshape(matprop(units,4),Nz,Nx); %heat capacity
kT     = reshape(matprop(units,2),Nz,Nx); %heat conductivity
Hr     = reshape(matprop(units,5),Nz,Nx); %radiogenic heating rate

%Test in constant coefficient unit value case (repeat for other variables)
%rho = 2400*ones(Nz,Nx);

% continue setting remaining model parameters, then call model routine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TGrad = 35;            %Reference Temperature Gradient from National Survey
Ttop  = 5;            % surface temperature
Tbot  = Ttop + TGrad*(D/1000);           %Temperature at the bottom of the domain using reference gradient


rho0  = 1000;         % reference density [kg/m3]
kT0   = 1e-7;         % reference heat diffusivity [m2/s]
cT    = 1e-9;         % kT T-dependence prefactor
mT    = 2;            % kT T-dependence powerlaw
g0    = 9.8;           % gravity [m/s2]
aT    = 1e-4;         % thermal expansivity [1/C]

ADVN  = 'WENO5';      % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')

yr    = 3600*24*365;  % seconds per year [s]
tend  = 1e4*yr;       % stopping time [s]
CFL   = 1/2;          % Time step limiter
nop   = 10;           % output figure produced every 'nop' steps
alpha = 0.99;         % iterative step size limiter
beta  = 0.95;         % iterative lag parameter
tol   = 1e-8;         % residual tolerance
nup   = 100;          % update T, check residual every nup iterations


run('./helmsdale_take_2.m');
