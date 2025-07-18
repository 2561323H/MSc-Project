%************   RUN FILE: GEOTHERMAL PROJECT (VERS. A)  *************

% clear workspace
clear all; close all; %clc;



%*** SETUP PARAMETERS*************************************************

W = 4e3; % width of the domain, should match image  [m]
dx = 40; % spacing between each grid point - resolution  [m]
Nx = W/dx; % number of column in x direction, in this case 800  [m]


%*** IMAGE PARAMETERS*************************************************

[units, D, Nz] = ModelFromImageDiscrete('image_faults_new.tiff', W, Nx);
% units = value of each pixel (colour)
% D = original depth
% Nz = target no. of rows in z-direction
% This function will take the width stated and scale to find the depth of the image.



%*** MATERIAL PROPERTIES*************************************************

matprop = [
        % unit  conductivity  density  heat capacity  heat production
          1	    2.6           2800	        1347       0.8244*10^-6     % yellow
          2	    2.6           2800	        1347       0.8244*10^-6     % orange
          3	    2.6           2800	        1347       0.8244*10^-6     % red
          4	    1e-6            1   	    1000	    0];                % air/water


% set these material properties across relevant rock type in grid
% reshaping the material properties matrix to fit the grid with relevant unit and number of columns by rows

sigma = reshape(matprop(units,2),Nz,Nx); % conductivity
rho = reshape(matprop(units,3),Nz,Nx); % density
Cp = reshape(matprop(units,4),Nz,Nx); % specific heat capacity
Hr = reshape(matprop(units,5),Nz,Nx); % heat rate
        
k0 = sigma*10^3 ./ rho ./ Cp; % calculate heat diffusivity

% sigma: Themal conductivity [W/(mK)]
% Cp: Specfic heat capacity  [J/(kgK)]
% rho: Density               [kg/m^3]
% Hr: Heat production rate   [W/m^2]
% ko: Diffusivity constant   [m^2/s]



%*** MODEL PARAMETERS*************************************************

dTdz = [0, 35/1000]; % set geothermal gradient to be either 0 or 35degC/km [deg C/m]
T0  = 10; % temperature of surface [deg C]
Tair = 10; % temperature of air [deg C]
nop   = 5000; % output figure produced every 'nop' steps
yr    = 3600*24*365; % seconds in a year [s]
tend  = 1e6*yr; % end time [s]
CFL   = 0.8; % Time step limiter



%***  RUN MODEL*************************************************
run('./geotherm_A.m');