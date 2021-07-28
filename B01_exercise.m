% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------

% ----------------------------------------
% Practice with the examples, by changing source/receiver position
% source type, velocity values


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
%       Acoustic wave equation finite diference simulator
% -------------------------------------------------------------------------
%
% This file is a manual and example script using the acoustic simulator
% The user defines the velocity model, the source parameters and 
% the simulation parameters
% The program displays the "evolving" acoustic pressure field amplitude
% in a movie-like figure (snapshots)
% The user can define an optional set of receivers positions: in this case
% the program shows and outputs the pressure recorded at the receivers
% (seismic traces). The receivers that are out of the model are
% automatically rejected, and so the seismic traces can be less than the
% input receivers: the output structure contains the actual position of the
% 'valid' receivers
%
% -------------------------------------------------------------------------
% 1. Model parameters in structure 'model'
%
% Compulsory parameters
% model.x: vector of x grid coordinates [m], Nx elements
% model.z: vector of z grid coordinates [m], Nz elements
% model.vel: matrix of velocity values [m/s], (Nz,Nx) elements
%
% Optional parameters
% model.recx: vector of x coordinates of receivers [m], Nr elements
% model.recz: vector of z coordinates of receivers [m], Nr elements
% model.dtrec: max time sampling interval for seismic traces [s]
%
% -------------------------------------------------------------------------
% 2. Source parameters in structure 'source'
% source coordinates are rounded to nearest grid point
% all of them can be **vectors**, in order to simulate multiple sources
%
% Compulsory parameters
% source.x:  x coordinate of source [m]
% source.z:  z coordinate of source [m]
% source.f0: central frequency of source Ricker wavelet [Hz]
% source.t0: time of source emission (referred to max peak of Ricker)
% source.type: 1 is Ricker, 2 is sinusoid at frequency source.f0
% source.amp: multiplier of source amplitude
%
% -------------------------------------------------------------------------
% 3. Simulation and graphic parameters in structure 'simul'
%
% simul.timeMax:     max simulation time [s]
% simul.borderAlg:   Absorbing boundaries (Yes:1, No:0)
% simul.printRatio:  pressure map shown every printRatio comput. time steps
% simul.higVal:  colormap between -highVal and + highVal (from 0 to 1)
% simul.lowVal:  values between -lowVal and +lowVal zeroed (from 0 to 1)
% simul.bkgVel:  velocity matrix as  a "shadow" in the images (1:yes, 0:no)
%
% Optional parameters
% simul.cmap:    colormap (default 'gray')
%
% -------------------------------------------------------------------------
% 4. Acoustic simulator program call
%
%  recfield=acu2Dpro(model,source,simul);
%
%  recfield.time: time axis of recorded signal [s], Nt elements
%  recfield.data: matrix of pressure at the receivers, (Nt,Nr1)
%  recfield.recx: vector of x grid coordinates of receivers [m], Nr1 elements
%  recfield.recz: vector of z grid coordinates of receivers [m], Nr1 elements
%
% -------------------------------------------------------------------------
% 5. Plotting seismic trace program call
%
%  seisplot2(recfield.data,recfield.time)
% 
%   [fact]=seisplot2(datain,t,tr,scal,pltflg,scfact,colour,clip)
%  
%   function for plotting seismic traces
%  
%   INPUT
%   datain  - input matrix of seismic traces
%   t       - time axis
%   tr      - trace axis
%   scal    - 1 for global max, 0 for global ave, 2 for trace max
%   pltflg  - 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
%             2 plots wiggle traces only, 3 imagesc gray, 4 pcolor gray
%   scfact  - scaling factor
%   colour  - trace colour, default is black
%   clip    - clipping of amplitudes (if <1); default no clipping
%  
%   OUPTPUT
%   fact    - factor that matrix was scaled by for plotting
%   if you want to plot several matrices using the same scaling factor,
%   capture 'fact' and pass it as input variable 'scal' to future matrices
%   with 'scfact' set to 1
% 
% -------------------------------------------------------------------------
% 6. Essential theory
%
%  Sampling interval for stability (computed automatically)
%  dt=0.95*sqrt(1/2)*min(dx,dz)/vMax
%
%  Courant-Friedrick-Levy (CFL) stability condition (display warning)
%  dx<(vMin/(f0*2*8)))
%  dz<(vMin/(f0*2*8)))
%
%  -----------
%  Ricker wavelet
%
%                  *                        !    ***
%                 * *                       !   * ! *
%    ____________*___*___________ T         !  *  !   *
%    ********   *     *   *******           ! *   !     *
%            ***       ***                  !*    !        * * * *
%             >   TD    <                   ------+---------------- F
%                                                F0
%
%    s(t) = (1-Y2*T*T)*exp(-Y2*T*T/2)   S(f) = 2*F^2/(f0^3*sqrt(pi))*
%    Y2 = 2*pi^2*f0^2                     *exp(-F*F/(f0^2))
%    TD = sqrt(6)/(pi*f0)
%
%  -----------

clear all

% ----------------------------------------
% CONSTANT MODEL AND RECEIVERS ON A CIRCLE
% ----------------------------------------

% ----------------------------------------
% 1. Model parameters

model.x   = 0:1:500;     % horizontal x axis sampling
model.z   = 0:1:500;     % vertical   z axis sampling

% temporary variables to compute size of velocity matrix
Nx = numel(model.x);
Nz = numel(model.z);

% velocity model assignement: constant velocity
model.vel=zeros(Nz,Nx)+1000;      % initialize matrix

% ----------------------------------------
% 2. Source parameters (in the center of the model)

source.x    = [100 ];
source.z    = [200 ]; 
source.f0   = [30 ];  %f of wavelet
source.t0   = [0.04  ]; % time delay
source.amp  = [1.5 ]; 
source.type = [1];    % 1: ricker, 2: sinusoidal  at f0


% optional receivers in (recx, recz)
% the program round their position on the nearest velocity grid
phi    = linspace(0,2*pi,100); %they are on circle
radius = 170;

model.recx  = radius * cos(phi) + model.x(end)/2;
model.recz  = radius * sin(phi) + model.z(end)/2;
model.dtrec = 0.004;

% ----------------------------------------
% 3. Simulation and graphic parameters in structure simul

simul.borderAlg=1;
simul.timeMax=0.5; % speed of simulation

simul.printRatio=5;
simul.higVal=.05; % higher than this will be saturated
simul.lowVal=0.01;% lower than this will be discarded
simul.bkgVel=1;

simul.cmap='pink';   % gray, cool, hot, parula, hsv

% ----------------------------------------
% 4. Program call

recfield=acu2Dpro(model,source,simul);

% Plot receivers traces

figure
scal   = 2;  % 1 for global max, 0 for global ave, 2 for trace max
pltflg = 2;  % 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
             % 2 plot wiggle traces only, 3 imagesc gray, 4 pcolor gray
scfact = 1;  % scaling factor
colour = ''; % trace colour, default is black
clip   = []; % clipping of amplitudes (if <1); default no clipping

seisplot2(recfield.data,recfield.time,[],scal,pltflg,scfact,colour,clip)
xlabel('receiver nr')




