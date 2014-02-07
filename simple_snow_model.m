function simple_snow_model

% This is a spatially distributed version of a simple one-layer
% snow mass/energy balance model.
%
% The key inputs are DEM data (elevation, slope, aspect) which are
% distributed in space (located in the input file DEM_data.mat) and the
% nominal meteorological forcings (located in the input file forcing.mat).  
% Here it is assumed that there is one timeseries of the met. forcings and 
% it is located at the mean elevation of the basin.  The domain size is
% tied directly to the DEM.
% 
% Several of the forcings (air temperature, specific humidity, 
% pressure, incoming solar radiation) are distributed in space using
% simple lapse rate and/or other topographic information.
%
% Currently none of the time-invariant snow parameters are distributed
% in space (i.e. assumed uniform).
%
% Author:   Bart Forman

% Define model simulation parameters
start_DOY=335;  % Starts at midnight (beginning) of this day
end_DOY=1;      % Ends at midnight (beginning) of this day
dt=1.0;         % Timestep, Hours -- CURRENTLY MUST MATCH FORCING

% In case simulations spans into new year
if end_DOY>start_DOY
    Nsteps=(end_DOY-start_DOY+1)*24/dt; % # of time steps
else
    Nsteps=((end_DOY-0)+(365-start_DOY))*24/dt;
end

% Define global variables
global Rd es0 T0 eps Ls elev slope aspect N_pix elev_nom Lv Rv lat_deg lon_deg albedo doy hour lat_vec lon_vec

% Physical constants
Ls=2.83e6;      % Latent heat of sublimation (J/kg)
Lv=2.5e6;       % Latent heat of vaporizaton (J/kg)
Lf=Ls-Lv;       % Latent heat of fusion (J/kg)
rhow=1000;      % Density of water (kg/m^3)
Rd=287;         % Dry air gas constant (J/kg/K)
Rv=461;         % Vapor gas constant (J/kg/K)
T0=273.15;      % Clausius-Clapeyron reference temp. (K)
es0=611;        % Clausius-Clapeyron reference vapor pressure (Pa)
eps=Rd/Rv;      % ratio of dry air to vapor gas constants
vonKarm=0.4;    % vonKarman constant
sigma=5.67e-8;  % Stefan-Boltzmann constant (W/m^2/K^4)
cp=1004;        % specific heat capacity of air (J/kg/K)
T_f=273.15;     % freezing temp.

% Parameters
albedo=0.75;    % snow albedo
emiss=0.91;     % snow thermal emissivity
h=0.03;         % height of roughness elements (bare snow only)
z0m=0.1*h;      % momentum roughness length (meters)
z0h=0.1*z0m;    % scalar roughness length (meters)
Csnow=0.55e+6;  % snow heat capacity (J/m^3/K)
zm=2.;          % Ref.-level (meters)
dispht=0.7*h;  % zero-plane disp. height (m)

% Load time-varying forcing file (this is already setup for this simulation
% period)
load forcing.mat

% Load topographic data; includes elevation, slope and aspect maps
% as well as lat/lon coordinates 
load DEM_data.mat

% Use DEM to setup domain size
N_lat=length(lat_deg);  % Dimension in latitude
N_lon=length(lon_deg);  % Dimension in longitude
N_pix=N_lat*N_lon;      % Number of model pixels

% Reshape topographic arrays into vectors:
elev=reshape(elev,N_pix,1);     % meters
slope=reshape(slope,N_pix,1);   % radians: (horiz. = 0, increases up to pi/2 (90 deg.)
aspect=reshape(aspect,N_pix,1); % radians: south facing =0, north facing: pi/-pi  
% Create coordinates for each pixel (needed in shortwave distribution)
[lon_matrix,lat_matrix]=meshgrid(lon_deg,lat_deg);
lon_vec=reshape(lon_matrix,N_pix,1);
lat_vec=reshape(lat_matrix,N_pix,1);

% Assume nominal forcings are at the mean elevation
elev_nom=mean(elev);

% Initial condition (assumed uniform in space)
SWE(:,1)=50.*ones(N_pix,1);     %  millimeters
Tsnow(:,1)=262.*ones(N_pix,1);   % deg. K


for istep=1:Nsteps-1
    if floor(start_DOY+istep*dt/24) <= 365
        doy=floor(start_DOY+istep*dt/24)
        hour=((start_DOY+istep*dt/24)-doy)*24
    else 
       doy=floor(start_DOY+istep*dt/24)-365
       hour=((start_DOY+istep*dt/24)-(doy+365))*24
    end

    % Load nominal forcing variables
    P_nom=forcing(istep,1)*3600;    % precip. rate (mm/hr)
    Psrf_nom=forcing(istep,2);      % suface pressure (Pa)
    Ta_nom=forcing(istep,3);        % reference-level (2m) air temperature
    qa_nom=forcing(istep,4);        % reference-level (2m) specific humidity
    Rs_down_nom=forcing(istep,5);   % incoming shortwave radiation (W/m^2)
    wind_nom=forcing(istep,6);      % reference-level (2m) wind speed (m/s)
    
    % Distribute forcing in space based on elevation:
    [P,Psrf,Ta,qa,Rs_down,wind,rho]=dist_forcings(P_nom,Psrf_nom,Ta_nom,qa_nom,Rs_down_nom,wind_nom);
    
    % Assume surface humidity is equal to sat. value with respect to ice
    qsurf=eps./Psrf*es0.*exp(Ls/Rv*(1/T0-1./Tsnow(:,istep)));  % kg/kg
    % aerodynamic resistance (s/m) -- Assumes neutral conditions
    ra=(log((zm-dispht)/z0m)*log((zm-dispht)/z0h))/vonKarm^2./wind;
    
    % Evaporation
    E=rho.*(qsurf-qa)./ra;    % evaporation rate (kg/m^2/s)
    LE=Ls*E;                % latent heat of sublimation
    E=E*1000/1000*3600;     % evaporation rate (mm/hr)
                      
    % Net radiation
    % Need to compute incoming longwave if it is not available
    ea=qa.*Psrf/eps/100;         % air vapor pressure (mb)
    emiss_a=0.74+0.0049*ea;     % air emissivity -- Idso model
    Rl_down=emiss_a*sigma.*Ta.^4; % incoming longwave radiation (W/m^2)    
    Rn=Rs_down*(1-albedo) + Rl_down - emiss*sigma*Tsnow(:,istep).^4;
    
    % Sensible heat flux
    H=rho*cp.*(Tsnow(:,istep)-Ta)./ra;
    
    % Classify precip. occuring at air temp. above freezing as rain; assume
    % all of this immediately runs  off (i.e. doesn't contribute to SWE)
    I=find(Ta>T_f & P>0);
    if ~isempty(I) 
        P(I)=0;
    end
    
    % Snow mass/energy balance
    SWE(:,istep+1) = SWE(:,istep) + dt*(P-E);   % Note this is done in mm
    Tsnow(:,istep+1) = Tsnow(:,istep) +dt/Csnow*(Rn-LE-H)*3600; % 3600 converts to fluxes to J/hr/m^2
    
    % Check for phase change
    I=find(Tsnow(:,istep+1)>T_f & SWE(:,istep+1)>0);
    if ~isempty(I)
        % Energy that would have gone into melting
        melt_energy=((Tsnow(I,istep+1)-T_f).*Csnow.*SWE(I,istep+1)/1000);
        % Subtract melt mass from SWE
        SWE(I,istep+1)=SWE(I,istep+1)-(melt_energy/rhow/Lf)*1000;
        Tsnow(I,istep+1)=T_f;
    end
        
end

% Maps of distributed snow states (N_lon, N_lat, time)
SWE_maps=reshape(SWE,N_lon,N_lat,Nsteps); 
Tsnow_maps=reshape(Tsnow,N_lon,N_lat,Nsteps);

save output.mat SWE_maps Tsnow_maps lat_deg lon_deg

return
    
%%
function [P,Psrf,Ta,qa,Rs_down,wind,rho]=dist_forcings(P_nom,Psrf_nom,Ta_nom,qa_nom,Rs_down_nom,wind_nom)
% This function distributes the nominal forcing in space
global Rd es0 T0 eps Ls elev N_pix elev_nom Lv Rv
g=9.81; % accel. of gravity m/s^2

% Adiabatic lapse rate
gamma_d=9.8e-3;

% elevation difference vector
del_z=elev-elev_nom;    
% Distribute temp.
Ta=Ta_nom*ones(N_pix,1)-gamma_d*del_z;

% Distribute humidity, pressure, and density
ea_nom=Psrf_nom/eps*qa_nom;                     % nominal vapor pressure (Pa)
e_s_nom=es0*exp(Lv/Rv*(1/T0-1/Ta_nom));         % nominal sat. vapor pressure
RH_nom=ea_nom/e_s_nom;                          % nominal RH                                   
e_s=es0.*exp(Ls/Rv*(1/T0-1./Ta));               % dist. sat. vapor pressure
ea=e_s*RH_nom;                                  % dist. vapor pressure

rho_nom=Psrf_nom/Rd/Ta_nom/(1+0.608*qa_nom);    % nominal air density (kg/m^3)
% Dist. pressure using nominal density
Psrf=Psrf_nom-rho_nom*g*del_z;
% Specific humidity
qa=eps./Psrf.*ea;
% Air density
rho=Psrf/Rd./Ta./(1+0.608*qa);

% Assume precip. is distributed uniformly (could use lapse rate)
P=P_nom*ones(N_pix,1);
% Assume wind is distributed uniformly
wind=wind_nom*ones(N_pix,1);

%Rs_down=Rs_down_nom*ones(N_pix,1);

[Rs_down]=dist_shortwave_radiation(Rs_down_nom);

return

%%
function [Rs_down]=dist_shortwave_radiation(Rs_down_nom)
% This "disaggregates" a measurement (on a horizontal plane) in space
% based on topographic variables -- Taken from Allen et al., 2006, Agric.
% Forest Meteorology

global slope aspect lat_vec lon_vec albedo doy hour N_pix
Gsc = 1367; % Global solar constant in W/m2

% Convert to radians
lat=lat_vec*pi/180;
lon=lon_vec*pi/180;

% sky-vieww factor (assuming isotropic surface)
f_i=(1+cos(slope))/2;
%   
% correction factor for relative Sun-Earth as a function of the day of
% year.
d = (1.0 + 0.017*cos(2*pi/365*(186-doy)))^2; 
% declination angle in radians
delta0 = 23.45*pi/180; 
delta=delta0*cos(2*pi/365*(172-doy)); 
% Compute sunrise/sunset hour angle
omega_s=acos(-tan(delta)*tan(mean(lat)));
% Compute hour angle
omega=2*pi*(hour-12)/24;  
% Compute cosing of solar zenith angle on a horiz. plane
cos_theta_hor=sin(delta)*sin(lat) + cos(delta)*cos(lat)*cos(omega);
I_hor=find(asin(cos_theta_hor)<0);
cos_theta_hor(I_hor)=0; % Set to zero in cases where it is negative
% Compute radiation on horizontal plane (this is assuming transmissivity
% =1.0)
Ra_hor=Gsc*cos_theta_hor/d;
% Compute cosine of solar zenith angle adjusted for slope/aspect
if omega > -omega_s && omega < omega_s  % Daytime
    % cos of solar zenith angle taking into account slope/aspect
    cos_theta =sin(delta).*sin(lat).*cos(slope) - ...
                    sin(delta).*cos(lat).*sin(slope).*cos(aspect) + ...
                    cos(delta).*cos(lat).*cos(slope).*cos(omega) + ...
                    cos(delta).*sin(lat).*sin(slope).*cos(aspect).*cos(omega) + ...
                    cos(delta).*sin(aspect).*sin(slope).*sin(omega);
else  % Nighttime
    cos_theta=zeros(length(lat),1);
end
I=find(asin(cos_theta)<0);
cos_theta(I)=0;
% Ratio of actual cos(theta) to that on horiz. plane
f_B=cos_theta./cos_theta_hor;
I=find(cos_theta_hor)==0;
if ~isempty(I)
    f_B(I)=0;
end
% Estimate atmospheric transmissivity based on upper limit and measured
% rad.
if mean(Ra_hor)==0
    tau_SW_hor=0 ;
else
    tau_SW_hor=Rs_down_nom/mean(Ra_hor);
end
% Empirical function for transmissivity index for actual direct beam
% radiation on a horizontal surface
if tau_SW_hor < 0.175
    K_B_hor=0.016*tau_SW_hor;
elseif tau_SW_hor < 0.42
    K_B_hor=0.022-0.280*tau_SW_hor+0.828*tau_SW_hor^2+0.765*tau_SW_hor^3;
else
    K_B_hor=1.56*tau_SW_hor-0.55;
end
% Index for actual diffuse radiation on horizontal surface
K_D_hor=tau_SW_hor-K_B_hor;

% Compute actual incoming radiation
if tau_SW_hor==0
    Rs_down=zeros(N_pix,1);
else
    Rs_down=Rs_down_nom*(f_B*K_B_hor/tau_SW_hor + f_i*K_D_hor/tau_SW_hor + albedo*(1-f_i));
end

return
