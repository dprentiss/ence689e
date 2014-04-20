% save_time_invar_params

function save_time_invar_params(elev,slope,aspect,lat,elev_nom)


% Model parameter inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal constants

%global cp Lv rhow g sigma vonKarm omega 
cp=1004;                % specific heat capacity of air (J/kg/K)
Lv=2.5e6;               % latent heat of vaporization (J/kg)
rhow=1000;              % density of water (kg/m^3)
g=9.81;                 % accel. of gravity m/s^2)
sigma=5.67e-8;          % Stefan-Boltzmann constant
vonKarm=0.4;            % von Karman constant
omega=1/86400;          % diurnal frequency

save time_invar_params.mat cp Lv rhow g sigma vonKarm omega

% Land surface parameters
%global zr emiss alb lat kBinv z0M zdh z0H
emiss=0.97;             % surface emissivity
alb=0.25;               % albedo
h_veg=1.0;              % surface height (m)
z0M=0.1*h_veg;         % surface momentum roughness
zdh=0.7*h_veg;         % zero-plane displacement height
kBinv=2.3;              % von Karman constant times inverse Stanton number
z0H=z0M/exp(kBinv);     % roughness length (heat) (m)
zr=2.;                  % reference (meaurement) height
elev_nom=1000;          % elevation of met. station

save time_invar_params.mat -append emiss alb h_veg z0M zdh kBinv z0H zr elev_nom

% Vegetation parameters
%global  Vc C_veg Rsmin Rsmax LAI RGL d_vpd a_T_veg b_T_veg
Vc=0.8;                 % vegetation fractional cover
C_veg=1e-3;             % Vegetation heat capacity
Rsmin=100;              % minimum stomatal resistance (s/m)
Rsmax=5000;             % maximum stomatal resistance (s/m)
LAI=3.5;                % leaf area index
RGL=100;                % canopy resistance (PAR) factor parameter
% Canopy resistance parameters
d_vpd=0.000238;         % vapor pressure deficit parameter (1/Pa)
a_T_veg=[-0.41; 1.18];  % Parameters in air temp. stress function
b_T_veg=[282.05;314];

save time_invar_params.mat -append Vc C_veg Rsmin Rsmax LAI RGL d_vpd a_T_veg b_T_veg

% Soil parameters
%global d1 d2 THETAs THETAl THETAfc Wfc Wwlt C_Gsat a_soil b_soil p_soil C_1sat C_2ref W_low K_sat
d1=0.05 ;               % surface layer thickness (m)
d2=0.5 ;                % total soil profile depth (m)
THETAs= 0.434 ;         % porosity
THETAl= 0.05 ;          % wilting point VSM
THETAfc=0.75*THETAs;    % field capacity VSM
Wfc=THETAfc/THETAs;     % field capacity rel. sat.
Wwlt=THETAl/THETAs;     % wilting point rel. sat.

save time_invar_params.mat -append d1 d2 THETAs THETAl THETAfc Wfc Wwlt

% Misc. parameters (see Noilhan and Planton, 1989 MWR paper)
% Silt loam parameters
C_Gsat=4.418e-6;
b_soil=5.30;
C_1sat=0.153;
C_2ref=0.8;
W_low=0.01;
a_soil=.105;
p_soil=6.;              
K_sat=1.66e-6;          % Sat. hydraulic conductivity (m/s)

save time_invar_params.mat -append C_Gsat b_soil C_1sat C_2ref W_low a_soil p_soil K_sat

% Topographic parameters
% lat=34;
% slope=0;
% elev=1500;
% elev_nom=1000;
% aspect=0;

save time_invar_params.mat -append lat slope elev aspect


return