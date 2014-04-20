% FRmodel.m

function FRmodel

% This is a Newton-Raphson solver for the Force-Restore land surface model

global rhoair
global cp Lv rhow g sigma vonKarm omega 
global zr emiss alb lat kBinv z0M zdh z0H
global  Vc C_veg Rsmin Rsmax LAI RGL d_vpd a_T_veg b_T_veg
global d1 d2 THETAs THETAl THETAfc Wfc Wwlt C_Gsat a_soil b_soil p_soil C_1sat C_2ref W_low K_sat

% Load control_parameters
load control_params.mat

% Load time-invariant parameters
load time_invar_params.mat

% Initial conditions
Ts(1)=Ts0; Td(1)=Td0; W1(1)=W10; W2(1)=W20;

global Tair qair Psurf Rsol PPT U_r LWdown

% Load forcing over simulation period
load forcing.mat

% Start numerical solution scheme
dt=0.5*(dt_max+dt_min);

t_beg=Day_beg;
t_end=Day_end;
time_index=1;
t=t_beg;
time(time_index)=t_beg;

global H LE Eg Etr
Ta_nom=interp1(forcing(:,1),forcing(:,2),t_beg);
Psrf_nom=interp1(forcing(:,1),forcing(:,3),t_beg);
qa_nom=interp1(forcing(:,1),forcing(:,4),t_beg);
Rs_down_nom=interp1(forcing(:,1),forcing(:,5),t_beg);
P_nom=interp1(forcing(:,1),forcing(:,6),t_beg);
wind_nom=interp1(forcing(:,1),forcing(:,7),t_beg);
LW_down_nom=interp1(forcing(:,1),forcing(:,8),t_beg);
%rhoair=Psurf/287./Tair/(1+0.61*qair);

DOY=floor(time(time_index));
HOUR=24*(time(time_index)-DOY);

% DISTRIBUTE FORCINGS -- based on topographic variability
[PPT,Psurf,Tair,qair,Rsol,U_r,rhoair,LWdown]= ...
            dist_forcings(P_nom,Psrf_nom,Ta_nom,qa_nom,Rs_down_nom,wind_nom, LW_down_nom, ...
            lat, elev, elev_nom, aspect, slope, alb, DOY, HOUR);

[ybeg]=initialize(Ts,Td,W1,W2);

Rnet(time_index)=Rsol(time_index)*(1-alb)+LWdown(time_index)-emiss*sigma*Ts^4;
[H,LE,Eg,Etr]=turb_flux_calc(Ts,W1,W2);

SH(time_index)=H;
LH(time_index)=LE;

% Main time-stepping loop
while t < Day_end

        Ta_nom=interp1(forcing(:,1),forcing(:,2),t+dt/86400);
        Psrf_nom=interp1(forcing(:,1),forcing(:,3),t+dt/86400);
        qa_nom=interp1(forcing(:,1),forcing(:,4),t+dt/86400);
        Rs_down_nom=interp1(forcing(:,1),forcing(:,5),t+dt/86400);
        P_nom=interp1(forcing(:,1),forcing(:,6),t+dt/86400);
        wind_nom=interp1(forcing(:,1),forcing(:,7),t+dt/86400);
        LW_down_nom=interp1(forcing(:,1),forcing(:,8),t+dt/86400);
        
        DOY=floor(time(time_index));
        HOUR=24*(time(time_index)-DOY);

        % DISTRIBUTE FORCINGS -- based on topographic variability
        [PPT,Psurf,Tair,qair,Rsol,U_r,rhoair,LWdown]= ...
            dist_forcings(P_nom,Psrf_nom,Ta_nom,qa_nom,Rs_down_nom,wind_nom, LW_down_nom, ...
            lat, elev, elev_nom, aspect, slope, alb, DOY, HOUR);
        
        % Reset beginning guess for state variable to previous time step converged value   
        if time(time_index) > time(1)
            ybeg=yend;
        end

        % Start Newton-Raphson iteration for given time step
        for m=1:Maxit 
      
            % Set intermediate solution ystar to beginning guess
            if m==1
                ystar=ybeg;
            end
      
            % Evaluate state tendency equations
            [FF]=prog_eval(ystar);
      
            % Evaluate Jacobian of state tendency equations
            [dFdy]=dFdyeval(ystar); 
      
            II = eye(4);
            AA= II -dt*dFdy;
            bb = dt*(FF - dFdy*ystar(1:4)) + ybeg(1:4);
      
            % Compute solution for linearized equation
            Ainv=inv(AA);
            ysol=Ainv*bb;
      
            if ysol(3)>=1.
                ysol(3)=1.;
    	    end
    	    if ysol(3) <= Wwlt
                ysol(3)=Wwlt;
            end        
            if ysol(4)>=1.
     		    ysol(4)=1.;
    	    end
    	    if ysol(4) <= Wwlt
                ysol(4)=Wwlt;
            end  
              
            DIFF=abs(ystar(1:4)-ysol);
      
            if DIFF(1)<DIFFtol(1) && DIFF(2)<DIFFtol(1) && DIFF(3)<DIFFtol(2) && DIFF(3)<DIFFtol(2)
                break 	% exit iterative solver loop
            end   
  
            ystar(1:4)=ysol;
      
        end % End iterative solver loop

        if m==Maxit
            dt=max(dt/2,dt_min);
        else            
            t = t + dt/86400. ;
        
            
            if m < 3 
                dt = min(2*dt,dt_max);
                time_index=time_index+1;
            elseif m >= 7 
                dt = max(0.5*dt,dt_min);
                time_index=time_index+1;
            else
                time_index=time_index+1;
            end
       
            yend=ysol(1:4);

            time(time_index)=t;
            Ts(time_index)=yend(1);
            Td(time_index)=yend(2);
            W1(time_index)=yend(3);
            W2(time_index)=yend(4);
           
            Rnet(time_index)=Rsol*(1-alb)+LWdown-emiss*sigma*Ts(time_index)^4;
            [H,LE,Eg,Etr]=turb_flux_calc(Ts(time_index),W1(time_index),W2(time_index));
            SH(time_index)=H;
            LH(time_index)=LE;

        end
        
end  % End time stepping loop
% Create output
t_reg=[t_beg:dt_output/86400:t_end]';

states=[time' Ts' Td' W1' W2'];
states=interp1(states(:,1),states,t_reg);
save states.out states -ascii

fluxes=[time' Rnet' SH' LH'];
fluxes=interp1(fluxes(:,1),fluxes,t_reg);
save fluxes.out fluxes -ascii

%toc
return

%%%%%%% MODEL FUNCTIONS AND SUBROUTINES SHOWN BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ybeg]=initialize(Ts,Td,W1,W2)
    global THETAl THETAs
    
    smin=THETAl/THETAs;
    ybeg(1)=Ts(1); ybeg(2)=Td(1); ybeg(3)=W1(1); ybeg(4)=W2(1);
    if ybeg(3)>=1.
        ybeg(3)=1.;
    end
    if ybeg(3) <= smin
        ybeg(3)=smin;
    end
    if ybeg(4)>=1.
        ybeg(4)=1.;
    end
    if ybeg(4) <= smin
        ybeg(4)=smin;
    end
    
    ybeg=ybeg';
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FF]=prog_eval(y)
    global rhow sigma omega 
    global alb emiss Vc C_veg
    global d1 d2 THETAs C_Gsat a_soil b_soil p_soil C_1sat C_2ref W_low K_sat
    global Rsol PPT LWdown
    
    global C_T C_G Rn c1 c2 W_eq
    
    Ts=y(1); Td=y(2); W1=y(3); W2=y(4);
   
    % Compute Ts tendency term
    C_G=C_Gsat*(W2)^(-b_soil/2/log(10));
    
    C_T=((1-Vc)/C_G + Vc/C_veg)^(-1);
    
    Rn=Rsol*(1-alb)+LWdown-emiss*sigma*Ts^4;   
 
    [H,LE,Eg,Etr]=turb_flux_calc(Ts,W1,W2);
    
    FF(1) = C_T*(Rn - H -LE) - 2*pi*omega*(Ts-Td) ;
   
    % Compute Td tendency
    FF(2) = omega*(Ts-Td) ;
   
   % Compute W1 tendency
   c1=C_1sat*W1^(-(b_soil/2+1));
   c2=C_2ref*(W2/(1-W2+W_low));
   
   W_eq=W2-a_soil*W2^p_soil*(1-W2^(8*p_soil));
   
   FF(3) = c1/rhow/d1/THETAs*(PPT-Eg) - c2*omega*(W1 - W_eq);
   
   % Compute W2 tendency
   qflux=K_sat*W2^(2*b_soil+3);
   FF(4) = 1./rhow/d2/THETAs*(PPT - Eg - Etr); - qflux/d2/THETAs ;
      
   FF=FF';
   
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dFdy]=dFdyeval(y)
    global rhoair cp Lv rhow g sigma omega 
    global zr  emiss  z0H Vc 
    global d1 d2 THETAs Wfc Wwlt a_soil b_soil p_soil C_2ref W_low K_sat 
    global Tair qair Psurf PPT U_r
    
    global C_T C_G Rn c1 c2 W_eq
    global C_B H LE Eg h_u h_v R_a  F2 F3 R_s
    global RiB C_BN f_RiB b_stab c_stab d_stab
    
	Ts=y(1); W1=y(3); W2=y(4);
  
%%% Compute dFTsdy derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dFTs/dTs
    dRndTs=4*emiss*sigma*Ts^3;
    %dRiBdTs=-g*zr/Tair/U_r^2;
    dRiBdTs=-g*zr/U_r^2*(Tair/Ts^2);
    %df_RiBdTs=-20*exp(10*RiB)*dRiBdTs;
    % Louis formulation
    if RiB > 0
        dfdRiB=-f_RiB^2*3*b_stab*sqrt(1+d_stab*RiB)*(1 + d_stab/2*RiB/(1+d_stab*RiB) );
    else
        dfdRiB=-(3*b_stab*(1+3*b_stab*C_BN*c_stab*sqrt(zr/z0H*abs(RiB))) - ...
                             3*b_stab*RiB*(3/2*b_stab*C_BN*c_stab/sqrt(zr/z0H*abs(RiB))*(-zr/z0H)))/...
                                   (1+3*b_stab*C_BN*c_stab*sqrt(zr/z0H*abs(RiB)))^2;
    end        
    df_RiBdTs=dfdRiB*dRiBdTs;
    dCBdTs=C_BN*df_RiBdTs;
    dHdTs=rhoair*cp*U_r*(dCBdTs*(Ts-Tair) + C_B);
    [e_s]=e_sat_calc(Ts);
    q_sat=0.622/Psurf*e_s;
    dqsatdTs=0.622/Psurf*Lv/461*e_s/Ts^2;
    dEgdTs=(1-Vc)*rhoair*U_r*(dCBdTs*h_u*(q_sat-qair) + C_B*h_u*dqsatdTs  );
    dRadTs=-U_r*R_a^2*dCBdTs;
    %dF3dTs=-d_vpd*Lv/461*e_s/Ts^2;
    dF3dTs=0;
    dRsdTs=-R_s*(1/F3*dF3dTs);
    dhvdTs=(dRadTs*R_s-R_a*dRsdTs)/(R_a+R_s)^2;
    dEtrdTs=Vc*rhoair*U_r*(dCBdTs*(h_v*(q_sat-qair)) + C_B*(dhvdTs*(q_sat-qair) + h_v*dqsatdTs) );
    dLEdTs=Lv*(dEgdTs+dEtrdTs);
    dFdy(1,1)=C_T*(dRndTs - dHdTs - dLEdTs) - 2*pi*omega; 
    % dFTs/dTd
    dFdy(1,2)=2*pi*omega;
    % dFTs/dW1
    dhudW1=1/2*(1-cos(W1*pi/Wfc))*sin(W1*pi/Wfc)*pi/Wfc;
    dEgdW1=(1-Vc)*rhoair*C_B*U_r*q_sat*dhudW1;
    dLEdW1=Lv*(dEgdW1);
    dFdy(1,3)=C_T*(-dLEdW1);
    % dFTs/dW2
    %dF2dW2=5/(Wfc-Wwlt)*F2^2*(exp(-5*(w_af-b_w)));
    if W2 > Wfc || W2 < Wwlt
        dF2dW2=0;
    else
        dF2dW2=1./(Wfc-Wwlt);
    end
    dRsdW2=-R_s/F2*dF2dW2;
    dhvdW2=-h_v/(R_a+R_s)*dRsdW2;
    dEtrdW2=Vc*rhoair*C_B*U_r*q_sat*dhvdW2;
    dLEdW2=Lv*dEtrdW2;
    dCGdW2=-C_G*b_soil/2/log(10)/W2;
    dCTdW2=C_T^2*(1-Vc)/C_G^2*dCGdW2;
    dFdy(1,4)=dCTdW2*(Rn-H-LE) + C_T*(-dLEdW2);
%%%% Compute dFTddy derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dFTd/dTs
    dFdy(2,1)= omega;
    % dFTd/dTd
    dFdy(2,2)=-omega;
    % dFTd/dW1
    dFdy(2,3)=0;
    % dFTd/dW2
    dFdy(2,4)=0;
%%%% Compute dFW1dy derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dFW1/dTs
    dFdy(3,1)=1/rhow/d1/THETAs*(c1*(-dEgdTs));  
    % dFW1/dTd
    dFdy(3,2)=0;
    % dFW1/dW1
    dc1dW1=-(b_soil/2 +1)*c1/W1;
    dFdy(3,3)=1/rhow/d1/THETAs*( dc1dW1*(PPT-Eg) + c1*(-dEgdW1) ) - (omega*c2);
    % dFW1/dW2
    dc2dW2=C_2ref*(1+W_low)/(1+W2+W_low)^2;
    dWeqdW2=1 - a_soil*p_soil*W2^(p_soil-1)*( 1-9*W2^(8*p_soil) ); 
    dFdy(3,4)= -omega*(dc2dW2*(W1-W_eq) + c2*(-dWeqdW2) );
%%%% Compute dFW2dy derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dFW2/dTs
    dFdy(4,1)= 1/rhow/d2/THETAs*(-dEgdTs - dEtrdTs);
    % dFW2/dTs
    dFdy(4,2)=0;
    % dFW2/dW1
    dFdy(4,3)= 1/rhow/d2/THETAs*(-dEgdW1);    
    % dFW2/dW2
    dqdW2=(2*b_soil+3)*K_sat*W2^(2*b_soil+2);
    dFdy(4,4)= 1/rhow/d2/THETAs*(- dEtrdW2) -1/d2/THETAs*(dqdW2);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CB_H]=drag_coeff(Ta,Ts,Ur)
    global vonKarm zr z0H g zdh z0M

    global RiB C_BN f_RiB b_stab c_stab d_stab
    
    %C_BN=vonKarm^2/((log(zr/z0H))^2-kBinv*log(zr/z0H))
    C_BN=vonKarm^2/log((zr-zdh)/z0M)/log((zr-zdh)/z0H);
    %RiB= g*zr*(Ta-Ts)/Ta/Ur^2;
    RiB= g*zr*(Ta-Ts)/Ts/Ur^2;
    % Capparini et al. (2003) formulation
    %f_RiB=1+2*(1-exp(10*RiB));
    % Louis et al. (1982) formulation
    b_stab=5; c_stab=5; d_stab=5;
    if RiB > 0      % stable
        f_RiB=1/(1+3*b_stab*RiB*sqrt(1+d_stab*RiB) );
    else            % unstable
        f_RiB= 1 - 3*b_stab*RiB/(1+3*b_stab*C_BN*c_stab*sqrt(zr/z0H*abs(RiB)) );
    end
    
    CB_H=C_BN * f_RiB;
    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e_sat]=e_sat_calc(T)

    Lv=2.5e6; Rv=467; es0=611; T0=273.15;
    
    e_sat=es0*exp(Lv/Rv*(1/T0 - 1./T));
    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,LE,Eg,Etr]=turb_flux_calc(Ts,W1,W2)

    global rhoair cp Lv 
    global Wfc Wwlt Vc Rsmin Rsmax LAI RGL d_vpd a_T_veg b_T_veg
    global Tair Psurf qair Rsol U_r

    global C_B h_u h_v R_a F1 F2 F3 F4 R_s
    
    [C_B]=drag_coeff(Tair,Ts,U_r);
    H=rhoair*cp*C_B*U_r*(Ts-Tair);
    
    if (W1<Wfc)
        h_u=0.25*(1-cos(W1/Wfc*pi))^2; % Taken from Xiu and Pleim, JAM, 2001
    else
        h_u=1;
    end
    [e_s]=e_sat_calc(Ts);
    q_sat=0.622/Psurf*e_s;
    Eg=(1-Vc)*rhoair*C_B*U_r*h_u*(q_sat-qair);

    % Aerodynamic resistance
    R_a=1/C_B/U_r;
    
    % Compute canopy resistance
    % PAR factor
    ff=0.55*2*Rsol/RGL/LAI;
    F1=(1+ff)/(ff+Rsmin/Rsmax);
    % Soil moisture factor
    if W2 > Wfc
        F2=1;
    elseif W2 < Wwlt
        F2=0;
    else
        F2=(W2-Wwlt)/(Wfc-Wwlt);
    end
    % VPD factor
    eair=Psurf/0.622*qair;
    F3=1-d_vpd*(e_s-eair);
    % If modifying F3, make sure to make corresponding changes in dFdyeval
    % to derivatives
    if F3 < 0.25    % As suggested by Pleim and Xiu, JAM, 2001
        F3 = 0.25;
    end
    if F3 > 1
        F3=1;
    end
    % Air temp. factor
    if Tair <302.15
        a_T=a_T_veg(1);
        b_T=b_T_veg(1);
    else
        a_T=a_T_veg(2);
        b_T=b_T_veg(2);
    end
    F4=1./(1 + exp(a_T*(Tair-b_T)) );
    % Canopy resistance
    R_s=Rsmin/LAI*F1/F2/F3/F4;
    
    % h_v ("rel. humidity" factor)
    h_v=R_a/(R_a + R_s);
    
    % Transpiration
    Etr=Vc*rhoair*C_B*U_r*h_v*(q_sat-qair);
    LE=Lv*(Eg+Etr);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,Psrf,Ta,qa,Rs_down,wind,rho,LW_down]= ...
    dist_forcings(P_nom,Psrf_nom,Ta_nom,qa_nom,Rs_down_nom,wind_nom, LW_down_nom, ...
    lat, elev, elev_nom, aspect, slope, albedo, DOY, HOUR)
% This function distributes the nominal forcing in space based on
% topographic characteristics

% Physical constants
g=9.81; % accel. of gravity m/s^2
Rd=287;
es0=611;
T0=273.15;
eps=0.622;
Lv=2.5e6;
Rv=461;
sigma=5.67e-8;

N_pix=1;

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
e_s=es0.*exp(Lv/Rv*(1/T0-1./Ta));               % dist. sat. vapor pressure
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

% Idso formula
emiss_a_nom=0.74+0.0049*ea_nom/100;
LW_nom_pred=emiss_a_nom*sigma*Ta_nom^4;
emiss_a=0.74+0.0049*ea/100;
LW_pred=emiss_a*sigma*Ta^4;
factor=LW_pred/LW_nom_pred;
LW_down=LW_down_nom*factor;

[Rs_down]=dist_shortwave_radiation(Rs_down_nom,slope,aspect,lat,albedo,DOY,HOUR);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rs_down]=dist_shortwave_radiation(Rs_down_nom,slope,aspect,lat, albedo, doy, hour)
% This "disaggregates" a measurement (on a horizontal plane) in space
% based on topographic variables -- Taken from Allen et al., 2006, Agric.
% Forest Meteorology

Gsc = 1367; % Global solar constant in W/m2

% Convert to radians
lat=lat*pi/180;
N_pix=1;

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
if Ra_hor==0
    tau_SW_hor=0;
else
    tau_SW_hor=Rs_down_nom/Ra_hor;
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
elseif tau_SW_hor >= 1.0 
    Rs_down=Rs_down_nom;
else
    Rs_down=Rs_down_nom*(f_B*K_B_hor/tau_SW_hor + f_i*K_D_hor/tau_SW_hor + albedo*(1-f_i));
end

return
