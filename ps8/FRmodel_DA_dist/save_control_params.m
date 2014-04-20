% save_control_params script

function save_control_params(Ts0,Td0,W10,W20,Day_beg,Day_end)

%Day_beg=160.;
%Day_end=161.;
% Numerial solution parameters
dt_forcing=1800.;            % Time step of forcing file (s)
dt_max= dt_forcing/2 ;      % Max. time step used in model (s)
dt_min=5*60 ;                 % Min. time step used in model (s)
dt_output=900.;              % Time step of output (s)
% Nonlinear iteration parameters
Maxit=10;                   % Max. number of iterations in time step
Tdifftol=1e-3;              % Temp. difference tolerance (K)
Sdifftol=1e-4;              % Soil moisture difference tolerance (-)
DIFFtol=[Tdifftol;Sdifftol];

% Initial Conditions:
%Ts0=295. ; Td0=296. ; W10= 0.5 ; W20=0.5;

save control_params.mat Day_beg Day_end dt_forcing dt_max dt_min dt_output Maxit DIFFtol Ts0 Td0 W10 W20

return
