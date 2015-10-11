%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lab Template for Project 1 TSFS09     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Revision: 1.0 $
% $Date: 2014/06/03 $

%  Decide here if the script generates the validation plots or not by
%  changing the binary varable doPlot
clear all
doPlot = 1;                                     % [-] doPlot==1 generate validation plots, doPlot==0 the contrary.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load the engine data from Project 1a %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the measurement files to the path with: File - Set Path
load('EngineMapTSFS09.mat');
load('AirStep.mat'); % Load here your dynamic measurements
load('FuelStep.mat'); %Load here your dynamic measurements
length=306;

%%
%Extra def
T_t=300;
rTurb=0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Extract Data From Engine Map     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_af            = EngineMap.p_af;               % [kPa]         Pressure air filter
p_bef_thr       = EngineMap.p_bef_thr;                        % [kPa]         Pressure before throttle
p_im            = EngineMap.p_im;%???                        % [kPa]         Pressure intake manifold
p_em            = EngineMap.p_em;%???;                       % [kPa]         Pressure exhaust manifold
p_es            = EngineMap.p_es;%???                        % [kPa]         Pressure after turbine
T_af            = EngineMap.T_af;%???                        % [K]           Temperature air filter 
T_bef_thr       = EngineMap.T_bef_thr;%???                        % [K]           Temperature before throttle
T_im            = EngineMap.T_im;%???                        % [K]           Temperature intake manifold
T_em            = EngineMap.T_em;%???                        % [K]           Temperature exhaust manifold
T_es            = EngineMap.T_es;%???                        % [K]           Temperature after turbine
N               = EngineMap.N;%???                        % [rps]         Engine speed
M_b             = EngineMap.M_e;%???                        % [N*m]         Engine Torque
alpha           = EngineMap.alpha;%???;                       % [-]           Throttle angle
m_dot_at        = EngineMap.m_dot_at;%???                        % [kg/s]        Mass flow after throttle
t_inj           = EngineMap.t_inj;%???                        % [s]           Fuel injection time
lambda_bc_cont  = EngineMap.lambda;%???                        % [-]           Continuous lambda signal before catalyst
lambda_bc_disc  = 0;%???                        % [-]           Discontinuous lambda signal before catalyst
lambda_ac_disc  = 0;%???                        % [-]           Discontinuous lambda signal after catalyst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Set up the engine data     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l               = EngineMap.engine.l;           % [m]           Connecting rod length
a               = EngineMap.engine.a;           % [m]           Crank radius 
B               = EngineMap.engine.B;           % [m]           Bore 
n_cyl           = EngineMap.engine.n_cyl;       % [-]           Number of Cylinders
V_d             = 2*a*pi*B^2/4;                 % [m^3]         Displaced Volume (per cylinder)
n_r             = EngineMap.engine.n_r;         % [-]           Crank revolutions per cycle
V_im            = EngineMap.engine.V_im;        % [m^3]         Intake manifold volume
V_em            = EngineMap.engine.V_em;        % [m^3]         Exhaust manifold volume
V_es            = EngineMap.engine.V_es;        % [m^3]         Exhaust system vomume
V_D             = EngineMap.engine.V_D;         % [m^3]         Total displaced volume
r_c             = EngineMap.engine.r_c;         % [-]           Compression ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Set up the vehicle data    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VehicleMass     = 1520;                         % [kg]          Vehicle mass
WheelRadius     = 0.3;                          % [m]           Wheel radius
VehicleArea     = 2.0;                          % [m^2]         Vehicle frontal area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Set upp other constants     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_exh           = EngineMap.airfuel.R_exh;      % [J/(kg*K)]    Exhaust gas constant
R_air           = EngineMap.airfuel.R_air;      % [J/(kg*K)]    Air gas constant
gamma_air       = EngineMap.airfuel.gamma_air;  % [J/(kg*K)]    Air ratio of specific heats
gamma_exh       = EngineMap.airfuel.gamma_exh;  % [J/(kg*K)]    Exhaust gas ratio of specific heats
cv_air          = R_air/(gamma_air-1);          % [J/(kg*K)]    Specific heat at constant volume, air @T=[-50..40]C
cv_exh          = R_exh/(gamma_exh-1);          % [J/(kg*K)]    Specific heat at constant volume, exhaust gas 
cp_air          = gamma_air*cv_air;             % [J/(kg*K)]    Specific heat at constant pressure, air @T=[-50..40]C
cp_exh          = gamma_exh*cv_exh;             % [J/(kg*K)]    Specific heat at constant pressure, exhaust gas

p_amb           = max(EngineMap.p_af);          % [kPa]         Ambient pressure
T_amb           = mean(EngineMap.T_af);         % [K]           Ambient temperature, mean of air filter temperature

q_LHV           = EngineMap.airfuel.q_LHV;      % [MJ/kg]       Fuel lower heating value
AFs             = EngineMap.airfuel.AFs;        % [-]           Stochiometric A/F for isooctane


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Set upp turbocharger constants     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_tc            = 2.55e-5;                      % [kg*m^2]      Turbocharger Inertia, (From Westin 2002 for MHI TD04HL-15T)

V_ic            = 10e-3;                        % [m^3]         Volume of the intercooler, and pipes between the compressor and the intake manifold
c_tc_fric       = 1e-6;                         % [N*m/(rad/s)] Turbo shaft friction constant
Cd_wg           = 0.8;                          % [-]           Assumed value for WG discharge constant
A_max_wg        = 0.035^2/4*pi;                 % [m^2]         Measured approximation of maximum opening area of WG valve
dp_thrREF       = 10e3;                   % [kPa]         Default desired pressure loss over the throttle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Set up the driver model   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KpDriver        = 0.04;                         % [-]           Driver model Kp value
KiDriver        = 0.02;                         % [-]           Driver model Ki value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Initialize Boost Control   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Boost_control   = 0;                            % [-]           0 -> Boost control disabled, 1st project (binary variable)

%The following constants are required to avoid conflict with the
%turbocharger subsystem, the values will be modified for project 2
KpThr           = 1e-6;                         % [-]           Throttle controller feedback setup
TiThr           = 0.1;                          % [-]           Throttle controller feedback setup

KpWg            = 1e-6;                         % [-]           Wastegate controller setup
TiWg            = 4;                            % [-]           Wastegate controller setup

tau_wg          = 0.1;                          % [s]           Response time of wastegate dynamics        
T_ic            = mean(T_bef_thr);              % [K]           Temperature in intercooler control volume (Assume isoterm model)
dC2             = 0;%???                        % [m]           Outer compressor impeller diameter, measured by the students.
dT1             = 0;%???                        % [m]           Turbine inlet diameter, measured by the students.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Computations for the Gas Pedal     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First plot the measured data to determine the time constant by visual
% inspection.

if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    figure(1),clf
    plot(AirStep.t,AirStep.alpha,'r')
    hold on;
    plot(AirStep.t,AirStep.alpha_ref,'b')
    axis([3 4 19 31])
    grid on
    xlabel('Time [s]')
    ylabel('Throttle position [%]')
    legend('Measured','Reference')    
end
% tau_th is determined to be 0.032 s with the given data (time to reach the
% 63% of the final step value)
tau_th      = 0.032;                             %[s]            Determined gas pedal time constant

% print the value
disp(' ')
disp('Gas pedal:')
disp(['tau_th = ' num2str(tau_th) ' [s]'])

% Proceed to validate the model using the simulink model comparation
% template (MeasModelCompare_throttle.slx)
sim('MeasModelCompare_throttle.slx')            % Run the simulink model

if doPlot
    figure(2)
    plot(AirStep.t,AirStep.alpha_ref,'b')
    hold on;
    plot(AirStep.t,AirStep.alpha,'r')
    plot(GasPedalValidation_Alpha.time,GasPedalValidation_Alpha.signals.values,'k')
    axis([3 4 19 31])
    xlabel('Time [s]')
    ylabel('Throttle position [%]')
    legend('Reference','Measured','Model')    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computations for the mass air flow at throttle %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PI = p_im./p_bef_thr;

est_points = [];
for i=1:length
    if PI(i)<0.73
        est_points = [est_points; i];
    end
end

PI_lim = max(PI(est_points), (2/(gamma_air+1))^(gamma_air/(gamma_air-1)));
Psi = sqrt(2*gamma_air/(gamma_air-1)*(PI_lim.^(2/gamma_air)-PI_lim.^((gamma_air+1)/gamma_air)));

%%%T_bef_thr%%%
T_hat_bef_thr = ones(length,1)\T_bef_thr;
%%%%%%%%%%%%%%%


b=m_dot_at(est_points).*sqrt(R_air*T_hat_bef_thr)./(p_bef_thr(est_points).*Psi);
A_eff=b;
A=[ones(numel(est_points),1), alpha(est_points), alpha(est_points).^2];
x = A\b;
a_0=x(1);
a_1=x(2);
a_2=x(3);
 
A_hat_eff = a_0 + a_1*alpha(est_points) + a_2*alpha(est_points).^2;
 
%m_hat_dot_at = p_bef_thr(est_points).*A_hat_eff.*Psi./sqrt(R_air*T_bef_thr(est_points));
m_hat_dot_at = p_bef_thr(est_points).*A_hat_eff.*Psi./sqrt(R_air*T_hat_bef_thr);


if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    figure(3); clf; hold on
    plot(alpha(est_points),m_dot_at(est_points),'ro')
    [ds,do] = sort(alpha(est_points));
    plot(ds,m_hat_dot_at(do),'k-')
    grid on
    xlabel('Throttle angle [rad]')
    ylabel('Mass air flow at throttle [kg/s]')
    legend('Measured','Model')
    title('When \Pi < 0.73')
    
    rel_error=100*abs(m_dot_at(est_points)-m_hat_dot_at)./m_dot_at(est_points);
    index = [1:numel(rel_error)];
    
    figure(4); clf; hold on
    %plot(index, rel_error, 'r*')
    plot(alpha(est_points), rel_error, 'r*')
    title('Mass air flow: relative error')
    xlabel('Throttle angle [%]')
    ylabel('Relative error [%]')
    
    figure(31); clf; hold on
    plot(alpha(est_points),A_eff,'ro')
    [ds,do] = sort(alpha(est_points));
    plot(ds,A_hat_eff(do),'k-')
    grid on
    xlabel('Throttle angle [rad]')
    ylabel('Effective area [m^2]')
    legend('Measured','Model')
    title('When \Pi < 0.73')
    
    rel_error=100*abs(A_eff-A_hat_eff)./A_eff;
    
    figure(41); clf; hold on
    %plot(index, rel_error, 'r*')
    plot(alpha(est_points), rel_error, 'r*')
    title('Effective area: relative error')
    xlabel('Throttle angle [%]')
    ylabel('Relative error [%]')
    
end

% print the value
disp(' ')
disp('Throttle:')
disp(['a_0 = ' num2str(a_0) ' [m^2]'])
disp(['a_1 = ' num2str(a_1) ' [m^2]'])
disp(['a_2 = ' num2str(a_2) ' [m^2]'])
disp(['T_bef_thr = ' num2str(T_hat_bef_thr) ' [K]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computations for the intake manifold  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=T_im;
A= ones(length,1);
x=A\b;
T_hat_im = x(1);

if 0  %Here doPlot is used, avoids the plot if it is set to 0

    rel_error=100*abs(T_im-T_hat_im)./T_im;
    index = [1:numel(rel_error)];
    
    figure(5); clf; hold on
    plot(T_im, rel_error, 'r*')
    title('Temperature, intake manifold: relative error')
    xlabel('Temperature')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computations for the volumetric efficency %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_vol = m_dot_at*n_r*R_air.*T_im./(p_im*V_d*n_cyl.*N);

b=n_vol;
A=[ones(length,1), sqrt(p_im), N]; 
x=A\b;

n_0=x(1);
n_1 = x(2);
n_2 = x(3);

n_hat_vol = n_0 + sqrt(p_im)*n_1 + n_2*N;

if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    nvol_plot(p_im, N, n_vol, 3);
    title('Volumetric efficiency - measured')
    nvol_plot(p_im, N, n_hat_vol, 4);
    title('Volumetric efficiency - model')
    %%%%%%%% Relative error:%%%%%%%%%%%%%%%
    nvol_plot(p_im, N, 100*abs((n_vol-n_hat_vol)./n_vol), 6);
    title('Volumetric efficiency - relative error')
    zlabel('\eta_{vol}: relative error [%]');
    
end
disp(' ')
disp('Intake manifold:')
disp(['n_0 = ' num2str(n_0) ' [1]'])
disp(['n_1 = ' num2str(n_1) ' [1]'])
disp(['n_2 = ' num2str(n_2) ' [1]'])
disp(['T_im = ' num2str(T_hat_im) ' [K]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation intake manifold          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim('MeasModelCompare_intake.slx')            % Run the simulink model

if doPlot
    figure(7)
    %plot(AirStep.t,AirStep.alpha_ref,'b')
    hold on;
    plot(AirStep.t,AirStep.p_im,'r')
    plot(IntakeManifoldValidation_p_im.time,IntakeManifoldValidation_p_im.signals.values,'k')
    %axis([3 4 19 31])
    xlabel('Time [s]')
    ylabel('Intake manifold pressure [kPa]')
    legend('Measured','Model')    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computations for the injection time %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m_fi=m_dot_at*n_r./(N*n_cyl);
m_fi=m_dot_at*n_r./(N*n_cyl.*lambda_bc_cont*AFs);

A=[t_inj -ones(numel(t_inj),1)];

X  = A\m_fi;
c_fi=X(1);
t_0 = X(2)/c_fi;


m_hat_fi=c_fi*(t_inj-t_0);

rel_error=100*abs(m_fi-m_hat_fi)./m_fi;
index = [1:numel(rel_error)];

if doPlot %Here doPlot is used, avoids the plot if it is set to 0
    figure(8); clf; hold on
    plot(t_inj,m_fi,'ro')
    plot(t_inj,m_hat_fi,'k-')
    legend('m_fi','m_hat_fi','location','NorthWest')
    grid on
    xlabel('Fuel injection time [s]')
    ylabel('Mass of fuel injected [kg]')
    legend('Measured','Reference')
    
    
    figure(9); clf; hold on
    plot(t_inj, rel_error, 'r*')
    title('Fuel injection time: relative error')
    xlabel('Fuel injection time [s]')
    ylabel('Relative error [%]')
end
disp(' ')
disp('Fuel injector:')
disp(['c_fi = ' num2str(c_fi) ' [kg/s]'])
disp(['t_0 = ' num2str(t_0) ' [s]'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computations for the cylinder %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_ip = V_D*(p_em-p_im);
m_dot_fc = m_dot_at./(AFs*lambda_bc_cont);

a1 = m_dot_fc*n_r*q_LHV./N*(1-1/r_c^(gamma_air-1)).*min(1,lambda_bc_cont);
a2 = V_D*ones(length,1);
 
A = [a1 -a2];
b = n_r*2*pi*M_b+W_ip;

x = A\b;
 
n_ig_ch = x(1);
c_fr_0 = x(2);

M_hat_b = (A*x - W_ip)./(n_r*2*pi);

if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    figure(10); clf; hold on
    plot(M_b,M_hat_b,'ro')
    plot(M_b,M_b,'k-')
    grid on
    xlabel('Torque [Nm]')
    ylabel('Torque [Nm]')
    legend('Measured vs model','Measured vs measured')
    
    rel_error=100*abs(M_b-M_hat_b)./M_b;
    %index = [1:numel(rel_error)];
    
    figure(11); clf; hold on
    plot(M_b, rel_error, 'r*')
    title('Torque, low load: relative error')
    xlabel('Torque [Nm]')
    ylabel('Relative error [%]')
    axis([5 40 0 305])
    %axis([40 inf 0 25])
    
    figure(12); clf; hold on
    plot(M_b, rel_error, 'r*')
    title('Torque, high load: relative error')
    xlabel('Torque [Nm]')
    ylabel('Relative error [%]')
    %axis([5 40 0 305])
    axis([40 inf 0 25])
end
disp(' ')
disp('Cylinder:')
disp(['n_ig_ch = ' num2str(n_ig_ch) ' [1]']) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['c_fr_0 = ' num2str(c_fr_0) ' [N/m^2]'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Computations of the exhaust system    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Light-off time               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doPlot
    load ('lightOff.mat');
    figure(13)
    plot(lightOff.time,lightOff.lambda_bc_disc)
    hold on
    plot(lightOff.time,lightOff.lambda_ac_disc)
    xlabel('Time [s]')
    ylabel('Lambda')
    title('Light-off time')
end

%%   Computations of exhaust Temperature    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_dot_fc = m_dot_at./(AFs*lambda_bc_cont);
m_dot_exh = m_dot_at+m_dot_fc;

b = T_em;
A = [ones(length,1) sqrt(m_dot_exh)]; % Ändrat till roten ur

x = A\b;

T_0 = x(1);
k = x(2);

T_hat_em = A*x; 

if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    figure(14); clf; hold on
    plot(m_dot_exh,T_em,'ro')
    [ds,do] = sort(m_dot_exh);
    plot(ds,T_hat_em(do),'k-')
    grid on
    xlabel('Exhaust mass flow [kg/s]')
    ylabel('Temperature in exhaust manifold [K]')
    legend('Measured','Reference')
    
    rel_error=100*abs(T_em-T_hat_em)./T_em;
    
    figure(15); clf; hold on
    plot(m_dot_exh, rel_error, 'r*')
    title('Exhaust mass flow: relative error')
    xlabel('Exhaust mass flow [kg/s]')
    ylabel('Relative error [%]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Computations of exhaust pressure      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_dot_fc = m_dot_at./(AFs*lambda_bc_cont);
m_dot_exh = m_dot_at+m_dot_fc;
m_dot_amb=m_dot_exh;

est_points = [];
for i=1:length
    if (p_es(i)-p_amb)>0
        est_points = [est_points; i];
    end
end

A=m_dot_amb(est_points).^2;
b=(p_es(est_points)-p_amb)./(R_exh*T_em(est_points)./p_es(est_points));


x=A\b;

C_2=x(1);

m_hat_dot_amb = sqrt(b/x); 

if doPlot  %Here doPlot is used, avoids the plot if it is set to 0
    figure(16); clf; hold on
    plot(m_dot_amb(est_points),m_hat_dot_amb,'ro')
    plot(m_dot_amb(est_points), m_dot_amb(est_points),'k-')
    grid on
    ylabel('Ambient mass flow [kg/s]')
    xlabel('Ambient mass flow [kg/s]')
    legend('Measured vs model','Measured vs measured')
    
    rel_error=100*abs(m_hat_dot_amb-m_dot_amb(est_points))./m_dot_amb(est_points);
    
    figure(17); clf; hold on
    plot(m_dot_amb(est_points), rel_error, 'r*')
    title('relative error')
    ylabel('Ambient mass flow: relative error [%]')
    xlabel('Ambient mass flow [kg/s]')
end
disp(' ')
disp('Exhaust system:')
disp(['C_2 = ' num2str(C_2) ' [1/m^4]'])
disp(['t_0 = ' num2str(t_0) ' [s]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Computations of lambda         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[G stepTime] = min(abs(FuelStep.t-3.93885));

fake_lambda_cyl = [ones(stepTime,1); 0.8*ones(numel(FuelStep.t)-stepTime,1)];

tau_lambda = 0.21;
tau_d_bc = 0.05;
sim('MeasModelCompare_lambda.mdl')            % Run the simulink model

if doPlot
    figure(18)
    plot(FuelStep.t, fake_lambda_cyl,'b')
    hold on;
    plot(FuelStep.t, FuelStep.lambda,'r')
    plot(lambda_validation.time, lambda_validation.signals.values,'k')
    xlabel('Time [s]')
    ylabel('lambda')
    legend('Reference','Measured','Model')    
end
disp(' ')
disp('Lambda sensor:')
disp(['tau_d_bc = ' num2str(tau_d_bc) ' [s]'])
disp(['tau_lambda = ' num2str(tau_lambda) ' [s]'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% <<<< Until here is for project 1b / Next lines need to be
%%%%%%%%%%%%% <<<< uncommented for project 1c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize I/O abstraction layer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more information about the I / O block look under Project1c/ECU / "I / O abstraction layer"
% This can be used to do special tests with manually selected parameters:
% Engine speed, throttle reference, wastegate (project2) and pedal position.
% When all {property}_manual = 0, this block is disconnected by default.
N_e_manual = 0; N_e_step = 1; NINI = 2000; NEND = NINI;  NeST=30; NeSlope = 1; NeStartTime = 60; NeRampInit = 800;
alpha_REF_manual = 0; alphaINI = 0.0; alphaEND = alphaINI; alphaST=30;
wg_REF_manual = 0; wgINI = 100; wgEND = wgINI; wgST=30;
pedPos_manual = 0; pedINI = 0.2; pedEND = 1.0; pedST=30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate rolling and air resistance %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [0 410 1050]';                  % [N] Force vector
v = [0 80 160]';                    % [km/h] Speed vector
rho_air=p_amb/(R_air*T_amb);        % [kg/m^3] Air density
A=[VehicleMass*ones(size(v)) VehicleMass*v/3.6 0.5*VehicleArea*rho_air*(v/3.6).^2];
x=A\F;
% Rolling resistance coeficcients
c_r1=x(1);
c_r2=x(2);
%Air drag coeficcient
C_d=x(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load cycle Data     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
init_drivecycle;

% Loaded data for the European driving cycle:
% Driving Scenario: Matrix from which the following information is extracted
% S: Reference Speed [km/h]
% U: Clutch position [-]
% G: Selected gear [0-5]
% T: Time Vector [s]

% Forming vectors for the driver model
Clutch = [T U];
Gear   = [T max(G,1)];
Speed  = [T S];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Parameters for regulation           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_p =0.015;
K_I = 0.1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---------------------------------------------- %
% % Aftertreatment when the simulation is complete %
% % ---------------------------------------------- %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%     Light-Off computation      %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% lightOff = 37;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%      Compute emissions        %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inputs to calcEmissions
% % tout: Time Vector from Simulink
% % lambda: Continuous lambda
% % Distance: Distance traveled in meters
% % dmacAct: Mass air flow to the cylinder in kg / s
% % dmfAct: Fuel flow to the cylinder in kg / s
% % lightoff: Time in seconds until the light-Off
% 
% calcEmissions(tout, lambda, Distance, dmacAct, dmfcAct, lightOff);
% 
% Calculate fuel consumption
% fuelCons = 
%   
% disp(sprintf('Fuel Consumption: %1.2f [l/(10 mil)]',fuelCons))
