%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Student code for project2B   %%%
%%% TSFS09 - Fordonssystem       %%%
%%% Vaheed Nezhadali 2015-10-22  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all
clc
% slCharacterEncoding('windows-1252')
figpath='Figures/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run Init_Project1.m;                                             %run init_project1.m to load model parameters
run Project2A.m;                                            %Load turbo model parameters from the second project by running project2A.m
load('TqEvsNeMAP.mat')                                      % Load the torque v.s speed map
doPlot = 0;                                                 % should figures be plotted or not
sim_model_name = 'Project_template_2b';                        %specify name of the simulink model for project2B


%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Driver model parameters (should not be changed)
KpDriver        = 0.8;                         
KiDriver        = 0.05;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Boost_control = 1;          % Activating the boost controller block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting the controller parameters in the boostcontroller block
%%%%%%%%%%%%
KpThr = 1*10^-6;     
TiThr = 0.1;  
%%%%%%%%%%%%%
 
KpWg = 10e-6;
TiWg = 5;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Many parameters have same values as in project 1:
% V_em, V_im, V_es, PI_bl, r_c
% These values do not need to be redefined here and are called by running "Init_Project1.m"
% But since the displacement volume changes, some  parameters should be recalculated:
V_D_downsized = 1.2e-3;              % Wi will model a 1.2L engine...
n_cyl = 4;                           % ...with 4 cylinders.
V_d_downsised = V_D_downsized/n_cyl; % Cylinders volume
s   = (V_d_downsised*4/pi)^(1/3);    
B   = s;                             % assuming a "cubic" engine ("square bore")
l   = 3*s;
a   = s/2;
V_d = B^2/4*pi*s; 
V_D = n_cyl*V_d;    
n_r = 2;           

% En extra kontrollvolym används mellan kompressorn och trotteln:
V_ic      = 10e-3;        % Volym från intercooler, och rörsystem mellan kompressorn och insugsröret
T_ic      = T_im;         % Antag isoterm modell med T_ic=T_im
J_tc      = 1e-5;         % Turbocharger inertia. From Westin:2002 for Mitsubishi Heavy Industry TD04HL-15T [kg m^2]
c_tc_fric = 1e-6;         % Turbo shaft friction constant
Cd_wg     = 0.9;          % Assumed value for WG discharge constant
A_max_wg  = 0.02^2/4*pi;  % Measured approximation of maximum opening area of WG valve
dP_thrREF = 10e3;         % Default desired pressure loss over the throttle
tau_wg = 0.1;             % Wastegate actuator dynamics, estimated from measurement data

% over-writing the throttle model parameter value for lower idle w_ice
a_0= 0.6e-05;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit stupid parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_es=mean(T_es);
p_es=mean(p_es);
T_ic=mean(T_ic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 4 %% Throttle controller test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECU I/O
N_e_manual = 1; N_e_step = 1; NINI = 2400; NEND = NINI;              %constant engine speed 2400rpm
alpha_REF_manual = 0;    
wg_REF_manual = 0; 
pedPos_manual = 1; pedINI = 0.1; pedEND = 0.35; pedST=5;              % Step in pedal position
% Set the Driver Model outputs on manual(All swithes at "upper" position)
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','1')
set_param([sim_model_name,'/ECU/Boost Control/wg_feedback_switch'],'sw','0')


%%%%%%%%%%%%%%%%%%%%
%% Exercise 4b %%
%%%%%%%%%%%%%%%%%%%%

set_param([sim_model_name,'/ECU/Boost Control/thr_feedback_switch'],'sw','1') % Include throttle feedback in the simulations
set_param(sim_model_name,'StopTime','10')                                     % final simulation time
sim(sim_model_name) % Simulera modellen


if doPlot
    close all
    h=figure; clf;
    plot(t,p_im/1000)
    title('Throttle controller test');
    xlabel('Time [s]')
    ylabel('P_{im }[kPa]')
    hold on
    plot(p_im_ref.time,p_im_ref.signals.values/1000,'r')
    legend('P_{im}','P_{im,ref}','Location','northwest')
    grid on
    saveas(h,'Figures\4b','png')
end

%%%%%%%%%%%%%%%%%%%%
%% Exercise 4c %%
%%%%%%%%%%%%%%%%%%%%
set_param([sim_model_name,'/ECU/Boost Control/thr_feedback_switch'],'sw','0') % Exclude throttle feedback from the simulations
sim(sim_model_name) % Simulera modellen

%Skapa plottar.
if doPlot
    %close all
    h=figure; clf;
    plot(t,p_im/1000)
    title('Throttle controller test');
    xlabel('Time [s]')
    ylabel('Pim [KPa]')
    hold on
    plot(p_im_ref.time,p_im_ref.signals.values/1000,'r')
    legend('P_{im}','P_{im,ref}','Location','northwest')
    grid on
    saveas(h,'Figures\4c','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 5 %% Wastegate controller test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ECU I/O
N_e_manual = 1; N_e_step = 1; NINI = 2400; NEND = NINI;  % Konstnt varvtal 2400rpm
alpha_REF_manual = 0; 
wg_REF_manual = 0; 
pedPos_manual = 1; pedINI = 0.6; pedEND = 0.8; pedST=10;  % Steg i pedal position
% Driver Model settings
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','1')
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','1')
set_param([sim_model_name,'/ECU/Boost Control/thr_feedback_switch'],'sw','1') % Include throttle feedback in the simulations
set_param([sim_model_name,'/ECU/Boost Control/wg_feedback_switch'],'sw','1')  % Include throttle feedback in the simulations

% Simmulera och plotta
set_param(sim_model_name,'StopTime','15') % 10s simulering
sim(sim_model_name) % Simulera modellen

%% Plotta
if doPlot
    close all;
    % create a subplot(4,1,-) and plot pedal position, Pim_ref vs Pim, Pic_ref vs Pic, Engine torque 
    h=figure;
    title('Wastegate controller test');
    subplot(4,1,1)
    plot(t,Pedal_Pos)
    ylabel('Pos [-]')
    axis([5, 15,0.5,0.9 ])
    grid on
    subplot(4,1,2)
    plot(t,p_im/1000)
    ylabel('P_{im} [KPa]')
    axis([5, 15,130,200 ])
    hold on
    plot(p_im_ref.time,p_im_ref.signals.values/1000,'r')
    legend('P_{im}','P_{im,ref}', 'Location','northwest')
    grid on
    subplot(4,1,3)
    plot(t,p_ic/1000)
    axis([5, 15,130,200 ])

    ylabel('P_{ic} [KPa]')
    hold on
    plot(p_ic_ref.time,p_ic_ref.signals.values/1000,'r')
    legend('P_{ic}','P_{ic,ref}', 'Location','northwest')
    grid on
    
    subplot(4,1,4)
    plot(t,Tq_e)
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    grid on
    axis([5, 15,100,200 ])
    saveas(h,'Figures\5','png')
    
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Exercise 6 %% Maxmomentkurva + PiC(m_dot_cCorr)  %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_e_manual = 1; N_e_step = 0; NeSlope = 1; NeStartTime = 60; NeRampInit = 800;
alpha_REF_manual = 0; 
wg_REF_manual = 0;
pedPos_manual = 1; pedINI = 1; pedEND = 1;
% Driver Model - ställer om alla switchar till manuell.
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','1')    % Manuell växel
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','1')  % Manuel koppling
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','1') % Manuell gas
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','1')   % Manuell broms
% Simmulera och plotta
set_param(sim_model_name,'StopTime','5360') % The engine speed ramp is defined in line 141, calculate the time required for following the whole ramp.
sim(sim_model_name)

%%
load N_e_ramp_data
close all
% Plotta resultat
if doPlot
% plot the torque from engine vs torque from the interpolation curves on
% same figure for different engine speeds
    h = figure;
    plot(N_e,Tq_e)
    title('Torque when ramp in engine speed');
    xlabel('Engine speed [Rpm]')
    ylabel('Torque [Nm]')
    hold on
    plot(N_e,Tq_max_inter,'r')
    legend('Modelled engine torque','Interpolation table engine torque', 'Location','southwest')
    grid on
    saveas(h,'Figures\6a','png')

end


%Plotta kompressornstryckkvot(WcCorr)
if doPlot
    h=figure;clf;
    % plot the compressor mass flow as function of pressure ratio
    plot(m_dot_cCorr_sim.Data,PiC_sim,'b-.')
    hold on
    plot(m_dot_cCorr_M,PiC_M,'r')   % These come from turbomap and are calculated in project2A.m
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'Interpreter', 'latex');
    ylabel('\Pi_c [-]');
    legend('Modelled', 'Engine map', 'Location','southwest')
    grid on
    title('Compressor mass flow and pressure ratio');
    saveas(h,'Figures\6b','png')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 7  % Turbo vs NA engine performance
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% %% Exercise 7a %% Accelerationstest 60s %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% %ECU I/O
N_e_manual = 0;
alpha_REF_manual = 0; 
wg_REF_manual = 0;
pedPos_manual = 0;
% Driver Model - ställer om alla switchar till manuell.
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','1')    % Manuell växel
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','1')  % Manuel koppling
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','1') % Manuell gas
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','1')   % Manuell broms
% Simmulera och plotta
set_param(sim_model_name,'StopTime','60') % 60s simuleringa
sim(sim_model_name)                       % Simulerar modellen

%%
%plotta hastighetsprofil
if doPlot
% plot the vehicle speed and show how much time is required for 70-110 km/h
% acceleration
    h = figure;
    plot(t,VehicleSpeed)
    xlabel('time [s]', 'Interpreter', 'latex');
    ylabel('Vehicle Speed [km/h]');
    title('Acceleration test');
    saveas(h,'Figures\7a','png')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Exercise 7b %% Körcykel 600s %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ECU I/O
N_e_manual = 0;
alpha_REF_manual = 0;
wg_REF_manual = 0;
pedPos_manual = 0;
% Driver Model
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','0')
% Simmulera och plotta
set_param(sim_model_name,'StopTime','600') % 600s simulering
sim(sim_model_name)                        % Simulerar modellen

%% plotta resultatet.
if doPlot
    
    % after running the drive cycle, plot: engine torque, engine speed,
    % intake manifold pressure, lambda, intake manifold pressure.
    h = figure;
    subplot(4,1,1)
    plot(t,Tq_e)    
    grid on
    ylabel('Engine Torque [Nm]')
    
    subplot(4,1,2)
    plot(t,N_e)    
    grid on
    ylabel('Engine speed [RPM]')
    
    subplot(4,1,3)
    plot(t,p_im)    
    grid on
    ylabel('Intake manifold pressure [kPa]')
    subplot(4,1,4)
    plot(t,lambda_cyl)    
    grid on
    ylabel('Lambda, cylinder [-]')
    xlabel('Time [s]')
    
     title('Drive cycle test');
    saveas(h,'Figures\7b','png')
       
end

%% Beräkan emission och bränsleförbrukning
lightOff = 38;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Compute emissions        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs to calcEmissions
% tout: Time Vector from Simulink
% lambda: Continuous lambda
% Distance: Distance traveled in meters
% dmacAct: Mass air flow to the cylinder in kg / s
% dmfAct: Fuel flow to the cylinder in kg / s
% lightoff: Time in seconds until the light-Off

calcEmissions(t, lambda_cyl, Distance, dmacAct, dmfcAct, lightOff);

%Calculate fuel consumption
fuelCons = fuel(end)/0.75/Distance(end)*100000;
  
disp(sprintf('Fuel Consumption: %1.2f [l/(10 mil)]',fuelCons))


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Exercise 8%% Maximalt utmoment %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECU I/O
N_e_manual = 1; N_e_step = 1; NINI = 2500; NEND = NINI;  % Konstnt varvtal 2500rpm
alpha_REF_manual = 0; 
wg_REF_manual = 0; 
pedPos_manual = 1; pedINI = 0.0; pedEND = 1; pedST=30;  % Steg i pedal position
% Driver Model
set_param([sim_model_name,'/Driver Model/Manual Gear Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual Clutch Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual acc ped Switch'],'sw','0')
set_param([sim_model_name,'/Driver Model/Manual Brake Switch'],'sw','0')
% Simmulera och plotta
set_param(sim_model_name,'StopTime','70')  
sim(sim_model_name) % Simulera modellen

%% Plotta resultat
if doPlot
% plot the engine torque, alpha-thr and turbo speed in a subplot
    h = figure;
    subplot(3,1,1)
    plot(t,Tq_e)    
    grid on
    title('Maximal torque');
    axis([-inf,inf,-30,300])
    ylabel('Engine Torque [Nm]')
    
    subplot(3,1,2)
    plot(t,alpha_thr)    
    grid on
    ylabel('Throttle agle [-]')
    
    subplot(3,1,3)
    plot(t,w_tc)    
    grid on
    ylabel('Turbo speed [rad/s]')
    xlabel('Time [s]')
    
    saveas(h,'Figures\8','png')
end




