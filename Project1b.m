%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Template file for Project 1b      %%%
%%% TSFS09 - Modeling and Control of  %%%
%%%          Engines and Drivelines   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This file is a template suitable for using in the first part of 
% Project 1b - the Analysis part. By writing commands heare instead
% of in the matlab prompt the calculations are easier to reproduce
% and change in case of error. To help you get started a few commands
% needed for Project 1 is already included.
%
%
% IMPORTANT: Make sure that you use SI-units for all calculations! 
% Other units can be better for plots to make them more easy to 
% interpret, for example its usally easier to have degrees in plots
% instead or radians.
%
%
% Don't forget to check that your results are reasonable!
% Something that is easy to control is for example the cylinder volume.
% 
%
%  Glöm inte att kontrollera att de resultat ni får är rimliga! Något som 
%  är ganska lätt att kontrollera är till exempel cylindervolymen.
%  Överensstämmer den slagvolym ni får med den slagvolym som anges i
%  laborationskompendiet och är kompressionsförhållandet rätt? 
%
%  (File uppdated 2015-08-28 by Andreas Thomasson)
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%          EXERCISE 1_high      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('pCylHigh.mat') % Load data file

pcyl   = pCylHigh.cyl_press;        % Cylinder pressure [Pa]   
M      = pCylHigh.M_b;        % Engine torque [Nm]
N      = pCylHigh.N;        % Engine speed [rpm]
theta  = pCylHigh.theta;        % Crank angle [rad]
r_c    = pCylHigh.r_c;        % Compression ratio
l      = 0.158;   % Vevstakens längd [m]
a      = 0.044;	  % Half the stroke length[m]
B      = 0.068;	  % Cylinder diameter [m]
ncyl   = 5;       % Number of cylinders
VD     = 1.6e-3;  % Total volume [m³] 
qLHV   = 44e6;    % Lower heating value of fuel 
R      = 280;     % Specific ga sconstant
lambda = 1;       % Normalized air/fuel ratio
gamma  = 1.3;     % Ratio of specific heats
AFs    = 15;      % Stochiometric air/fuel ratio
n_r    = 2;       % Strokes per cycle


%Calculate cylinder volume, V(theta)
V_c=VD/((r_c-1)*ncyl);
s_theta=a*cos(theta)+sqrt(l^2-a^2*sin(theta).^2);
V=V_c+pi*B^2/4*(l+a-s_theta);

%Calculate the avrage cycles of mesured cycles
pcyl_max=zeros(20:1);
for n=1:20
  pcyl_max(n)=max(pcyl(:,n));  
end
pcyl_max_mean=mean(pcyl_max);
[tmp, avrage_cycle]=min(abs(pcyl_max_mean-pcyl_max));

Wcalc = ncyl*trapz(V,pcyl(:,avrage_cycle)); % Calculate work for entire engine
Mcalc = Wcalc/(n_r*2*pi); % Calculate instantaneous torque for entire engine
Mmeas = mean(M(:,avrage_cycle)); % Calculate mean torque (measured) for entire engine
 
disp(['W - beräknat: ', num2str(Wcalc)]) 
disp(['M - beräknat: ', num2str(Mcalc)]) 
disp(['M - uppmätt:  ', num2str(Mmeas)]);

figure(1); clf; 
plot(V*1e3,pcyl*1e-5)
title('pV-diagram (all cycles, high load)');
xlabel('Volume [dm^3]');
ylabel('Pressure [bar]');


% To store the figure as an .eps you can use the command:
% print -f1 -depsc2 filnamn.eps
%    -f1    : Figure nummber 1
%    -depsc2: Format eps 2 with color (c)
% To store the figure as a jpeg you can use:
% print -djpeg filnamn.jpg
% Are you unsure or want more information, type "help print" 
% in the matlab  prompt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%          EXERCISE 1_low       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('pCylLow.mat') % Load data file

pcyl   = pCylLow.cyl_press;        % Cylinder pressure [Pa]   
M      = pCylLow.M_b;        % Engine torque [Nm]
N      = pCylLow.N;        % Engine speed [rpm]
theta  = pCylLow.theta;        % Crank angle [rad]
r_c    = pCylLow.r_c;        % Compression ratio
l      = 0.158;   % Vevstakens längd [m]
a      = 0.044;	  % Half the stroke length[m]
B      = 0.068;	  % Cylinder diameter [m]
ncyl   = 5;       % Number of cylinders
VD     = 1.6e-3;  % Total volume [m³] 
qLHV   = 44e6;    % Lower heating value of fuel 
R      = 280;     % Specific ga sconstant
lambda = 1;       % Normalized air/fuel ratio
gamma  = 1.3;     % Ratio of specific heats
AFs    = 15;      % Stochiometric air/fuel ratio
n_r    = 2;       % Strokes per cycle

%Calculate cylinder volume, V(theta)
V_c=VD/((r_c-1)*ncyl);
s_theta=a*cos(theta)+sqrt(l^2-a^2*sin(theta).^2);
V=V_c+pi*B^2/4*(l+a-s_theta);

%Calculate the avrage cycles of mesured cycles
pcyl_max=zeros(20:1);
for n=1:20
  pcyl_max(n)=max(pcyl(:,n));  
end
pcyl_max_mean=mean(pcyl_max);
[tmp, avrage_cycle]=min(abs(pcyl_max_mean-pcyl_max));

Wcalc = ncyl*trapz(V,pcyl(:,avrage_cycle)); % Calculate work for entire engine
Mcalc = Wcalc/(n_r*2*pi); % Calculate instantaneous torque for entire engine
Mmeas = mean(M(:,avrage_cycle)); % Calculate mean torque (measured) for entire engine

disp(['W - beräknat: ', num2str(Wcalc)]) 
disp(['M - beräknat: ', num2str(Mcalc)]) 
disp(['M - uppmätt:  ', num2str(Mmeas)]);

figure(2); clf; 
plot(V*1e3,pcyl*1e-5)
title('pV-diagram (all cycles, low load)');
xlabel('Vloume [dm^3]');
ylabel('Pressure [bar]');


% To store the figure as an .eps you can use the command:
% print -f1 -depsc2 filnamn.eps
%    -f1    : Figure nummber 1
%    -depsc2: Format eps 2 with color (c)
% To store the figure as a jpeg you can use:
% print -djpeg filnamn.jpg
% Are you unsure or want more information, type "help print" 
% in the matlab  prompt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          EXERCISE 2           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Calculations done on pCylHigh
%load('EngineMapTSFS09.mat')
load('pCylHigh.mat') % Load data file
pcyl   = pCylHigh.cyl_press;        % Cylinder pressure [Pa] 
r_c=pCylHigh.r_c;
theta  = pCylHigh.theta;        % Crank angle [rad]
AF  = 14;        % Crank angle [rad]
qLHV   = 44e6;    % Lower heating value of fuel 
R      = 280;     % Specific ga sconstant
gamma  = 1.3;     % Ratio of specific heats
VD     = 1.6e-3;  % Total volume [m³] 
ncyl   = 5;       % Number of cylinders
a      = 0.044;	  % Half the stroke length[m]
B      = 0.068;	  % Cylinder diameter [m]
l      = 0.158;   % Vevstakens längd [m]
m_dot_at=pCylHigh.m_dot_at;
avrage_cycle=14;


Ratio_m_f_m_t = 1/(1+AF);
c_v=R/(gamma-1);
%q_lhv*(mf/mt * cv)=E
E=qLHV*Ratio_m_f_m_t/c_v;

%qin=mf*q_lhv/mt, qin=cv(T3-T2)  T3-T2=q_lhv*mf/(mt * cv)
%ideal gas cv=R/(gamma-1)

%Calculate cylinder volume, V(theta)
V_c=VD/((r_c-1)*ncyl);
s_theta=a*cos(theta)+sqrt(l^2-a^2*sin(theta).^2);
V=V_c+pi*B^2/4*(l+a-s_theta);

V_theta_high_low=V(length(V)/4:3*length(V)/4+1);

T_1=pCylHigh.T_im;
p_1=pCylHigh.p_im(720/4,1)*0.8;

T_2 = r_c^(gamma-1)*T_1;
p_2=r_c^gamma*p_1;

T_3=E+T_2;
p_3=T_3*p_2/T_2;

T_4=T_3/r_c^(gamma-1);
p_4=p_3/r_c^gamma;


%Vector for pressure ideal otto
p_theta=zeros(1,ceil(length(V_theta_high_low)));
%Calc lower pressure curve
k_1=p_1*V_theta_high_low(1)^gamma;
p_theta(1:ceil(length(V_theta_high_low)/2))=k_1./(V_theta_high_low(1:ceil(length(V_theta_high_low)/2)).^gamma);

%Calc higher pressure curve
k_1=p_3*V_theta_high_low(ceil(length(V_theta_high_low)/2))^gamma;
p_theta(ceil(length(V_theta_high_low)/2):end-1)=k_1./(V_theta_high_low(ceil(length(V_theta_high_low)/2):end-1).^gamma);

p_theta(end)=p_theta(1);

figure; clf;
plot(V_theta_high_low*1e3,p_theta*1e-5)
hold on
plot(V*1e3,pcyl(:,avrage_cycle)*1e-5)
hold off
title('pV-diagram');
xlabel('Volume [dm^3]');
ylabel('Pressure [bar]');

%%
%c
qLHV_fit=qLHV*0.624;

E=qLHV_fit*Ratio_m_f_m_t/c_v;

T_1=pCylHigh.T_im;
p_1=pCylHigh.p_im(720/4,1)*0.8;

T_2 = r_c^(gamma-1)*T_1;
p_2=r_c^gamma*p_1;

T_3=E+T_2;
p_3=T_3*p_2/T_2;

T_4=T_3/r_c^(gamma-1);
p_4=p_3/r_c^gamma;

%Q_loss=(T_4-T_1)*m_tot*c_v;

%Vector for pressure ideal otto
p_theta=zeros(1,ceil(length(V_theta_high_low)));
%Calc lower pressure curve
k_1=p_1*V_theta_high_low(1)^gamma;
p_theta(1:ceil(length(V_theta_high_low)/2))=k_1./(V_theta_high_low(1:ceil(length(V_theta_high_low)/2)).^gamma);

%Calc higher pressure curve
k_1=p_3*V_theta_high_low(ceil(length(V_theta_high_low)/2))^gamma;
p_theta(ceil(length(V_theta_high_low)/2):end-1)=k_1./(V_theta_high_low(ceil(length(V_theta_high_low)/2):end-1).^gamma);

p_theta(end)=p_theta(1);

figure; clf;
plot(V_theta_high_low*1e3,p_theta*1e-5)
hold on
plot(V*1e3,pcyl(:,avrage_cycle)*1e-5)
hold off
title('pV-diagram');
xlabel('Volym [dm^3]');
ylabel('Tryck [bar]');

Ratio_ideal_real=trapz(V,pcyl(:,avrage_cycle))/trapz(V_theta_high_low,p_theta);
disp(['Real_Otto/Ideal_Otto: ', num2str(Ratio_ideal_real)]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          EXERCISE 3           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%% Load the engine map data
load ('EngineMapTSFS09.mat')
m_dot_at=EngineMap.m_dot_at;
lambda=EngineMap.lambda;
AFs=EngineMap.airfuel.AFs;
M_e=EngineMap.M_e;
N=EngineMap.N;
qLHV=EngineMap.airfuel.q_LHV;

m_dot_f=m_dot_at./(lambda*AFs);
SFC=m_dot_f./(M_e*2*pi.*N);
%FC=1/(eta_f*qLHV);

sfc_plot(M_e,N,SFC,lambda,1);

% Tips: Use the function sfc_plot


