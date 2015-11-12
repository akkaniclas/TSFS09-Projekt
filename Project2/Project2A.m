%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Student code for project2A   %%%
%%% TSFS09 - Fordonssystem       %%%
%%% Vaheed Nezhadali 2015-10-22  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figNr=1;
figpath='Figures/';

doPlot=1;
doExpFig=0;
doSimulinkFig=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      S�tter upp turbomotordata med antaganden     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  L�s in modellparametrar fr�n Projekt1 f�r att f� basmotorn alla parametrar

%clear all; close all;
load turboMap
%load Enginemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ladda ny m�tdata och parametrar %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NcCorr = comp.NcCorr;
m_dot_Ccorr = comp.m_dot_cCorr;
p01    = comp.p01;
T01    = comp.T01;
pCref  = comp.pCref;
TCref  = comp.TCref;
PiC    = comp.PiC;
etaC   = comp.etaC;
etaCmin=min(etaC/2); %will be used in Simulink models

p04  = turb.p04;
T03  = turb.T03;
TFP  = turb.TFP;
PiT  = turb.PiT;
TSP  = turb.TSP;
etaT = turb.etaT;
p03  = p04*PiT;
etaTmin=min(turb.etaT)/2; %will be used in Simulink models
% Compressor och Turbindiameter
dComp = 56e-3;
rComp = dComp/2;
dTurb = 52.2e-3;
rTurb = dTurb/2;
cp_exh = 1.2133e+03;            % [J/(kg*K)]
gamma_exh = 1.3;                % [J/(kg*K)]

cp_air =  980.0000;
gamma_air = 1.4000;

%Ber�knade storheter
Nc  = NcCorr*sqrt(T01/TCref);
m_dot_c  = m_dot_Ccorr*(p01/pCref)/(sqrt(T01/TCref));
Uc2 = rComp*2*pi*Nc/60;

m_dot_t  = TFP.*p04*1e-3.*PiT./sqrt(T03); %or  TFP.*p03*1e-3./sqrt(T03);
Nt  = TSP.*sqrt(T03);
BSR = 2*pi*Nt*rTurb./sqrt(2*cp_exh*T03*(1-PiT.^(-(gamma_exh-1)/gamma_exh)))/60;

% create matrices where every column represents the data from same turbo
% speed
m_dot_cCorr_M = reshape([m_dot_Ccorr ; NaN],7,5);
PiC_M         = reshape([PiC ; NaN],7,5);
NcCorr_M      = reshape([NcCorr ; NaN],7,5);
etaC_M        = reshape([etaC ; NaN],7,5);
Nc_M          = reshape([Nc ; NaN],7,5);
m_dot_c_M     = reshape([m_dot_c ; NaN],7,5);
Uc2_M         = reshape([Uc2 ; NaN],7,5);
% create matrices where every column represents the data from same turbo
% speed
TSP_M  = reshape(TSP,6,5);
TFP_M  = reshape(TFP,6,5);
PiT_M  = reshape(PiT,6,5);
etaT_M = reshape(etaT,6,5);
BSR_M  = reshape(BSR,6,5);
m_dot_t_M   = reshape(m_dot_t,6,5);
Nt_M   = reshape(Nt,6,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compressor modell    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPRESSOR Mass flow model
x=[PiC Uc2];

y=m_dot_Ccorr;

% Define reasonable start values on the parameters
x0 = [1 1];

% Define the nonlinear function
f_m_dot_Ccorr_mod = @(a,x)(a(1).*sqrt(1-(x(:,1)./((x(:,2).^2.*a(2)/(2*cp_air*T01)+1).^(gamma_air/(gamma_air-1)))).^2));
func = f_m_dot_Ccorr_mod;

par = lsqcurvefit(func, x0, x, y);
WcCorrMax = par(1);
PsiMax = par(2);

m_dot_Ccorr_mod = func(par,x);

m_dot_Ccorr_mod_M = reshape([m_dot_Ccorr_mod ; NaN],7,5);



if doPlot
    close all
    h = figure
    plot(m_dot_cCorr_M, PiC_M,'b-o',m_dot_Ccorr_mod_M, PiC_M, 'r--s')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'interpreter', 'latex')
    ylabel('\Pi_c [-]')    
    saveas(h,'Figures\compressor_mass_flow','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPRESSOR efficiency model

%par=[PiC_at_etaCmax WcCorr_at_etaCmax etaCmax Q11 Q22 Q12];
x=[PiC m_dot_Ccorr];
y=etaC;

x0=[1.977,0.08,0.8,90,-1,-6];

func=@f_etaC_mod;
par = lsqcurvefit(func, x0, x, y);

PiC_at_etaCmax = par(1);
WcCorr_at_etaCmax = par(2);
etaCmax = par(3);
Q11 = par(4);
Q22 = par(5);
Q12 = par(6);

etaC_mod=func(par,x);

etaC_mod_M = reshape([etaC_mod ; NaN],7,5);

if doPlot
    close all
    h = figure;
    plot(m_dot_cCorr_M, etaC_M,'b-o',m_dot_cCorr_M, etaC_mod_M, 'r--s')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'interpreter', 'latex')
    ylabel('\Pi_c [-]')
    saveas(h,'Figures\compressor_efficiensy','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Turbin modell    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURBIN Mass flow model


x=PiT;
y=TFP;

% Define reasonable start values on the parameters
c0ini = 0.0052; % Max(TFP) found at Pi=inf
c1ini = 2; % From FS-book, p. 159
x0 = [c0ini c1ini];

% Define the nonlinear function
f_TFPmod = @(a,x)(a(1).*sqrt(1-1./x.^a(2)));
func = f_TFPmod;

par = lsqcurvefit(func, x0, x, y);
k0 = par(1);
k1 = par(2);

TFPmod = func(par, x);

TFPmod_M = reshape(TFPmod,6,5);

if doPlot
    close all
    h = figure
    plot(PiT_M, TFP_M,'b-o',PiT_M, TFPmod_M, 'r--s')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('\Pi_t [-]')
    ylabel('TFP [kg/s K^{0.5}/kPa]')
    saveas(h,'Figures\turbine_mass_flow','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURBIN efficiency model
x=BSR;
y=etaT;

% Define reasonable start values on the parameters
c0ini = 0.8; % 
c1ini = 50; % 
x0 = [c0ini c1ini];

% Define the nonlinear function
f_etaT_BSR = @(a,x)(a(1).*(1 - ((x-a(2))./a(2)).^2));
func = f_etaT_BSR;

par = lsqcurvefit(func, x0, x, y);
etaTmax = par(1);
BSRmax = par(2)

% par=[BSRmax etaTmax];
% x=[Nt PiT]

% ----- find unknown model parameters using nonlinear least squares method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BMEP model    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par=[Cp0 Cp1];
% x=[EngineMap.p_im EngineMap.M_e]

% ----- find unknown model parameters using least squares method



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the validation figures for all models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doPlot
    %% your validation figure must be same as Figure 5.2 in the compendium.
    %% A validation plot for BMEP model is also required
end




