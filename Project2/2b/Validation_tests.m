%Prodject 2b model validtion tests
doPlot=1;

sim('Compressor.slx')
    close all

m_dot_Ccorr_val_mod_M = reshape([m_dot_cCorr_val_mod(1,:)' ; NaN],7,5);

if doPlot
    h = figure
    plot(m_dot_cCorr_M, PiC_M,'b-o',m_dot_Ccorr_mod_M, PiC_M, 'r--s',m_dot_Ccorr_val_mod_M, PiC_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'interpreter', 'latex')
    ylabel('\Pi_c [-]')    
    saveas(h,'Figures\compressor_massflow_mod_val','png')
end

etaC_val_mod_M = reshape([etaC_val_mod(1,:)' ; NaN],7,5);

if doPlot
    h = figure;
    plot(m_dot_cCorr_M, etaC_M,'b-o',m_dot_cCorr_M, etaC_mod_M, 'r--s', m_dot_cCorr_M, etaC_val_mod_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'interpreter', 'latex')
    ylabel('\eta_c [-]')
    saveas(h,'Figures\compressor_efficiensy','png')
end

if doPlot
    h = figure;
    plot(Tq_c_val_mod(1,:),'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('-')
    ylabel('Torque')
    saveas(h,'Figures\compressor_torque','png')
end


%% Turbine
sim('turbine.slx')
    close all

TFPmod_val_mod_M = reshape(TFPmod_val_mod(1,:)',6,5);

if doPlot
    h = figure
    plot(PiT_M, TFP_M,'b-o',PiT_M, TFPmod_M, 'r--s', PiT_M, TFPmod_val_mod_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('\Pi_t [-]')
    ylabel('TFP [kg/s K^{0.5}/kPa]')
    saveas(h,'Figures\turbine_mass_flow_validation','png')
end


etaT_val_mod_M = reshape(etaT_val_mod(1,:)',6,5);

if doPlot
    h = figure
    plot(BSR_M, etaT_M,'b-o',BSR_M, etaT_mod_M, 'r--s',BSR_M, etaT_val_mod_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('BSR [-]')
    ylabel('\eta_t [-]')
    saveas(h,'Figures\turbine_efficiency','png')
end

if doPlot
    h = figure;
    plot(Tq_t_val_mod(1,:),'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('-')
    ylabel('Torque')
    saveas(h,'Figures\turbine_torque','png')
end
