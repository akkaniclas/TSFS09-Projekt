%Prodject 2b model validtion tests

sim('Compressor.slx')

m_dot_Ccorr_val_mod_M = reshape([m_dot_cCorr_val_mod(1,:)' ; NaN],7,5);

if doPlot
    close all
    h = figure
    plot(m_dot_cCorr_M, PiC_M,'b-o',m_dot_Ccorr_mod_M, PiC_M, 'r--s',m_dot_Ccorr_val_mod_M, PiC_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('$\dot{m}_\textrm{c,corr}$ [kg/s]', 'interpreter', 'latex')
    ylabel('\Pi_c [-]')    
    saveas(h,'Figures\compressor_massflow_mod_val','png')
end

%% Turbine
sim('turbine.slx')

TFPmod_val_mod_M = reshape(TFPmod_val_mod(1,:)',6,5);

if doPlot
    close all
    h = figure
    plot(PiT_M, TFP_M,'b-o',PiT_M, TFPmod_M, 'r--s', PiT_M, TFPmod_val_mod_M, 'g*')
   % legend('Measured', 'Model', 'Location','northwest')
    xlabel('\Pi_t [-]')
    ylabel('TFP [kg/s K^{0.5}/kPa]')
    saveas(h,'Figures\turbine_mass_flow_validation','png')
end