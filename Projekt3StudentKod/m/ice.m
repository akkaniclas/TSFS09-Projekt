% ice.m
%
% Sets all mask parameters in the ice
% model to their default workspace variable names.
%

% engine mass [kg]
set_param([modelname '/ice'],'mass','0')
% engine inertia [kg m^2]
set_param([modelname '/ice'],'inertia','J_e')

% Engine parameters
set_param([modelname '/ice'],'tau_d','tau_d');
set_param([modelname '/ice'],'tau_e','tau_e');
set_param([modelname '/ice'],'init_T','init_T');
set_param([modelname '/local controller'],'T_max','T_max');
set_param([modelname '/local controller'],'N_max','N_max');
set_param([modelname '/local controller'],'N_idle','N_idle');
set_param([modelname '/local controller'],'P','IdleControllerP');
set_param([modelname '/local controller'],'I','IdleControllerI');
