% propellershaft.m
%
% Sets all mask parameters in the propeller shaft
% model to their default workspace variable names.
%

% mass [kg]
set_param([modelname],'mass','0');
% inertia [kg m^2]
set_param([modelname],'inertia','0');

% stiffness
set_param([modelname],'k_p','k_p');
% damping
set_param([modelname],'c_p','c_p');
