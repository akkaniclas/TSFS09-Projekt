% flywheel.m
%
% Sets all mask parameters in the flywheel
% model to their default workspace variable names.
%

% mass [kg]
set_param([modelname],'mass','0')
% inertia [kg m^2]
set_param([modelname],'inertia','0')
% initial engine speed [rad/s]
set_param([modelname],'speedinit','init_e')
% initial angle [rad]
set_param([modelname],'angleinit','0')
