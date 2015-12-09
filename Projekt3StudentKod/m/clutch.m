% clutch.m
%
% Sets all mask parameters in the clutch model
% to their default workspace variable names.
%

% mass [kg]
set_param([modelname],'mass','0');
% inertia [kg m^2]
set_param([modelname],'inertia','0');

% stiffness
set_param([modelname],'k_c','k_c');
% damping
set_param([modelname],'c_c','c_c');
