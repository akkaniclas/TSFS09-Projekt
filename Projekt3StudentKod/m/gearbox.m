% gearbox.m
%
% Sets all mask parameters in the gearbox model
% to their default workspace variable names.
%

% mass [kg]
set_param([modelname '/gearbox'],'mass','0');
% inertia [kg m^2]
set_param([modelname '/gearbox'],'inertia','0');
% Damping [Nm s/rad]
set_param([modelname '/gearbox'],'b_t','b_t');

% gear number to ratio and inertia table
set_param([modelname '/local controller'],'gearNumber','gearNumber'); 
set_param([modelname '/local controller'],'gearRatio','gearRatio');
set_param([modelname '/local controller'],'gearInertia','gearInertia');
