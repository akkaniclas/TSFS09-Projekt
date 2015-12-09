% finaldrive.m
%
% Sets all mask parameters in the final drive
% model to their default workspace variable names.
%

% mass [kg]
set_param([modelname '/final_drive'],'mass','0');
% inertia [kg m^2]
set_param([modelname '/final_drive'],'inertia','J_f');
% conversion ratio[-]
set_param([modelname '/final_drive'],'i_f','i_f');
% damping [Nm s/rad]
set_param([modelname '/final_drive'],'b_f','b_f');
