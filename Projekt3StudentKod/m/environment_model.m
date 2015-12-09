% envirnoment_model.m
%
% Sets all mask parameters in the environment
% model to their default workspace variable names.
%


% Environment model
%%%%%%%%%%%%%%%%%%%

% positionvector in road slope map [m]
set_param([modelname],'pos','position');
% slopevector in road slope map [rad]
set_param([modelname],'slope','slope');
