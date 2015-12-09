% restoremask.m
%
% Sets all mask parameters to their default workspace variable names.
%

% The simulink model considered
% are the one in the following string.
name = 'car';

%
% Set vehicle model parameters
% ---------------------------

close_system(name);
open_system(name);

modelname = [name '/car_model/ice model'];
ice;
modelname = [name '/car_model/flywheel'];
flywheel;
modelname = [name '/car_model/clutch'];
clutch;
modelname = [name '/car_model/gearbox model'];
gearbox;
modelname = [name '/car_model/gear inertia'];
gearinertia;
modelname = [name '/car_model/propeller shaft'];
propellershaft;
modelname = [name '/car_model/final drive model'];
finaldrive;
modelname = [name '/car_model/final gear inertia'];
finalgearinertia;
modelname = [name '/car_model/drive shafts'];
driveshafts;
modelname = [name '/car_model/chassis model'];
chassis;

%
% Save and close system
% ---------------------
save_system(name);

% the wait is needed to allow save_system() to finish saving and
% removing backup files before close_system() is issued
pause(0.5);

close_system(name);
