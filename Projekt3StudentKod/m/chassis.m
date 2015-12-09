% chassis.m

%

% Sets all mask parameters in the chassis

% model to their default workspace variable names.

%



% gravity constant [m/s^2]

set_param([modelname '/chassis'],'g','g');

% mass [kg]

set_param([modelname '/chassis'],'mass','mass');

% inertia [kg m^2]

set_param([modelname '/chassis'],'J_w','J_w');

% wheel radius [m]

set_param([modelname '/chassis'],'r_w','r_w');

% inital drive shaft speed [rad/s]

set_param([modelname '/chassis'],'init_d','init_d');

% air drag coefficient [-]

set_param([modelname '/chassis'],'c_w','c_w');

% maximum vehicle cross section area [m^2]

set_param([modelname '/chassis'],'A_a','A_a');

% Air density [kg/m^3]

set_param([modelname '/chassis'],'Rho_a','Rho_a');

% rolling resistance coefficient [-]

set_param([modelname '/chassis'],'c_r','c_r');

% brake gain [Nm]

% brake torque T is calculated as T = - brake_gain*B

% where B is the control signal [0,1]

set_param([modelname '/local controller'],'brake_gain', ...
                        'brake_gain');

% magic Bx

set_param([modelname '/chassis'],'Bx','Bx');

% magic Cx

set_param([modelname '/chassis'],'Cx','Cx');

% magic Dx

set_param([modelname '/chassis'],'Dx','Dx');

% magic Ex

set_param([modelname '/chassis'],'Ex','Ex');

