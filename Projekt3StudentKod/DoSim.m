% DoSim.m

close all

% simulate the reference model

% 


%

% Setup

% -----

% add path to scripts

addpath m;

% set simulation paramters

simopts = simset('Solver','ode23tb', ...
                 'DstWorkSpace','current','SrcWorkSpace','current');

% time to simulate

tspan = 30;

% load parameters

modelparam;

% if desired, restore the mask parameters with the following command

% restoremask;

% -----



%

% Simulate

% --------

sim('car',tspan,simopts);

% --------



%

% Store result

% ------------

% Time [s]

t = tout(end);

% Wheel speed [rad/s]

th_w = chassisout(:,2);

% Velocity [m/s]

v = chassisout(:,3);

% Acceleration [m/s^2]

a = chassisout(:,4);

% Slip ratio, SAE def.

kappaSAE = chassisout(:,5);

% converted to TSFS05 def.

kappa = 1-1./(1+kappaSAE);

% Longitudinal force [N]

Flong = chassisout(:,6);

% Position [m]

p = cumtrapz(tout,v);

% Acceleration pedal level [0,1]

pedal = controllerout(:,1);

% Brake pedal level [0,1]

brake = controllerout(:,2);

% Torque demand driver

pedal_driver = controllerout(:,5);

% Gear number

gear = controllerout(:,3);

% Engine speed [rpm]

rpm = (30/pi)*iceout(:,2);

% Engine torque [Nm]

torque = iceout(:,1);

% Drive shaft torque [Nm]

T_d = driveshaftsout(:,1);

% Drive shaft torsion [rad]

thDiff_d = driveshaftsout(:,2);

% ------------



%

% Plot

% ----

DoPlot;

% ----



% Plot observer performance



figure(3)
subplot(221); grid on; hold on;
plot(tout,thDiff_d,'b--','LineWidth',2);
plot(tout,obsstates(:,1),'r-')
title('Torsion')
xlabel('Time [s]')
ylabel('Torsion [rad]')

subplot(222); grid on; hold on;
plot(tout,rpm,'b--','LineWidth',2)
plot(tout,60/2/pi*obsstates(:,2),'r-')
title('Engine speed')
xlabel('Time [s]')
ylabel('Speed [rpm]')

subplot(223); grid on; hold on;
plot(tout,th_w,'b--','LineWidth',2);
title('Wheel speed');
xlabel('Time [s]');
ylabel('Wheel speed [rad/s]');
plot(tout,obsstates(:,3),'r-')

subplot(224); grid on; hold on;
plot(tout,torque,'b--','LineWidth',2)
plot(tout,obsstates(:,4),'r-')
title('Torque')
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('Measured','Observed')

% Plot Torque demand from driver and controller
figure(4)
subplot(121); hold on; grid on;
plot(tout,pedal,'blue-','LineWidth',2);
plot(tout,pedal_driver,'red-','LineWidth',2);
title('Torque demand from controller and driver');
xlabel('Time [s]');
ylabel('Torque [Nm]');
legend('Controller','Driver')
axis([0 inf -50 150]);

% Acceleration
subplot(122); hold on; grid on;
plot(tout,a,'LineWidth',2);
title('Vehicle');
ylabel('Acceleration [m/s^2]');
%axis([StartTime inf -inf inf]);

