% DoPlot.m
%
% Plots simulation results
%

% Show data from StartTime
StartTime = t1-2;

% create figure 1
h = figure(1); %clf;
%scrsz = get(0,'ScreenSize');
%ppos = [0.6345 0.6345 19.72 28.41];
% Position in [left bottom width height]
%pos = [0 scrsz(4)*0.2 scrsz(3)*0.5 scrsz(4)*0.8];
%set(h,'Color','white','Position',pos,'PaperUnits','centimeters','PaperType','a4','PaperOrientation','Portrait','PaperPosition',ppos)

% Accelerator and brake levels
subplot(321); hold on; grid on;
plot(tout,pedal,'blue-','LineWidth',2);
plot(tout,brake,'red-','LineWidth',1);
title('Controller');
xlabel('Time [s]');
ylabel('Level');
legend('acc','brake');
%axis([StartTime inf -inf inf]);

% Gear
subplot(322); hold on; grid on;
plot(tout,gear,'LineWidth',2);
title('Controller');
xlabel('Time [s]');
ylabel('Gear number');
%axis([StartTime inf min(gear)-1 max(gear)+1]);

% Engine torque
subplot(323); hold on; grid on;
plot(tout,torque,'LineWidth',2);
title('Engine');
ylabel('Torque [Nm]');
%axis([StartTime inf -inf inf]);

% Engine speed
subplot(324); hold on; grid on;
plot(tout,rpm,'LineWidth',2);
title('Engine');
ylabel('Speed [rpm]');
%axis([StartTime inf -inf inf]);

% Velocity
subplot(325); hold on;  grid on;
plot(tout,3.6*v,'black','LineWidth',2);
title('Vehicle');
ylabel('Velocity [km/h]');
%axis([StartTime inf -inf inf]);

% Acceleration
subplot(326); hold on; grid on;
plot(tout,a,'LineWidth',2);
title('Vehicle');
ylabel('Acceleration [m/s^2]');
%axis([StartTime inf -inf inf]);

% create figure 2
h = figure(2); % clf;
%scrsz = get(0,'ScreenSize');
%ppos = [0.6345 0.6345 19.72 28.41];
% Position in [left bottom width height]
%pos = [scrsz(3)*0.5 scrsz(4)*0.2 scrsz(3)*0.5 scrsz(4)*0.8];
%set(h,'Color','white','Position',pos,'PaperUnits','centimeters','PaperType','a4','PaperOrientation','Portrait','PaperPosition',ppos)

% Wheel speed
subplot(321); hold on; grid on;
plot(tout,th_w,'LineWidth',2);
title('Wheel speed');
xlabel('Time [s]');
ylabel('[rad/s]');
%axis([StartTime inf -inf inf]);

% Slip ratio
subplot(322); hold on; grid on;
plot(tout,kappa,'LineWidth',2);
title('Slip ratio');
xlabel('Time [s]');
ylabel('[-]');
%axis([StartTime inf -inf inf]);

% Long.force
subplot(323); hold on; grid on;
plot(tout,Flong,'LineWidth',2);
title('Longitudinal force');
xlabel('Time [s]');
ylabel('[N]');
%axis([StartTime inf -inf inf]);

% Longitudinal force characteristic
subplot(324); hold on; grid on;
uSAE = min(-0.5,min(kappaSAE)):0.05:max(10,max(kappaSAE));
Fx = (1/(muX*Fz))*Dx*sin(Cx*atan(Bx*uSAE-Ex.*(Bx*uSAE-atan(Bx*uSAE))));
u = 1-1./(1+uSAE);
plot(u,Fx,'LineWidth',1); plot(kappa(1),interp1(u,Fx,kappa(1)),'.');
h = legend('Characteristic','Operation points');
set(h,'Location','SouthEast','FontSize',6);
plot(kappa,interp1(u,Fx,kappa),'.');
title('Longitudinal force characteristic');
xlabel('\kappa');
ylabel('F_x/\mu F_z');

% Drive shaft torque
subplot(325); hold on; grid on;
plot(tout,T_d,'LineWidth',2);
title('Drive shafts');
xlabel('Time [s]');
ylabel('Torque [Nm]');
%axis([StartTime inf -inf inf]);

% Drive shaft torsion
subplot(326); hold on; grid on;
plot(tout,thDiff_d,'LineWidth',2);
title('Drive shafts');
xlabel('Time [s]');
ylabel('Torsion [rad]');
%axis([StartTime inf -inf inf]);
