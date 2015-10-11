%% Accelerator position step for torque%%
load acc_step.mat;

if doPlot
    figure(1)
    subplot(2,1,1)
    plot(t,Pedal_Pos)
    axis([29.9 30.1 -inf inf])
    title('Pedal position step')
    ylabel('Pedal position [-]')
    subplot(2,1,2)
    plot(t,Tq_e)
    axis([29.9 30.1 -inf inf])
    xlabel('Time [s]') 
    ylabel('Torque [Nm]')
end

%% 
load full_throttle.mat
if doPlot
    figure(2)
    plot(t,VehicleSpeed)
    axis([-inf inf 70 110])
    xlabel('Time [s]')
    ylabel('Speed [km/h]')
    title('Full throttle, 4th gear')
end
%%
load drive_cycle_30
%%
K_p=0.002;
K_I=0.03;
%%
if doPlot
    figure(3)
    subplot(4,1,1)
    plot(t,VehicleSpeed)
    
    grid on
    ylabel('Vehicle speed [km/h]')
    subplot(4,1,2)
    plot(t,p_im)    
    grid on
    ylabel('Pressure, intake manifold [kPa]')
    subplot(4,1,3)
    plot(t,Tq_e)    
    grid on
    ylabel('Torque [Nm]')
    subplot(4,1,4)
    plot(t,lambda_cyl)    
    grid on
    ylabel('Lambda, cylinder [-]')
    xlabel('Time [s]')
    
    
    %hold on
    %plot(brake_pedal.time,brake_pedal.data)
    %plot(t,EngineSpeed)
end

%%
load drive_cycle_600
%%
if doPlot
    figure(3)
    subplot(4,1,1)
    plot(t,VehicleSpeed)
    
    grid on
    ylabel('Vehicle speed [km/h]')
    subplot(4,1,2)
    plot(t,p_im)    
    grid on
    ylabel('Pressure, intake manifold [kPa]')
    subplot(4,1,3)
    plot(t,Tq_e)    
    grid on
    ylabel('Torque [Nm]')
    subplot(4,1,4)
    plot(t,lambda_cyl)    
    grid on
    ylabel('Lambda, cylinder [-]')
    xlabel('Time [s]')
    
    e=emissions(lambda_cyl);
    %hold on
    %plot(brake_pedal.time,brake_pedal.data)
    %plot(t,EngineSpeed)
end