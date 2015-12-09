function [cons, vol]= consumption(time, velocity, flow, dispflag)

%CONSUMPTION Calculates fuel consumption
% 
% [cons,vol] = consumption(time,velocity,flow,dispflag)
%
%  time       - time       [s]
%  velocity   - velocity [m/s]
%  flow       - fuel flow [g/s]
%  dispflag   - do (1) or do not (0) display result
%
%  cons - fuel consumption per 100 km [l/100km]
%  vol -  absolute fuel consumption [l]
%
%  The fuel density is set to 820 g/l

% Density of diesel [g/l]
rho_diesel = 820;

dist = trapz(time, velocity)*1e-4;
vol  = trapz(time, flow)/rho_diesel;
cons = vol/(dist*1e-1);

if(dispflag == 1)
    disp(sprintf('\nDistance        : %1.2f km\nFuel consumption: %1.2f l  or  %1.2f l/100km\n',10*dist,vol,cons));
end