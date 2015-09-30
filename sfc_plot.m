function sfc_plot(torque,rps,sfc,lambda,fignmbr);
%
%
%          Musseldiagram och konturplot �ver sfc
%
%          funktion sfc_plot(torque,rpsec,sfc,lambda,fignmbr)
%
%               torque  = utmoment [Nm]
%               rps     = varv per sekund [rps]
%               sfc     = specifik br�nslef�rbrukning [kg/J]
%               lambda  = normerat luft-/br�nslef�rh�llande
%               fignmbr = figurnummer   
%
%
%          Funktionen skapar ett musseldiagram och en konturplot  
%          �ver specifik br�nslef�rbrukningen som funktion av 
%          moment och varvtal. Ol�mpliga m�tv�rden (de f�r vilka 
%	   lambda < 0.96 tas bort.
%
%          Figurerna ritas i figur 'fignmbr' och 'fignmbr'+1.
%
%          (Uppdaterad 2002-09-20)
%
%
%          

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Skriven av YN    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plockar bort element ur datam�ngden som �r op�litliga
% pga d�liga arbetsomr�den f�r sensorerna
len = length(torque);
index = setdiff(1:len,union(find(torque<15),find(lambda<0.95)));
M     = torque(index);
rpsec = rps(index);
sfc = sfc(index);


% Skapar ett rutn�t, och interpolerar sedan fram en 
% tredimensionell yta ovanp� detta
[X,Y] = meshgrid(min(M):(max(M)-min(M))/35:max(M),...
    min(rpsec):(max(rpsec)-min(rpsec))/35:max(rpsec));
Z = griddata(M,rpsec,sfc,X,Y);
Xmin = min(min(X));  Xmax = max(max(X));
Ymin = min(min(Y));  Ymax = max(max(Y));
Zmin = min(min(Z));  Zmax = max(max(Z));


% Plottar ett musseldiagram
if nargin<5
  fignmbr = 4;
end

figure(fignmbr); clf; 
mesh(X,Y*60,Z*1e6*3600);
view(42.5,25)
xlabel('Moment [Nm]');
ylabel('Varvtal [rpm]');
zlabel('sfc [g/kWh]');
axis([0 Xmax 0 Ymax*60 Zmin*3600e6 Zmax*3600e6]);
title('Specifik br�nslef�rbrukning');

% Ritar en figurplot
figure(fignmbr+1); clf;
sfcmin = min(sfc)*1e6*3600;
sfcmax = max(sfc)*1e6*3600;
v       = linspace(sfcmin,sfcmax,75);
vvalue  = v([1:15 20 30 50 75]);

%clabel(contour(X,Y*60,Z*1e6*3600,v),vvalue);
clabel(contour(Y*60,X,Z*1e6*3600,v),vvalue);
ylabel('Moment [Nm]');
xlabel('Varvtal [rpm]');
title('Specifik br�nslef�rbrukning [g/kWh]- niv�kurvor');



