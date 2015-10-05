function [EmissionmassBeforeCat,EmissionmassAfterCat] = ...
    calcemissions(tp, lambda3, dist, mairc,mfc, lightOffTime)
% [EmissionmassBeforeCat,EmissionmassAfterCat] = 
% calcemissions(t, lambda, dist, mairc, mfc, lightOffTime)
%
% Funktionen calcEmissions anv�nder f�ljande  indata f�r att
% best�mma emissionsniv�n f�re och efter cat.
% Indata:
%	t       tillh�rande tidsvektor med tidpunkter f�r varje sampel [s]
%   lambda  en vektor med labmdav�rden f�r varje tidpunkt i en k�rcykeln [s]
%	dist	den k�rda distansen i m vid varje sampel.
%	mairc	vektor med luftfl�de in i cylinder f�r varje tidpunkt [kg/s]
%	mfc	    -"-        br�nslef�de   -"- [kg/s]
%   lightOffTime Tiden i sekunder tills katalysatorn startar [s]
% Utdata:
%	emissionmassBC	emissionsniv� f�r [CO NOx HC] i [g/km] f�re cat
%   emissionmassAC	emissionsniv� f�r [CO NOx HC] i [g/km] efter cat
%
%  Av : Ingemar Andersson
%       Per Andersson
%
%  $Date: 2004/07/06 09:37:50 $
%  $Revision: 1.1.1.1 $

% Konvertera fr�n SI-enheter
dist  = dist/1e3;
mairc = mairc * 1e3;
mfc   = mfc * 1e3;

tic;

%disp('Ber�knar emissioner')
result = emissions(lambda3);
%disp('Utv�rderar')

% Utv�rderar data

loff = tp > lightOffTime; % Ger alla tidpunkter st�rre �n light off

MM = [28 14+16 16]; % Molmassor. Kompenserar f�r att emissions
% returnerar tal i procent f�r CO samt promille f�r NOx och HC.
% Observera att O2 och H2 inte anv�nds

MCO  = 12+16;
MNOx = 14+16;
MHC  = 12+4;
MO2  = 32;

% Total massa emissioner
EMF =  mairc + mfc;

Mtot = 28.5; % Ur Heywood.

len = length(result.emCO);

resultCO  = EMF .* (result.emCO*MCO)/(Mtot)/100;
resultNOx = EMF .* (result.emNOx*MNOx)/(Mtot)/1000;
resultHC  = EMF .* (result.emHC*MHC)/(Mtot)/1000;
resultO2  = EMF .* (result.emO2*MO2/Mtot);

EmissionmassBeforeCat = [resultCO resultNOx resultHC]'*[diff(tp); 0]/dist(length(dist));

disp('Emissioner f�re katalysator:')
disp(sprintf('\t CO  : %1.2f [g/km]', EmissionmassBeforeCat(1)))
disp(sprintf('\t HC  : %1.2f [g/km]', EmissionmassBeforeCat(3)))
disp(sprintf('\t NOx : %1.2f [g/km]', EmissionmassBeforeCat(2)))


%result2(i,:) = EMF(i) * (XvolCat(i,:).*MM)/(Mtot);
resultCO  = EMF .* (result.emCO*MCO)/(Mtot)/100.*(1-result.catEffCO.*loff);
resultNOx = EMF .* (result.emNOx*MNOx)/(Mtot).*(1-result.catEffNOx.*loff)/1000;
resultHC  = EMF .* (result.emHC*MHC)/(Mtot).*(1-result.catEffHC.*loff)/1000;
resultO2  = EMF .* (result.emO2*MO2/Mtot);
 
EmissionmassAfterCat = [resultCO resultNOx resultHC]'*[diff(tp); ...
0]/dist(length(dist));
disp(' ')

disp('Emissioner efter katalysator:')
disp(sprintf('\t CO  : %1.2f [g/km]\tGr�nsv�rde EURO 3 : 2.3  [g/km]   EURO 4 : 1.0  [g/km]', EmissionmassAfterCat(1)))
disp(sprintf('\t HC  : %1.2f [g/km]\tGr�nsv�rde EURO 3 : 0.20 [g/km]   EURO 4 : 0.10 [g/km]', EmissionmassAfterCat(3)))
disp(sprintf('\t NOx : %1.2f [g/km]\tGr�nsv�rde EURO 3 : 0.15 [g/km]   EURO 4 : 0.08 [g/km]', EmissionmassAfterCat(2)))
%toc


