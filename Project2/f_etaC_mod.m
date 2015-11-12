function [etaC_mod] = f_etaC_mod(a,x)
% PiC_at_etaCmax = par(1);
% WcCorr_at_etaCmax = par(2);
% etaCmax = par(3);
% Q11 = par(4);
% Q22 = par(5);
% Q12 = par(6);

%x=[PiC m_dot_Ccorr];

etaC_mod = [];
for i=1:length(x)
    
Q = [a(4) a(6); a(6) a(5)];
X = [x(i,2)-a(2);sqrt(x(i,1)-1)-(a(1)-1)];   

tmp=a(3) - X'*Q*X;
%if tmp>0
    etaC_mod = [etaC_mod; tmp];
%else
 %   etaC_mod = [etaC_mod; 0.3];
%end
end

