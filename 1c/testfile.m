k=200;

PI = p_im./p_bef_thr;

est_points = [];
for i=1:length
    if PI(i)<0.73
        est_points = [est_points; i];
    end
end


p_im_test = p_im(est_points);
p_bef_test = p_bef_thr(est_points);
alpha_test = alpha(est_points);

m_dot_at_test = m_dot_at(est_points);

%sim('throttle.slx')


if  0 %Here doPlot is used, avoids the plot if it is set to 0
   
    rel_error=100*abs(m_dot_at(est_points)-simout.signals.values(1,:)')./m_dot_at(est_points);
    index = [1:numel(rel_error)];
    
    figure(4); clf; hold on
    %plot(index, rel_error, 'r*')
    plot(alpha(est_points), rel_error, 'r*')
    title('Mass air flow: relative error')
    xlabel('Throttle angle [%]')
    ylabel('Relative error [%]')
    
end
%%
k=20;

N_test = N(k);
t_inj_test = t_inj(k);
m_dot_fc_test = N(k)*n_cyl*m_fi(k)/n_r;

sim('fuel_injector.slx')
