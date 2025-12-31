clear all; close all ; clc;
% usinig function
%%
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);

%%

e_r =[1; 1];
sigma0 = [0;0];
ms0 = [3e9 0 0];
alphe= 20*0.00125;
max_feled = 3e10;
pulse_width = 10;
fea= pi/2;
%%
[E,H,ms,mat_vec, dz] = gassian_simulation(e_r,sigma0,ms0,alphe, max_feled ...
    , pulse_width, fea, ...
    "time_cut",120);
dt= dz/2/c;
z1= find(mat_vec==1,1,'first');
z2= find(mat_vec==1,1,'last');

B= H(:, z1:z2,:) +ms;
%%

H_p= squeeze((H(1,:,:)).^2 + (H(2,:,:)).^2 + (H(3,:,:)).^2).^0.5;
E_p= squeeze((E(1,:,:)).^2 + (E(2,:,:)).^2 + (E(3,:,:)).^2).^0.5;
B_p= squeeze((B(1,:,:)).^2 + (B(2,:,:)).^2 + (B(3,:,:)).^2).^0.5;
miu_r_size= H_p(z1:z2,:)./B_p;
miu_r_vec= H(:,z1:z2,:)./B;
miu_r_p= squeeze((miu_r_vec(1,:,:)).^2 + (miu_r_vec(2,:,:)).^2 + (miu_r_vec(3,:,:)).^2).^0.5;
% plot(t(1:12000),mio_r(10,:))

% time_plot_2D_fields(E/max(max(max(E))),dz,z(z1),z(z2),2);
% time_plot_2_fields(squeeze(H(1,:,:))',squeeze(H(2,:,:))',dz,mat_vec*2.02-1.01,12,"H_x","H_y");
% time_plot_2_fields(H_p',E_p',dz,mat_vec*2.02-1.01,12,"|H|","|E|");



power_fator_to_measuer= 1;
Dx= 10;
direction=  "first" ; % 'last
t_void= measure_peak_moving_time([10 10+Dx],H,power_fator_to_measuer, direction);
t_moment= measure_peak_moving_time([z1+10 z1+10+Dx],H, power_fator_to_measuer, direction);
t_in= measure_peak_moving_time(z1+1: z2,H, 0.1, 'first');
t_out= measure_peak_moving_time(z1+1: z2,H, 0.1, 'last');
replction_index= (t_moment(2)-t_moment(1))/(t_void(2)- t_void(1));


hold on
plot(miu_r_size(10,t_in(10):t_out(10)))
% plot(miu_r_p(10,t_in(10):t_out(10)))
plot(squeeze(miu_r_vec(1,10,t_in(10):t_out(10))))
plot(squeeze(miu_r_vec(2,10,t_in(10):t_out(10))))
% plot(squeeze(miu_r_vec(3,10,t_in(10):t_out(10))))

hold off
miu_r_mean= mean(miu_r_size(10,t_in(10):t_out(10)));
mir_r_mean_vec= mean(miu_r_vec(:,:,t_in(10):t_out(end)),3);
mm= mean(mir_r_mean_vec,2);
% plot_moments(ms,dt,5)
% plot_moments(miu_r_vec(:,:,t_in(1):t_out(end)),dt,5)

%%
lamda0 =800e-9;
W=2*pi*c/lamda0;
lamda_void= measure_wave_length(H(:,:,t_in(1)),[1 (z1-1)],1)*dz;
lamda_in_x= measure_wave_length(H(:,:,t_out(1)),[z1 z2],1)*dz;
lamda_in_y= measure_wave_length(H(:,:,t_out(1)),[z1 z2],2)*dz;
C_measure= W*lamda_void/2/pi;
V_matral_x= W*lamda_in_x/2/pi;
V_matral_y= W*lamda_in_y/2/pi;