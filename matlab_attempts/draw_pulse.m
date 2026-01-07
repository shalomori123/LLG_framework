clc
clear all

path= 'C:\maxwell-LLG\ws\simulation_result\posibilty_invr\';
conset_file= 'const_for_all_running_for_date_2024_1_11.mat';
conset= open([path conset_file]);


eps_r = 1*  [1, 1] ;
sigma =0* [1 1];


ms0 = [0 0 1];
alpha =0.0001;
H_ext = [0 0 0];
close_system = 0 ;
edge_reflection = 1;



[pulse ,dt, lamda0] = gassian_gen(10000, pi/2, "relation_time" , 3e-15,"pulse_width", 0.2* 2.6685e-14);
[E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);

t= 1:length(pulse);
%plot3(t,pulse(1,:),pulse(2,:),'linewidth',2);



%plot3(t,squeeze(real(E(1,10,:))),squeeze(real(E(2,10,:))),"b","LineWidth",2.1)
%hold on
%plot3(t,squeeze(real(H(1,10,:))),squeeze(real(H(2,10,:))),"r","LineWidth",2.1)
time_plot_2D_fields(real(E),dz,z1,z2,5);

%%
path= 'C:\maxwell-LLG\ws\simulation_result\posibilty_invr\';
conset_file= 'const_for_all_running_for_date_2024_1_11.mat';

conset= open([path conset_file]);

polarstion_name = ["right", "left" ,"linar"];



ind_ms= 4;
ind_alp= 4;
ind_amp= 1;

file_name_mf_w_a_open= [path 'gasuian_mf_' num2str(conset.max_fild_run(ind_amp)) '_w_' num2str(conset.pulse_width_run) '_a_' num2str(conset.alpha_run(ind_alp)) '_ms0_' num2str(conset.ms_run(ind_ms)) '_close_0.mat'];
file_name_mf_w_a_close= [path 'gasuian_mf_' num2str(conset.max_fild_run(ind_amp)) '_w_' num2str(conset.pulse_width_run) '_a_' num2str(conset.alpha_run(ind_alp)) '_ms0_' num2str(conset.ms_run(ind_ms)) '_close_1.mat'];


simlation_op= open(file_name_mf_w_a_open);
simlation_cl= open(file_name_mf_w_a_close);
E1 = squeeze(simlation_op.H(1,:,:))';
E2 = squeeze(simlation_cl.H(1,:,:))';
E1_name = "Neglected Magnetic Light reaction";
E2_name = "Considered Magnetic Light reaction";
len_z_use= size(simlation_cl.S);
z1= floor(len_z_use(1)*conset.z1/conset.z_len);
z2= ceil(len_z_use(1)*conset.z2/conset.z_len);
dz= conset.dz;
Y_axis1=max(max(E1))*1.01;
Y_axis2= max(max(E2))*1.01;
%Y_axis=5;


n = 3000;
[N ,I]=size(E1);
z=0:dz:dz*(I-1);
z_lin= [z(1:z1) z(z1:z2) z(z2:end)];
line_matrial= 2.1*[zeros(1,z1) ones(1,z2-z1+1) zeros(1,I-z2+1)]-1.05;

ax = tiledlayout(2,1);
frame1= nexttile;
frame2= nexttile;


plot(frame1,z,E1(n,:),z_lin,Y_axis1 * line_matrial,'linewidth',2);
axis(frame1, [ 0 z(end)   -Y_axis1 Y_axis1])
plot(frame2,z,E2(n,:),z_lin, Y_axis2 *line_matrial,'linewidth',2);
axis(frame2, [ 0 z(end)   -Y_axis2 Y_axis2])
title(frame1, E1_name, "FontSize",15)
title(frame2, E2_name, "FontSize",15)
xlabel(frame1,"Z [m]")
ylabel(frame1, "Magnetic Field [A/m]")
xlabel(frame2,"Z [m]")
ylabel(frame2, "Magnetic Field [A/m]")
%%
z_lin_2= [z(1:z1) z(z1:z2-1) z(z2:end)];

figure
plot(z_lin,line_matrial)

        