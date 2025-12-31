clc
clear all

%% maxwell ferfication - using faraday simaltion
%silcon prparty
e_r = 15.813 *[1 1];
sigma = 1.56e-3*[ 1 1];
ms0 = 1 ;
alpha = 0;

%lighe proprty
amplitude = 300;
run_fea= 0;
wave_lens = 1e-6*linspace(0.21, 0.82,10);

close_system = 0;
edge_reflection = 0;

%exepted result
complex_ind =3.9766 +0.030209i;

lamda_in = zeros(10,1);
lamda_out = zeros(10,1);

for ind = 1:10;
    [E,H,ms, z1, z2, z_len,dz,W, Tsteps] =faraday_simulation_complex(e_r,sigma ,ms0,alpha, amplitude, run_fea, 0 , close_system ...
                                ,"time_cut", 40 , "start_matiral",0.55,"len_matiral",3*wave_lens(ind), ...
                                "len_in_lemda",8 ,"simlation_time",1e-11 ,"lamd0", wave_lens(ind));
    
    
    lamda_in(ind)= measure_wave_length(real(E(:,:,end)),[z1 z2])*dz;
    lamda_out(ind)= measure_wave_length(real(E(:,:,end)),[5 z1-1])*dz;
end
%%
maesser_ind = wave_lens'./lamda_in;
plot(wave_lens, maesser_ind)

% verfication_result = abs( maesser_ind - abs(complex_ind))
%%
plot(squeeze(E))
%% glass

eps_r = 3.7*  [1, 1] ;
sigma =1e-8* [1 1];


ms0 = [0 0 1];
alpha =0.0001;
H_ext = [0 0 0];
close_system = 0 ;
edge_reflection = 0;
%%
[pulse ,dt, lamda0] = gassian_gen(10000, 0, "relation_time" , 3e-13,"pulse_width", 1* 2.6685e-14);
[E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection, ...
    "len_simulation", 50* lamda0, "start_matiral", 35.1* lamda0, "len_matiral", 10*lamda0);
%%
hold on
plot(squeeze(E(1,landz,:)),DisplayName="end")
plot(squeeze(E(1,landz-1,:)),DisplayName="end-1")
plot(squeeze(E(1,landz,:) - E(1,landz-1,:)),DisplayName="end - (end -1)")

%%
time_plot_2_fields(squeeze(E(1,:,:))',squeeze(E(2,:,:))',dz,z1, z2,9,"Ex","Ey");
%%
% phisical constants;
eps=8.854e-12; %[F/m]
miu=4*pi*1e-7; %[H/m]
c=1/sqrt(eps*miu);
heta = sqrt(miu/eps);

E_true = E/ heta; 
S= cross(E_true,H);
S_transrotm = squeeze(S(3,z1+1,:));
S_return = squeeze(S(3,5,:));

hold on
grid on
plot(S_return)
plot(S_transrotm)
S_in = sum(S_return(S_return > 0));
S_r = -sum(S_return(S_return < 0));
S_t = sum(S_transrotm);
% R = abs(min(S_return))/max(S_return)
% T = max(S_transrotm)/max(S_return)
% R = S_r /S_in
% T= S_t /S_in
% R+T
%%
S= cross(E,H);
P_t = sum(S_transrotm);
P_T_E = sum(sum(E(:, z1+1, :).^2));
P_in = sum(sum(pulse.^2));
% S_in = sum(sum(S(:,6,:)));
A= S_in/P_in
% T= P_t/P_in

n_entring = sqrt(4);
t_ = 2/(1+n_entring);
T_anlit= n_entring *(t_^2);
% T_n = 2*sqrt(3.7)/(1+sqrt(3.7))