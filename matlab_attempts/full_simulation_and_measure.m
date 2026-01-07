
eps=8.854e-12;
miu=4*pi*1e-7;
c=1/sqrt(eps*miu);

%%
e_r= [5;5]; %define eps_r
max_field= 5e9;  % [A/m]
pulse_width= 7;
sigma= [ 5.8e7 ;  5.8e7]*1e-3; %define condatory  [S/m]
ms0= [ 0 3e5 0 ];
alpha= 0.00125;
fea = pi/2;
close_system = false;

path= 'C:\maxwell-LLG\ws\simulation_result\';


[E,H,ms, z1, z2, z_len,dz,W] =gassian_simulation(e_r,sigma,ms0,alpha, ...
    max_field, pulse_width, fea, close_system,"time_cut",800,"len_matiral",400e-9,"len_in_lemda",1.5);

%
S= cross(E,H);
%%
% measere the optical linar ksi
ksi= zeros(1,z2- z1);

for k= 1:z2-z1
    ksi(k)= measure_optical_ksi(H,ms,W,dz/c/2,pulse_width,k,z1);
end

% measere the wavelen and velocity 
t_in= find(S(3,z2-z1,:)> 0.3 * max(max(max(S))),1 ,'first'); %entring to matiral time index
t_out= find(S(3,z2,:)> 0.3 * max(max(max(S))),1 ,'last'); %exit matiral time index
lamda_void= zeros(1,t_out-t_in);
lamda_in= zeros(1,t_out-t_in);
for t = 1: t_out-t_in 
    lamda_void(t)= measure_wave_length(H(:,:,t),[z2-z1-1 (z1-1)])*dz;
    lamda_in(t)= measure_wave_length(H(:,:,t),[z1 z2])*dz;
end
%%
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,z_len-z2)];
time_plot_2_fields(squeeze(E(1,:,:))',squeeze(E(2,:,:))',dz,mat_vec*2.02-1.01,20,"E_x","E_z");

%%

file_name= [path 'gasuian_mf_' num2str(max_field) '_w_' num2str(pulse_width) '_a_' num2str(alpha) '_ms0_' num2str(max(ms0)) '_close_' char(close_system+'0') '.mat'];
save(file_name, 'E' ,'H', 'ms', 'ksi' ,'lamda_in',"lamda_void","S")


