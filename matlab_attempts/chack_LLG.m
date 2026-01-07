clc
close all
% clear all
%%

eps=8.854e-12; %[F/m]
miu=4*pi*1e-7; %[H/m]
c=1/sqrt(eps*miu);


% void porprty
er= [1, 1];
sigma = [0 , 0];


ms = 3e5*[1, 0 , 0];
alpha = 0.7;
H_ext = [0, 0, 0];


lamda0 = 8e-7;
time_sycel = lamda0/c;
time_cut = 100;


max_feald = 9e9;
bwhp= 675e-12;
bwhp_to_us= 20*time_sycel;

%%
tic
    [puls ,dt, lamda_0] = gassian_gen_new(max_feald, bwhp_to_us, pi/2, ...
        time_cut= time_cut, ansition_factor= 5, relation_factor= 7, teta= pi/4);
    [E ,H ,ms_i ,z1 ,z2 ,lenz ,dz]  = simulation(er,sigma,ms, alpha, H_ext, ...
        puls, dt, false, false, ...
        len_matiral= 2*c*dt, len_simulation= lamda_0, start_matiral=0.6*lamda_0);

timeElapsed = toc
%%

t = dt*(0:(length(puls)-1));
hold on
plot(t, squeeze(real(ms_i(1,1,:))), DisplayName="x")
plot(t, squeeze(real(1.001*ms_i(2,1,:))), DisplayName="y")
plot(t, squeeze(real(ms_i(3,1,:))), DisplayName="z")
legend()


