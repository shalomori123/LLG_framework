clear all; close all ; clc;
%our first try of integration between maxwell eq to LLG(runge kutte)

%% phisical constants;
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);


%%
% sistem constents 

lada0=800e-9;
dt=(1/40)*(lada0/c); % the recomended sampeling freq is ten point per wave length is (1/40)*(lada0/c);

dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c
W=2*pi*c/lada0;


max_field= 5e9;  % [A/m]
pulse_width_cycal= 10;
pulse_width= pulse_width_cycal*2*pi/W;
Tsteps= 10*ceil(pulse_width/dt);
close_system= 0;


% grid and setup
% the matrix should be as (value X z_loction X time_step)



% the time spece defintion
L=20*lada0;
z=0:dz:L;
lenz=length(z);
H=zeros(3,lenz,Tsteps);
E=zeros(3,lenz,Tsteps);


% the matial defniton
z1=ceil(0.8*lenz); %location of material start
z2=ceil(0.9*lenz); %material end
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)];


e_r= [1.2;1.2]; %define eps_r
epsr= ones(2,lenz)+(e_r-1)*mat_vec;

sigma= [0 ; 0]; %define condatory 
sigma= sigma*mat_vec;
eaf= dt/2/eps*sigma./epsr;
const_for_les_E= ((1-eaf)./(1+eaf));
const_for_H= (0.5./epsr./(1+eaf));


ms0= [ 0 3e5 0 ]; % magnetization difniton
ms= zeros(3, z2-z1,Tsteps);
ms(:,:,1)= (ones(z2-z1,1)* ms0)' ;
alphe= 20*0.00125;



t=0:dt:dt*(Tsteps);
% Pulse= max_field*rectangularPulse(20,Tsteps-40,0:1:Tsteps);
Pulse =create_gassian_envelope(Tsteps,dt,pulse_width,max_field,0.27);
Pulse= amplified(W,dt,Pulse,pi/2, 0) ;
Pulse= [Pulse(:,1)';Pulse(:,2)'; zeros(1,Tsteps+1)];

lim= max(max(Pulse));
plot3(t,Pulse(1,:), Pulse(2,:));
axis([ 0 t(end) -lim lim -lim lim])
%% propagate

for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
    %n
    E(:,4,n)=E(:,4,n)+Pulse(:,n); %inject the input signal

    H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
    H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));

    ms(:,:,n+1) = runge_kutte_LLG_1(ms(:,:,n),H(:,z1+1:z2,n),dt,alphe); % solving LLG step
    if close_system == 1
        H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms(:,:,n)- ms(:,:,n+1); % integration line
    end
    E(2,1:lenz-1,n+1)=const_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+const_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
    E(1,1:lenz-1,n+1)=const_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-const_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
end

for n=3:1:Tsteps-1
    %n
   E(:,4,n)=E(:,4,n)+Pulse(:,n); %inject the input signal
   
   E(:,lenz,n)=E(:,lenz-1,n-2);%right boundery condition
   E(:,1,n)=E(:,2,n-2);%left boundery condition

   % instnd i loope

   H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
   H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));

   ms(:,:,n+1) = runge_kutte_LLG_1(ms(:,:,n),H(:,z1+1:z2,n),dt,alphe); % solving LLG step
   if close_system == 1
       H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms(:,:,n)- ms(:,:,n+1); % integration line
   end
   E(2,1:lenz-1,n+1)=const_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+const_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
   E(1,1:lenz-1,n+1)=const_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-const_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
        
end
%% senty check on masser gma using larmor
q=1.60217646e-19;                        %[cb]
miu=4*pi*1e-7;                           %[H/m]
g=2;                                     %[Landau factor]
me=9.1093821545*1e-31;                   %[gram] electron mass
gma_factor=1;
gma = gma_factor*g*q/(2*me);             %[rad/sec*T]

plot_moments(ms/max(max(max(ms))),dt,2)
plot(t(1:end-1),squeeze( H(1,z1+1,:)))

[pks,locs] = findpeaks(squeeze(ms(1,2,:)));
T_larmor= (locs(2:end)-locs(1:end-1))*dt;
W_larmor= 2*pi./T_larmor;
gamma_f= W_larmor/max_field/miu;
[mean(gamma_f) gma_factor*g*q/(2*me)]
%%

B= H(:,z1+1:z2,:)+ms; 
S= cross(E,H);
S_pulse= cross(Pulse(:,1:end-1),[Pulse(2,2:end) ;Pulse(1,2:end); zeros(1,Tsteps)]);
P_pulse= abs(2*Pulse(1,:).*Pulse(2,:));
H_p= squeeze((H(1,:,:)).^2 + (H(2,:,:)).^2 + (H(3,:,:)).^2).^0.5;
E_p= squeeze((E(1,:,:)).^2 + (E(2,:,:)).^2 + (E(3,:,:)).^2).^0.5;
B_p= squeeze((B(1,:,:)).^2 + (B(2,:,:)).^2 + (B(3,:,:)).^2).^0.5;
hold on
plot(squeeze(S(3,5,1:end)))
plot(squeeze(S(3,380,1:end)))
% plot(squeeze(S_pulse(3,:)))
% plot(squeeze(S(3,159,2000:end)))

% time_plot_2D_fields(E/max(max(max(E))),dz,z(z1),z(z2),2);
% time_plot_2_fields(squeeze(S(3,:,:))',squeeze(H(2,:,:))',dz,mat_vec*2.02-1.01,10,"S_z","H_y");
% time_plot_2_fields(H_p',E_p',dz,mat_vec*2.02-1.01,12,"|H|","|E|");

%% maseer the valosty off weve

% i masser in the power (grop velosty)
power_fator_to_measuer= 1;
Dx= 10;
direction=  "first" ; % 'last
t_void= measure_peak_moving_time([10 10+Dx],H,power_fator_to_measuer, direction);
t_moment= measure_peak_moving_time([z1+10 z1+10+Dx],H, power_fator_to_measuer, direction);
replction_index= (t_moment(2)-t_moment(1))/(t_void(2)- t_void(1));


t_in= measure_peak_moving_time(z1+1: z2,H, 0.1, 'first');
t_out= measure_peak_moving_time(z1+1: z2,H, 0.1, 'last');




% masser in the weve lenght (phassa valosity)
lamda0 =800e-9;
W=2*pi*c/lamda0;
lamda_void= measure_wave_length(H(:,:,t_in(1)),[1 (z1-1)],1)*dz;
lamda_in_x= measure_wave_length(H(:,:,t_out(1)),[z1 z2],1)*dz;
lamda_in_y= measure_wave_length(H(:,:,t_out(1)),[z1 z2],2)*dz;
C_measure= W*lamda_void/2/pi;
V_matral_x= W*lamda_in_x/2/pi;
V_matral_y= W*lamda_in_y/2/pi;



%% try to masser in the miu_r

miu_r_size= H_p(z1+1:z2,:)./B_p; % using power of fildes
miu_r_vec= H(:,z1+1:z2,:)./B; % as dignoz opertion
% miu_r_p= squeeze((miu_r_vec(1,:,:)).^2 + (miu_r_vec(2,:,:)).^2 + (miu_r_vec(3,:,:)).^2).^0.5;


% plot(t(1:12000),mio_r(10,:))

hold on
% plot(miu_r_size(10,t_in(10):t_out(10)))
% plot(miu_r_p(10,t_in(10):t_out(10)))
% plot(squeeze(miu_r_vec(1,10,t_in(10):t_out(10))))
% plot(squeeze(miu_r_vec(2,10,t_in(10):t_out(10))))
% plot(squeeze(miu_r_vec(3,10,t_in(10):t_out(10))))

hold off
miu_r_mean= mean(miu_r_size(10,t_in(10):t_out(10)));
mir_r_mean_vec= mean(miu_r_vec(:,:,t_in(10):t_out(end)),3);
mm= mean(mir_r_mean_vec,2);
% plot_moments(ms,dt,5)
% plot_moments(miu_r_vec(:,:,t_in(1):t_out(end)),dt,5)
%% trying to measure something close to miu ( miu = 1+ksi) using fft 
ksi_2= measure_optical_ksi(H,ms,W,dt,pulse_width_cycal,1,z1);

z_m= z1+1;
t_pik= measure_peak_moving_time(z_m,H,1, 'first');
t_pess= ceil(t_pik + 2.35*pulse_width/dt);

ms_fft= fft(ms(1,1,t_pik:t_pess));
ms_fft=ms_fft/length(ms_fft);
fs= 1/dt;
f = (0:length(ms_fft)-1)*fs/length(ms_fft);
f_optiy= W/2/pi;
k_op= find(f < f_optiy,1,"last");
ms_op= abs(ms_fft(k_op))+abs(ms_fft(k_op+1));
ksi= ms_op/max(max(max(H)));
% 
% plot(f,abs(squeeze(ms_fft)))
% 
% figure
% plot(squeeze(ms(1,1,t_pik:t_pess)))
