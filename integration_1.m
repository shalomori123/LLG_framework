clear all; close all ; clc;
%our first try of integration between maxwell eq to LLG(runge kutte)

%% phisical constants;
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);


%%
% signal constents

lada0=800e-9;
W=2*pi*c/lada0;
S= -7;  % E ampllitod in dB
pulse_width_cycal= 5;
pulse_width= pulse_width_cycal*2*pi/W;
Tsteps=10*pulse_width/dt;



% grid and setup
% the matrix should be as (value X z_loction X time_step)

dt=(1/160)*(lada0/c); % the recomended sampeling freq is ten point per wave length is (1/40)*(lada0/c);
dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c

% the time spece defintion
L=17*lada0;
z=0:dz:L;
lenz=length(z);
H=zeros(3,lenz,Tsteps);
E=zeros(3,lenz,Tsteps);


% the matial defniton
z1=ceil(0.4*lenz); %location of material start
z2=ceil(0.9*lenz); %material end
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)];


e_r= [1;1]; %define eps_r
epsr= ones(2,lenz)+(e_r-1)*mat_vec;

sigma= [0 ; 0]; %define condatory 
sigma= sigma*mat_vec;
eaf= dt/2/eps*sigma./epsr;
conset_for_les_E= ((1-eaf)./(1+eaf));
conset_for_H= (0.5./epsr./(1+eaf));


ms= 3e5; % magnetization difniton
ms= [ms*ones(1,z2-z1,Tsteps); zeros(2,z2-z1,Tsteps)];
alphe= 80*0.00125;



t=0:dt:dt*(Tsteps);
% Pulse= 5e8*rectangularPulse(20,Tsteps,0:1:Tsteps);
Pulse =create_gassian_envelope(Tsteps,dt,pulse_width,S,0.3);
lim= max(Pulse);
% Pulse= [Pulse.' Pulse.';];

Pulse= amplified(W,dt,Pulse,pi/2, 0) ;
Pulse= [Pulse(:,1)';Pulse(:,2)'; zeros(1,Tsteps+1)];
plot3(t,Pulse(1,:), Pulse(2,:));
axis([ 0 t(end)   -lim lim -lim lim])
%% propagate

for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
    %n
    E(:,4,n)=E(:,4,n)+Pulse(:,n); %inject the input signal

    H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
    H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));

    ms(:,:,n+1) = runge_kutte_LLG_1(ms(:,:,n),H(:,z1+1:z2,n),dt,alphe); % solving LLG step
    H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms(:,:,n)- ms(:,:,n+1); % integration line

    E(2,1:lenz-1,n+1)=conset_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+conset_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
    E(1,1:lenz-1,n+1)=conset_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-conset_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
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
   H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms(:,:,n)- ms(:,:,n+1); % integration line

   E(2,1:lenz-1,n+1)=conset_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+conset_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
   E(1,1:lenz-1,n+1)=conset_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-conset_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
        
end

%%


H_p= squeeze((H(1,:,:)).^2 + (H(2,:,:)).^2 + (H(3,:,:)).^2).^0.5;
E_p= squeeze((E(1,:,:)).^2 + (E(2,:,:)).^2 + (E(3,:,:)).^2).^0.5;

% time_plot_2D_fields(E/max(max(max(E))),dz,z(z1),z(z2),2);
% time_plot_2_fields(squeeze(E(1,:,:))',squeeze(E(2,:,:))',dz,mat_vec*2.02-1.01,12,"E_x","E_y");
% time_plot_2_fields(H_p',E_p',dz,mat_vec*2.02-1.01,12,"|H|","|E|");

% plot_moments(ms,dt,5)

power_fator_to_measuer= 1;
Dx= 530;
direction=  "first" ; % 'last
Dt_void= measure_peak_moving_time(10,10+Dx,H,power_fator_to_measuer, direction);
Dt_moment= measure_peak_moving_time(z1+10,z1+10+Dx,H, power_fator_to_measuer, direction);
reflction_index= Dt_moment/Dt_void;





