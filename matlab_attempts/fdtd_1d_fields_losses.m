clear all; close all ; clc;

%% phisical constants;
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);

%% signal constents

lada0=800e-9;
W=2*pi*c/lada0;
S= -10;  % E ampllitod in dB
pulse_width=0.1e-12;




%% grid and setup

time_cut=1500;
Tsteps=60*time_cut;
dt=(1/time_cut)*(lada0/c); % the recomended sampeling freq is ten point per wave length.
dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c

% the time spece defintion
L=1.5*lada0;
z=0:dz:L;
lenz=length(z);
Hy=zeros(Tsteps,lenz);
Ex=zeros(Tsteps,lenz);


% the matial defniton
z1=ceil(0.7*lenz); %location of material start
z2=ceil(z1+ (40e-9/dz)); %material end
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)];

e_r=6; %define eps_r
epsr= ones(1,lenz)+(e_r-1)*mat_vec;

sigma= 5.97e7; %define condatory [S/m] copper
sigma= sigma*mat_vec;
eaf= dt/2/eps*sigma./epsr;
conset_for_les_E= (1-eaf)./(1+eaf);
conset_for_H= 0.5./epsr./(1+eaf);



t=0:dt:dt*(Tsteps);
% Pulse= 10^(S/10)*rectangularPulse(20,420,0:1:Tsteps);
Pulse =create_gassian_envelope(Tsteps,dt,0.2e-13,S,0.25);
Pulse= Pulse.*sin(W*t);
plot(Pulse)

%% propagate

for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
    %n
    Ex(n,4)=Ex(n,4)+Pulse(n); %inject the input signal
    for i=1:1:lenz-1
        %i
        Hy(n+1,i+1)=Hy(n,i+1)-0.5*(Ex(n,i+1)-Ex(n,i));
         Ex(n+1,i)=conset_for_les_E(i)*Ex(n,i)-conset_for_H(i)*(Hy(n+1,i+1)-Hy(n+1,i));
    end
end

for n=3:1:Tsteps-1
    %n
    Ex(n,4)=Ex(n,4)+Pulse(n); %inject the input signal
    for i=1:1:lenz-1
        %i
        Ex(n,lenz)=Ex(n-2,lenz-1);%right boundery condition
        Ex(n,1)=Ex(n-2,2);%left boundery condition
        Hy(n+1,i+1)=Hy(n,i+1)-0.5*(Ex(n,i+1)-Ex(n,i));
        Ex(n+1,i)=conset_for_les_E(i)*Ex(n,i)-conset_for_H(i)*(Hy(n+1,i+1)-Hy(n+1,i));
        
    end
end

%%


time_plot_epsr(Ex/max(max(Ex)),dz,2*epsr-5,40);
%time_plot(H,dz,'testing123');
%%

E_in= Ex(:,z1:z2);
[pks1,locs1]= findpeaks(E_in(Tsteps/2,:));
plot(E_in(Tsteps/2,:))


