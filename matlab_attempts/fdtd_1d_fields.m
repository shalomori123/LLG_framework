clear all; close all ; clc;

%% phisical constants;
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);

%% signal constents

lada0=1.5e-6;
W=2*pi*c/lada0;
S= -10;  % E ampllitod in dB
pulse_width=0.1e-12;
Tsteps=1800;



%% grid and setup


dt=(1/40)*(lada0/c); % the recomended sampeling freq is ten point per wave length.
dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c

% the time spece defintion
L=20*lada0;
z=0:dz:L;
lenz=length(z);
Hy=zeros(Tsteps,lenz);
Ex=zeros(Tsteps,lenz);


% the matial defniton
z1=ceil(0.4*lenz); %location of material start
z2=ceil(0.6*lenz); %material end
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)];

e_r=5; %define eps_r
epsr= ones(1,lenz)+(e_r-1)*mat_vec;

sigma= 0; %define condatory 
sigma= sigma*mat_vec;




% Pulse= 10^(S/10)*rectangularPulse(4,404,1:1:Tsteps);
Pulse =create_gassian_envelope(Tsteps,dt,0.2e-13,S,0.2);
plot(0:dt:dt*Tsteps,Pulse)
%% propagate

for t=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
    %n
    Ex(t,4)=Ex(t,4)+Pulse(t); %inject the input signal
    for i=1:1:lenz-1
        %i
        Hy(t+1,i+1)=Hy(t,i+1)-0.5*(Ex(t,i+1)-Ex(t,i));
        Ex(t+1,i)=Ex(t,i)-0.5/epsr(i)*(Hy(t+1,i+1)-Hy(t+1,i));
    end
end

for t=3:1:Tsteps-1
    %n
    Ex(t,4)=Ex(t,4)+Pulse(t); %inject the input signal
    for i=1:1:lenz-1
        %i
%         Ex(t,lenz)=Ex(t-2,lenz-1);%right boundery condition
%         Ex(t,1)=Ex(t-2,2);%left boundery condition
        Hy(t+1,i+1)=Hy(t,i+1)-0.5*(Ex(t,i+1)-Ex(t,i));
        Ex(t+1,i)=Ex(t,i)-0.5/epsr(i)*(Hy(t+1,i+1)-Hy(t+1,i));
        
    end
end

%%


time_plot_epsr(Ex/max(max(Ex)),dz,2*epsr-5,1);
%time_plot(H,dz,'testing123');



