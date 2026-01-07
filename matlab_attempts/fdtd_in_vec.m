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
Tsteps=1500;



%% grid and setup


dt=(1/40)*(lada0/c); % the recomended sampeling freq is ten point per wave length.
dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c

% the time spece defintion
L=17*lada0;
z=0:dz:L;
lenz=length(z);
H=zeros(Tsteps,lenz,2);
E=zeros(Tsteps,lenz,2);


% the matial defniton
z1=ceil(0.4*lenz); %location of material start
z2=ceil(0.7*lenz); %material end
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)];


e_r= [1;1]; %define eps_r
epsr= ones(2,lenz)+(e_r-1)*mat_vec;

sigma= 4000; %define condatory 
sigma= sigma*mat_vec;
eaf= dt/2/eps*sigma./epsr;
conset_for_les_E= ((1-eaf)./(1+eaf));
conset_for_H= (0.5./epsr./(1+eaf));



t=0:dt:dt*(Tsteps);
Pulse= 10^(S/10)*rectangularPulse(20,420,0:1:Tsteps);
% Pulse =create_gassian_envelope(Tsteps,dt,pulse_width,S,0.2);
lim= max(Pulse);
% Pulse= [Pulse.' Pulse.';];

Pulse= amplified(0,dt,Pulse,0, 0);
plot3(t,Pulse(:,1), Pulse(:,2));
axis([ 0 t(end)   -lim lim -lim lim])
%% propagate

for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
    %n
    E(n,4,1)=E(n,4,1)+Pulse(n,1); %inject the input signal
    E(n,4,2)=E(n,4,2)+Pulse(n,2);

    H(n+1,2:lenz,1)=H(n,2:end,1)+0.5*(E(n,2:end,2)-E(n,1:lenz-1,2));
    E(n+1,1:lenz-1,2)=conset_for_les_E(2,1:lenz-1).*E(n,1:lenz-1,2)+conset_for_H(2,1:lenz-1).*(H(n+1,2:lenz,1)-H(n+1,1:lenz-1,1));
        
    H(n+1,2:lenz,2)=H(n,2:lenz,2)-0.5*(E(n,2:lenz,1)-E(n,1:lenz-1,1));
    E(n+1,1:lenz-1,1)=conset_for_les_E(1,1:lenz-1).*E(n,1:lenz-1,1)-conset_for_H(1,1:lenz-1).*(H(n+1,2:lenz,2)-H(n+1,1:lenz-1,2));
    
%     for i=1:1:lenz-1
%         %i
%         H(n+1,i+1,1)=H(n,i+1,1)+0.5*(E(n,i+1,2)-E(n,i,2));
%         E(n+1,i,2)=conset_for_les_E(i)*E(n,i,2)+conset_for_H(i)*(H(n+1,i+1,1)-H(n+1,i,1));
%         
%         H(n+1,i+1,2)=H(n,i+1,2)-0.5*(E(n,i+1,1)-E(n,i,1));
%         E(n+1,i,1)=conset_for_les_E(i)*E(n,i,1)-conset_for_H(i)*(H(n+1,i+1,2)-H(n+1,i,2));
%         
%     end
end

for n=3:1:Tsteps-1
    %n
   E(n,4,1)=E(n,4,1)+Pulse(n,1); %inject the input signal
   E(n,4,2)=E(n,4,2)+Pulse(n,2);
   E(n,lenz,:)=E(n-2,lenz-1,:);%right boundery condition
   E(n,1,:)=E(n-2,2,:);%left boundery condition

   % instnd i loope

    H(n+1,2:lenz,1)=H(n,2:end,1)+0.5*(E(n,2:end,2)-E(n,1:lenz-1,2));
    E(n+1,1:lenz-1,2)=conset_for_les_E(2,1:lenz-1).*E(n,1:lenz-1,2)+conset_for_H(2,1:lenz-1).*(H(n+1,2:lenz,1)-H(n+1,1:lenz-1,1));
        
    H(n+1,2:lenz,2)=H(n,2:lenz,2)-0.5*(E(n,2:lenz,1)-E(n,1:lenz-1,1));
    E(n+1,1:lenz-1,1)=conset_for_les_E(1,1:lenz-1).*E(n,1:lenz-1,1)-conset_for_H(1,1:lenz-1).*(H(n+1,2:lenz,2)-H(n+1,1:lenz-1,2));
    
%     for i=1:1:lenz-1
%         %i
%         
%         H(n+1,i+1,1)=H(n,i+1,1)+0.5*(E(n,i+1,2)-E(n,i,2));
%         E(n+1,i,2)=conset_for_les_E(i)*E(n,i,2)+conset_for_H(i)*(H(n+1,i+1,1)-H(n+1,i,1));
%         H(n+1,i+1,2)=H(n,i+1,2)-0.5*(E(n,i+1,1)-E(n,i,1));
%         E(n+1,i,1)=conset_for_les_E(i)*E(n,i,1)-conset_for_H(i)*(H(n+1,i+1,2)-H(n+1,i,2));
%         
%     end
end

%%


% time_plot_2D_fields(E/max(max(max(E))),dz,z(z1),z(z2),2);
time_plot_2_fields(E(:,:,1), E(:,:,2),dz,mat_vec*2.02-1.01,1,"E_x","E_y",'rect pulse sigma=0 epsilon=0 w=0');



