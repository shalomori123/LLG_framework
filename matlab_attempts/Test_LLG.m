clear all; close all ; clc;
%%
eps=8.854e-12;
miu=4*pi*1e-7;
eta=sqrt(miu/eps);
c=1/sqrt(eps*miu);
%%
lada0=800e-9;
W=2*pi*c/lada0;


Tsteps=200000;
dt= 4*2.5e-18;
H0= 4e8;
H= [zeros(2,Tsteps); H0*ones(1,Tsteps)];
M0= 0.3e6;
M= [  M0*ones(1,4,Tsteps); zeros(2,4,Tsteps) ];
T_Global= 0:dt:(Tsteps-1)*dt;

Pulse= amplified(0,dt,H0*ones(1,Tsteps),pi/2, pi/2);
% Pulse= [Pulse ;Pulse];

Pulse_1 = zeros(size(M));
for i= 1:1:4
    Pulse_1(:,i,:)= [Pulse  zeros(Tsteps,1)].';
end
plot3(T_Global ,squeeze(Pulse_1(1,2,:)), squeeze(Pulse_1(2,2,:)));
axis([ 0 T_Global(end)   -H0 H0 -H0 H0])

%%
for n = 1:1:Tsteps-1
    M(:,:,n+1)= runge_kutte_LLG_1(M(:,:,n),Pulse_1(:,:,n),dt,20*0.0125);
end

% H= H/H0;
% M= M/ max(max(M));
% plot3([0 M(1,1)],[0 M(2,1)], [0 M(3,1)])
% hold on
% plot3([0 M(1,end)],[0 M(2,end)], [0 M(3,end)])
% 
% plot3([0 H(1,1)],[0 H(2,1)], [0 H(3,1)])

%%
plot_moments(M,dt,4)
