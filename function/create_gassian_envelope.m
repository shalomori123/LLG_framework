function [pulse]=create_gassian_envelope(steps,dt,tou_p,max_felds,center)
%create a gaussian pulse centered around meu and width of sig.
%create_pulse(steps,dt,W,FWHM,E)
%E is the energy [J/m^2] of pulse in dB units. 

eps=8.854e-12;
miu=4*pi*1e-7;
% eta=sqrt(miu/eps);              %medium impedance

% max_felds=10^(max_felds/10);                    %energy in J

%t=0:dt:dt*(steps-1);
t=0:dt:dt*(steps);
meu=t(floor(steps*center));
%FWHM= 2.35*tou_p;
sig = tou_p;
% A=max_felds/(eta*sqrt(pi)*sig);
pulse=max_felds*exp(-((t-meu).^2)/(2*sig^2));
% plot(t,pulse)
