function [Pulse]= amplified(w,dt, envelope, fie, teta)
% This function gets the envelope of pulse and returns the amplified pulse
% with polarization
% w = frequency
% fie = the phase differnce between ampified x and y
% teta = turning angle of main polarization


turn= [ cos(teta) , -sin(teta); sin(teta) cos(teta)];
Tsteps= length(envelope);
t= 0:dt:dt*(Tsteps-1);

Pulse= [envelope; envelope] .*[exp(1i*w*t); exp(1i*(w*t-fie))];
Pulse= turn * Pulse;
Pulse= Pulse.';