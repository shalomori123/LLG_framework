function [moving_time]= measure_power_moving_time(z1, z2, H,dt,power_persent)
%This fuction measures how much time the peak takes to move from z1 to z2.


H_power= squeeze(cumsum( (H(1,[z1 z2],:)).^2 + (H(2,[z1 z2],:)).^2 + (H(3,[z1 z2],:)).^2));
M=max(H_power');




t1 = find(H_power(1,:) > power_persent*M(1),1);
t2=  find(H_power(2,:) > power_persent*M(2),1);
moving_time= t2- t1;