function [moving_time]= measure_peak_moving_time(z, H,prasent, direction)
%This fuction measures the time wiche present off the peak move in z cordinton.
% z - doubel off cordintion

moving_time= zeros(size(z));
H_power=  squeeze((H(1,:,:)).^2 + (H(2,:,:)).^2 + (H(3,:,:)).^2);
H_peak= max(H_power');

for i = 1:1:length(z)
    moving_time(i) = find(H_power(z(i),:) >= prasent* H_peak(z(i)),1, direction);
end
