function [Pulse, dt, lada0] = gassian_gen( power, fea, NameValueArgs) 

    arguments
            power (1,1) double % toatl power in the pulse [J] ,p= heta/2 *E^2 bacause normalization of E
            fea (1,1) double % phase differention of the pulse [rad]
            NameValueArgs.pulse_width = 2.6685e-14  % time to get lowe to e^-1 of max [sec], default 10 cycal in defult lamda
            NameValueArgs.relation_time = 10*2.6685e-14 % time add to pulse after low e^-2.5 of max [sec], default 5 cycal in defult lamda
            NameValueArgs.time_cut=40 % resulution parameter of how much to cut T/dt
            NameValueArgs.lamd0= 800e-9 % wavelength [m]
            NameValueArgs.teta= 0 % initial angle of the polarization [rad]
    end
    
   % phisical constants;
    eps=8.854e-12; %[F/m]
    miu=4*pi*1e-7; %[H/m]
    c=1/sqrt(eps*miu);
    heta = sqrt(miu/eps);
    
    % sistem constents 
    
    lada0= NameValueArgs.lamd0;
    dt=(1/NameValueArgs.time_cut)*(lada0/c); % the recomended sampeling freq is ten point per wave length is (1/40)*(lada0/c);
    
    W=2*pi*c/lada0;
    pulse_width = NameValueArgs.pulse_width;
    relation_time = NameValueArgs.relation_time;

    Tsteps= ceil((relation_time + 5 * pulse_width)/dt) ;
    peak_location =2.5*pulse_width/dt/Tsteps;
    
    Pulse =create_gassian_envelope(Tsteps,dt,pulse_width,1e7,peak_location);
    Pulse= amplified(W,dt,Pulse,fea, NameValueArgs.teta) ;
    

    norm = sum(sum(abs(Pulse).^2)); % nead power fctor!!! 2 hate
    Pulse = sqrt(1/(norm))*Pulse;
    Pulse = sqrt(2*power/(heta))*Pulse;
    

    %plot(Pulse)
    Pulse= [Pulse(:,1)';Pulse(:,2)'; zeros(1,Tsteps+1)];
    


end