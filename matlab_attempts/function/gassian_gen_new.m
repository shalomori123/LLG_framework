function [Pulse, dt, lada0] = gassian_gen_new( max_fealed, bwhp, fea, NameValueArgs) 

    arguments
            max_fealed (1,1) double % max magntic filled of pouls [A/m]
            bwhp (1,1) double % the total time of power {= 1/heta*E^2} is biger the 0.5*max_ffeal 
            fea (1,1) double % phase differention of the pulse [rad]
            NameValueArgs.relation_factor = 5 % time add to pulse after pike [in sigma]
            NameValueArgs.ansition_factor = 5 % time add to pulse before pike [in sigma]
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
    pulse_width = 0.5*sqrt(1/log(2))*bwhp;
    pike_loction = NameValueArgs.ansition_factor* pulse_width;
    total_time = pike_loction + NameValueArgs.relation_factor* pulse_width;
    t = 0:dt:total_time;

    
    Pulse =max_fealed* exp(- (t- pike_loction).^2/(2*pulse_width^2));
    Pulse= amplified(W,dt,Pulse,fea, NameValueArgs.teta) ;
        

    %plot(Pulse)
    Pulse= [Pulse(:,1)';Pulse(:,2)'; zeros(1,length(t))];
    


end