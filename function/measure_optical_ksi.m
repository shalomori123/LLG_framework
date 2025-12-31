function [ksi]= measure_optical_ksi(H,ms,W,dt,pulse_tow,moment_num,z1)
%works right now with single moment
    z_m= z1+moment_num;
    t_pik= find(H(1,z_m,:) == max(H(1,z_m,:)),1,"first"); 
    t_pess= ceil(t_pik + 2.35*pulse_tow*2*pi/W/dt);
    
    ms_fft= fft(ms(1,moment_num,t_pik:t_pess));
    ms_fft=ms_fft/length(ms_fft);
    fs= 1/dt;
    f = (0:length(ms_fft)-1)*fs/length(ms_fft);
    f_optiy= W/2/pi;
    k_op= find(f < f_optiy,1,"last");
    ms_op= abs(ms_fft(1,:,k_op))+abs(ms_fft(1,:,k_op+1));
    ksi= ms_op/max(max(max(H)));
end