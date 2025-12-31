function plot_polarization(H_set,Hext_set, dt,W,z1,z2)
    
    T= ceil(2*pi/W/dt);
    Time_us= (length(H_set{1})-T):length(H_set{1}) - 1;


    nexttile;
    hold on
    grid on

    for ind= 1: length(H_set)
        H=H_set{ind};
        H_x= squeeze(H(1,z2+10,Time_us));
        H_y =squeeze(H(2,z2+10,Time_us));
        plot(H_x,H_y,"DisplayName","H_{ext}= "+ num2str(Hext_set(ind))+ "_T")
    end
    plot(squeeze(H(1,z1-10,Time_us)),squeeze(H(2,z1-10,Time_us)),"DisplayName","befor matral")


 




end