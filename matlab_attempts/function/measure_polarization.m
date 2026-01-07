function [fie_mean, theta_mean, fie_steady , theta_steady] = measure_polarization(field,z_loc)
    % return :
    % theta - the angle between the polarization direction and the x axis
    % fie - the phase difference between x and y
    x=squeeze(field(1,z_loc,:));
    y= squeeze(field(2,z_loc,:));

    [pks_x,locs_x] = findpeaks(x);
    [pks_,locs_y] = findpeaks(y);
    
    T_vec= [locs_x(2:end)- locs_x(1:end-1) ; locs_y(2:end)- locs_y(1:end-1)];
    T= mean(T_vec);
    max_ind= min(length(locs_y), length(locs_x));
    fie_vec= locs_x(1:max_ind)- locs_y(1:max_ind);
    fie_mean= 2*pi*mean(fie_vec)/T;
    fie_steady= 2*pi*(locs_x(end)- locs_y(end))/T;
    if fie_steady< -pi
        fie_steady= fie_steady+ 2*pi;
    elseif fie_steady> pi
        fie_steady= fie_steady- 2*pi;
    end

    if fie_mean< -pi
        fie_mean= fie_mean+ 2*pi;
    elseif fie_mean> pi
        fie_mean= fie_mean- 2*pi;
    end

    x_mean= mean(x(y>0));
    y_mean= mean(y(x>0));
    theta_mean= atan(y_mean/x_mean);
    x_last_cycle= x(end-ceil(T):end);
    y_last_cycle= y(end-ceil(T):end);

    theta_steady= atan(mean(y_last_cycle(x_last_cycle>0))/mean(x_last_cycle(y_last_cycle>0)));
end