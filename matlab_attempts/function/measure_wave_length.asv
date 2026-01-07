function [wave_length] = measure_wave_length(field,z_lim)
    [pks1,locs1] = findpeaks(field(1,z_lim(1):z_lim(2)),z_lim(1):z_lim(2));
    [pks2,locs2] = findpeaks(field(2,z_lim(1):z_lim(2)),z_lim(1):z_lim(2));
    
    wave_length_vec= [locs1(2:end)- locs1(1:end-1) locs2(2:end)- locs2(1:end-1)];
    wave_length= mean(wave_length_vec);
    if (length(wave_length_vec) == 0)
        wave_length =0 ;
    end

end