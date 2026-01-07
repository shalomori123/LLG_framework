
eps=8.854e-12;
miu=4*pi*1e-7;
c=1/sqrt(eps*miu);

%%
[clk_s tf] = clock;
e_r= [5;5]; %define eps_r
% max_field= 5e9;  % [A/m]
% pulse_width= 10;
sigma= [ 5.8e7 ;  5.8e7]*1e-3; %define condatory  [S/m]
% ms0= [ 0 3e5 0 ];
% alpha= 20*0.00125;
fea = pi/2;
% close_system = false;
max_fild_run= logspace(7,11,4);
pulse_width_run= 9:12;
ms_run= 1e5:7e5:4e6;
alpha_run= linspace(0.00125,40*0.00125,5);

simltion_statos= zeros(length(alpha_run),length(ms_run),length(pulse_width_run),length(max_fild_run));

path= 'C:\maxwell-LLG\ws\simulation_result\sigma_58000_er_5\';
%%
sm_nm=1;
for max_field= max_fild_run
    for pulse_width= pulse_width_run
        for ms0= ms_run
            for alpha= alpha_run
                for close_system= [false true]
                     file_name= [path 'gasuian_mf_' num2str(max_field) '_w_' num2str(pulse_width) '_a_' num2str(alpha) '_ms0_' num2str(ms0) '_close_' char(close_system+'0') '.mat'];
                     try
                        [E,H,ms, z1, z2, z_len,dz,W] =gassian_simulation(e_r,sigma,ms0,alpha, ...
                            max_field, pulse_width, fea, close_system ...
                            ,"time_cut",800,"len_matiral",400e-9,"len_in_lemda",1.5);
                        
                        %
                        S= cross(E,H);
                        
                        % measere the optical linar ksi
                        ksi= zeros(1,z2- z1);
                        
                        for k= 1:z2-z1
                            ksi(k)= measure_optical_ksi(H,ms,W,dz/c/2,pulse_width,k,z1);
                        end
                        
                        % measere the wavelen and velocity 
                        t_in= find(S(3,z2-z1,:)> 0.3 * max(max(S(3,1:z1,:))),1 ,'first'); %entring to matiral time index
                        t_out= find(S(3,z2,:)> 0.3 * max(max(S(3,z2:end,:))),1 ,'last'); %exit matiral time index
                        lamda_void= zeros(1,t_out-t_in);
                        lamda_in= zeros(1,t_out-t_in);
                        for t = 1: t_out-t_in 
                            lamda_void(t)= measure_wave_length(H(:,:,t),[z2-z1-1 (z1-1)])*dz;
                            lamda_in(t)= measure_wave_length(H(:,:,t),[z1 z2])*dz;
                        end                      
                        save(file_name, 'E' ,'H', 'ms', 'ksi' ,'lamda_in',"lamda_void","S", 'z_len', 'z1', 'z2', "dz", "t_in","t_out")
                        simltion_statos(sm_nm)=1;
                        sm_nm= sm_nm+1;

                     catch
                        Error = ['fail in runing file:  ' file_name]
                     end
                end
            end
        end
    end
end


[clk_f tf] = clock;

file_name= [path 'const_for_all_running_for_date_' num2str(clk_f(1)) '_' num2str(clk_f(2)) '_' num2str(clk_f(3)) '.mat' ];
save(file_name, 'e_r' , 'sigma', "z1", "z2", "z_len", 'dz', 'W',"clk_f","clk_s", "max_fild_run", "ms_run",'alpha_run',"pulse_width_run") 


