clear all, clc
eps=8.854e-12;
miu=4*pi*1e-7;
c=1/sqrt(eps*miu);

%%
% [clk_s tf] = clock;
e_r= [1;1]; %define eps_r
% max_field= 5e9;  % [A/m]
% pulse_width= 10;
sigma= [ 0 ;  0]*1e-3; %define condatory  [S/m]
% ms0= [ 0 3e5 0 ];
% alpha= 20*0.00125;
fea = pi/2;
% close_system = false;
max_fild_run= linspace(1e9,1e11,3);
pulse_width= 10;
ms_run= linspace(3e10,1e12,6);
alpha_run= linspace(0.02,0.2,5);
Teslta= 795774.7; %[A/m] = 1[T]
H_ext= linspace(1,3,4); % in [T]
ms_size_ratio= zeros(length(alpha_run),length(ms_run),length(H_ext),length(max_fild_run));
evr_lamda = zeros(length(alpha_run),length(ms_run),length(H_ext),length(max_fild_run));
path= 'C:\maxwell-LLG\ws\simulation_result\gaussian_faraday\';
%%
resultion= 4;

sm_nm=1;
for max_field= max_fild_run
    for H_i= H_ext
        H_z= [0 ,0 ,H_i*Teslta];
        for ms0= ms_run
            ms_u= [  0 ,0 , ms0 ];
            
            for alpha= alpha_run
                
                     file_name= [path 'gasuian_mf_' num2str(max_field) '_Hext_' num2str(H_i) '_a_' num2str(alpha) '_ms0_' num2str(ms0) '.mat'];
                   
                    [E,H,ms, z1, z2, z_len,dz,W ,Tsteps] =gassian_simulation(e_r,sigma,ms_u,alpha, max_field, pulse_width, fea, H_z, 1 ...
                        ,"time_cut",80*resultion ,"len_matiral",1800e-9,"len_in_lemda",3);
                    
                    %
                  
                    S= cross(E,H);
                    
                    
                    % measere the optical linar ksi
%                         ksi= zeros(1,z2- z1);
%                         
%                         for k= 1:z2-z1
%                             ksi(k)= measure_optical_ksi(H,ms,W,dz/c/2,pulse_width,k,z1);
%                         end
%                         
                    % measere the wavelen and velocity 

                    t_in= find(S(3,z1,:)> 0.3 * max(max(S(3,1:z1,:))),1 ,'first'); %entring to matiral time index
                    t_out= find(S(3,z2,:)> 0.3 * max(max(S(3,z2:end,:))),1 ,'last'); %exit matiral time index
                    lamda_in= zeros(1,t_out-t_in);
                    for t = t_in : t_out 
                        l= measure_wave_length(H(:,:,t),[z1 z2])*dz;
                        if  ~(isnan(l))
                            lamda_in(t)= l;
                        end
                    end
                    x= lamda_in> 0 ;
                    lamda= mean(lamda_in(x));
                   

                    ms_finsh_s = sqrt(ms(1,:,end).^2 + ms(2,:,end).^2 +ms (3,:,end).^2);
                    ms_start = sqrt(ms(1,:,1).^2 + ms(2,:,1).^2 +ms (3,:,1).^2);
                    ms_ratio = abs(sum(ms_finsh_s)-sum(ms_start))/sum(ms_start);
                    ms_size_ratio(sm_nm)=  ms_ratio;
                    evr_lamda(sm_nm)= lamda;
                    sm_nm= sm_nm+1;


                    H= H(:,1:resultion:end,1:resultion:end);
                    S= squeeze(S(3, 1:resultion:end, 1:resultion:end));
                    lamda_in= lamda_in(1: resultion :end);
                    ms= ms(:,:,1:resultion:end);
                      try
                        save(file_name, 'H', 'ms','lamda_in',"S", "t_in","t_out", 'ms_ratio')
                        if (sm_nm==2)
                            dt= dz/2/c;
                            run_time= dt*Tsteps;
                            dz_saved_file= dz*resultion;
                            [clk_f tf] = clock;
                            
                            file_name= [path 'const_for_all_running_for_date_' num2str(clk_f(1)) '_' num2str(clk_f(2)) '_' num2str(clk_f(3)) '.mat' ];
                            save(file_name, 'e_r' , 'sigma', "z1", "z2", "z_len", 'dt', 'dz_saved_file', 'W',"clk_s", 'run_time', 'H_ext', ...
                                "max_fild_run", "ms_run",'alpha_run',"pulse_width",'evr_lamda', "ms_size_ratio", 'resultion') 
                        end
                     catch
                        Error = ['fail saving file:  ' file_name]
                     end
                
            end
        end
    end
end

dt= dz/2/c;
run_time= dt*Tsteps;
dz_saved_file= dz*resultion;
[clk_f tf] = clock;


file_name= [path 'const_for_all_running_for_date_' num2str(clk_f(1)) '_' num2str(clk_f(2)) '_' num2str(clk_f(3)) '.mat' ];
save(file_name, 'e_r' , 'sigma', "z1", "z2", "z_len", 'dt', 'dz_saved_file', 'W',"clk_s", 'run_time', 'H_ext', ...
    "max_fild_run", "ms_run",'alpha_run',"pulse_width",'evr_lamda', "ms_size_ratio", 'resultion') 

%%
d =sum(sum(sum(ms_size_ratio ,5),3))
% size(simltion_statos);


ind_mf=3;
ind_al=3;

plot(ms_run,ev_close(ind_al,:,ind_mf),ms_run,lamda_open(ind_al,:,ind_mf))

%%

