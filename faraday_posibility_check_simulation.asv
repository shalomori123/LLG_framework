clear all, clc
eps=8.854e-12;
miu=4*pi*1e-7;
c=1/sqrt(eps*miu);

q=1.60217646e-19;                        %[cb]
g=2;                                     %[Landau factor]
me=9.1093821545*1e-31;                   %[gram] electron mass
gma_factor=1;
gma = gma_factor*g*q/(2*me);  
lamda= 800e-9;


%%
[clk_s tf] = clock;
e_r= [1;1]; %define eps_r
% max_field= 5e9;  % [A/m]
% pulse_width= 10;
sigma= [ 0 ;  0]*1e-3; %define condatory  [S/m]
% ms0= [ 0 3e5 0 ];
% alpha= 20*0.00125;
% close_system = false;
Teslta= 795774.7; %[A/m] = 1[T]
Bex_res= c/lamda/gma;

run_amplitude= logspace(4, 6,3);
run_ms= linspace(1e3,1e5,6);
run_alpha= linspace(0.05,0.2,3);
run_fea= [pi/2 -pi/2 0 ];
run_external_H= [1:7:22 Bex_res]; % [T]

polarstion_name = [ "left" , "right","linar"];

ms_size_ratio=         zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H),length(run_amplitude));
stedy_lamda_in =       zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));
stedy_lamda_out =      zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));

post_material_tau_stedy=   zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));
post_material_epsilon_stedy=     zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));
pro_material_tau_stedy=    zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));
pro_material_epsilon_stedy=      zeros(length(run_fea),length(run_ms), length(run_alpha), length(run_external_H) ,length(run_amplitude));



path= 'C:\maxwell-LLG\ws\simulation_result\posibilty_faradey_strong_H_ext\';
%%
sm_nm=1;

space_resolution= 1;
time_to_save= 800;
lamda= 800e-9;

for H_ext= run_external_H
    for ind_amp= 1:length(run_amplitude)    
        for alpha= run_alpha
            for ms0= run_ms
                for index_fie= 1:3
                    if index_fie ==3
                        amplitude= run_amplitude(ind_amp)/sqrt(2);
                    else
                        amplitude= run_amplitude(ind_amp);
                    end                   
                     file_name= path + "faradey_amp_"+ num2str(run_amplitude(ind_amp))+ "_Hext_" + num2str(H_ext)   + "_a_" + num2str(alpha) + "_ms0_" + num2str(ms0) + "_polriztion_"+ polarstion_name(index_fie)  + ".mat";
    %                  if


                        [E,H,ms, z1, z2, z_len,dz,W, Tsteps] =faraday_simulation_complex(e_r,sigma,ms0,alpha, amplitude, run_fea(index_fie), H_ext*Teslta , 1 ...
                            ,"time_cut",40*space_resolution,"start_matiral",0.1,"len_matiral",12*lamda,"len_in_lemda",14,"simlation_time",5e-11 ...
                            ,"lamd0", lamda);

                        lamda_in= measure_wave_length(real(E(:,:,end)),[z1 z2])*dz;
                        lamda_out= measure_wave_length(real(E(:,:,end)),[(z2+1) z_len])*dz;
                        
                        
                        ms_start_size= sum(sqrt(ms(1,:,1).^2 + ms(2,:,1).^2 + ms(3,:,1).^2));
                        ms_finish_size= sum(sqrt(ms(1,:,end).^2 + ms(2,:,end).^2 + ms(3,:,end).^2));
                        ms_ratio= abs(ms_finish_size - ms_start_size)/ ms_start_size ;
                        
                        [tau_pro_stedy  ,epsilon_pro_stedy] = polellip([E(1,z1-1,end);E(2,z1-1,end)]);
                        [tau_post_stedy , epsilon_post_stedy] = polellip([E(1,z2+1,end);E(2,z2+1,end)]);



                        ms_size_ratio(sm_nm)= ms_ratio;
                        stedy_lamda_in(sm_nm)= lamda_in;
                        stedy_lamda_out(sm_nm)= lamda_out;

                     
                        post_material_epsilon_stedy(sm_nm)= epsilon_post_stedy;
                        pro_material_epsilon_stedy(sm_nm)= epsilon_pro_stedy;
                        pro_material_tau_stedy(sm_nm)= tau_pro_stedy;
                        post_material_tau_stedy(sm_nm)= tau_post_stedy;

                        sm_nm= sm_nm+1;
                        
                        if ms_ratio<1
                            H= H(:,1:space_resolution:end,(end-time_to_save):end);
                            E= E(:,1:space_resolution:end,(end-time_to_save):end);

                            S= cross(E,H);
                            S= squeeze(S(3,:,:));
                            ms= ms(:,1:space_resolution:end,(end-time_to_save):end)*space_resolution;
                           
                                
                            save(file_name, 'H', 'ms',"S", "ms_ratio")
                        end
                        if sm_nm==2
                            [clk_f tf] = clock;
                            dz_save_file= dz*space_resolution;
                            dt= dz/2/c;
                           
                            simlation_time= dt*Tsteps;
                            
                            file_name= [path 'safety_for_all_running_for_date_' num2str(clk_f(1)) '_' num2str(clk_f(2)) '_' num2str(clk_f(3)) '.mat' ];
                           save(file_name, 'e_r' , 'sigma', "run_ms","run_fea", "run_alpha", "run_external_H","run_amplitude",...
                                "z1", "z2", "z_len", 'dz', "dt","dz_save_file", "space_resolution", 'W', "simlation_time") 
                        end
    %                  catch
    %                     Error = ['fail in runing file:  ' file_name]
    %                  end
                    
                end
            end
        end
    end
end
[clk_f tf] = clock;
dz_save_file= dz*space_resolution;
dt= dz/2/c;
dt_save_file= dt*time_to_save;
simlation_time= dt*Tsteps;

file_name= [path 'const_for_all_running_for_date_' num2str(clk_f(1)) '_' num2str(clk_f(2)) '_' num2str(clk_f(3)) '.mat' ];
save(file_name, 'e_r' , 'sigma',"run_ms","run_fea", "run_alpha", "run_external_H","polarstion_name", "run_amplitude", ...
     "z1", "z2", "z_len", 'dz', "dt", "dz_save_file", "space_resolution", 'W',"clk_s", "clk_f", "simlation_time", ...
    "stedy_lamda_in","stedy_lamda_out", "ms_size_ratio", "pro_material_theta_mean","pro_material_fie_mean","post_material_theta_mean","post_material_fie_mean", ...
     "pro_material_tau_stedy","pro_material_epsilon_stedy","post_material_tau_stedy","post_material_epsilon_stedy") 


%%
d =sum(sum(sum(kiping_ms_size,5),3))
% size(simltion_statos);


ind_mf=5;
ind_al=5;

plot(run_ms,lemad_close(ind_al,:,ind_mf),run_ms,lemad_open(ind_al,:,ind_mf))
%%

z=0:dz:z_len*dz;
mat_vec= [zeros(1,z1) ones(1,z2-z1) zeros(1,z_len-z2)];
time_plot_2_fields(squeeze(H(1,:,:))',squeeze(H(2,:,:))',dz,mat_vec*2.02-1.01,10,"H_x open","H_x close");


% time_plot_2D_fields(H,dz,z(z1), z(z2),2,'')

%%

plot((0:(z_len-1))*dz,E(1,:,end))