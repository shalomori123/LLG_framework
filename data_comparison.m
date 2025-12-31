clear all 
clc

% phisical constants;
eps=8.854e-12; %[F/m]
miu=4*pi*1e-7; %[H/m]
c=1/sqrt(eps*miu);
Teslta= 795774.7; %[A/m] = 1[T]
%%

path= 'C:\maxwell-LLG\ws\simulation_result\posibilty_invr\';
conset_file= 'const_for_all_running_for_date_2024_1_11.mat';

conset= open([path conset_file]);

polarstion_name = ["right", "left" ,"linar"];




%% loding lamda averg data

lamda_open= conset.lamda_open; %[alpha]* [ms] * [max_fild]
lamda_close= conset.lamda_close;
ms_value= conset.ms_run; % [fea]*[ms]*[alpha]*[H_ext]*[amplitou]

ms_valid= conset.ms_size_ratio<1e-3;
max_ind_mf= length(conset.max_fild_run);

max_ind_alp= length(conset.alpha_run);

%%
for ind_mf = 1: max_ind_mf
    figure
    hold on
    grid on
    for ind_a =1:max_ind_alp            
        plot(ms_value, squeeze(lamda_close(ind_a,:,ind_mf)),DisplayName="alpha =" + num2str(conset.alpha_run(ind_a)))
    end
    title( "max field: " + num2str( floor(conset.max_fild_run(ind_mf)/1e9))+ "_{GA/m}")
    xlabel('M_s [A/m]')
    ylabel('wave len [m]')
    hold off
end



%%  compering the alphe gref
ind_amp=1;

for ind_Hext= 1: max_ind_Hext 
    figure
    hold on
    grid on
    tiledlayout(3,1)
    for p=1:3
        nexttile
        for ind_alpha= 1:max_ind_alp
            name =  " a: "+ num2str(conset.run_alpha(ind_alpha)) ;
            valid= ms_valid(p,:,ind_alpha, ind_Hext,ind_amp);
            plot(ms_value(valid),lamda_open(p,valid,ind_alpha,ind_Hext,ind_amp), 'DisplayName',name )
            
        end
%         plot(ms_value(valid), lamda(ind_alpha,valid,ind_Hext,ind_mf),'k*')
        
    end
    title(polarstion_name(p) + "  H external: " + num2str(conset.run_extrnal_H(ind_Hext)) + "_T" +   " amplitod: " + num2str( floor(conset.run_amplitou(ind_amp)/1e6))+ "_{MA/m}")
    xlabel('M_s [A/m]')
    ylabel('wave len [m]')
    hold off
end

%% compering the H esternal grap

for ind_alpha= 1:length(conset.alpha_run) 
    figure
    hold on
    grid on
    for ind_Hext= 1: length(H_ext)
        name =   "H_0 = " + num2str(H_ext(ind_Hext)/Teslta) + "_T";
        valid= ms_valid(ind_alpha,:, ind_Hext,ind_amp);
        plot(ms_value(valid),lamda(ind_alpha,valid,ind_Hext,ind_amp), 'DisplayName',name )
        
%     
%         plot(ms_value(valid), lamda(ind_alpha,valid,ind_Hext,ind_mf),'k*')
        
    end
    title("alpha: "+ num2str(conset.alpha_run(ind_alpha)) +  " max field: " + num2str( floor(conset.max_fild_run(ind_amp)/1e6))+ "_{MA/m}" )
    xlabel('M_s [A/m]')
    ylabel('wave len [m]')
    hold off
end
%% compering the max power of the light


ind_Hext= 1;

for ind_alpha= 1:length(conset.alpha_run) 
    figure
    hold on
    grid on
    for ind_amp= 1: length(conset.max_fild_run)
        name =  " max field: " + num2str( floor(conset.max_fild_run(ind_amp)/1e6))+ "_{MA/m}";
        valid= ms_valid(ind_alpha,:, ind_Hext,ind_amp);
        plot(ms_value(valid),lamda(ind_alpha,valid,ind_Hext,ind_amp), 'DisplayName',name )
        
%     
%         plot(ms_value(valid), lamda(ind_alpha,valid,ind_Hext,ind_mf),'k*')
        
    end
    title("alpha: "+ num2str(conset.alpha_run(ind_alpha)) +   "  H_0 = " + num2str(H_ext(ind_Hext)/Teslta) + "_T" )
    xlabel('M_s [A/m]')
    ylabel('wave len [m]')
    hold off
end
%% lodind and runing time plot comretion



ind_ms= 4;
ind_alp= 4;
ind_amp= 1;

file_name_mf_w_a_open= [path 'gasuian_mf_' num2str(conset.max_fild_run(ind_amp)) '_w_' num2str(conset.pulse_width_run) '_a_' num2str(conset.alpha_run(ind_alp)) '_ms0_' num2str(conset.ms_run(ind_ms)) '_close_0.mat'];
file_name_mf_w_a_close= [path 'gasuian_mf_' num2str(conset.max_fild_run(ind_amp)) '_w_' num2str(conset.pulse_width_run) '_a_' num2str(conset.alpha_run(ind_alp)) '_ms0_' num2str(conset.ms_run(ind_ms)) '_close_1.mat'];


simlation_op= open(file_name_mf_w_a_open);
simlation_cl= open(file_name_mf_w_a_close);


% file_name_c= [path 'gasuian_mf_' num2str(conset.max_fild_run(ind_mf)) '_w_' num2str(conset.pulse_width_run) '_a_' num2str(conset.alpha_run(ind_alp)) '_ms0_' num2str(conset.ms_run(ind_ms)) '_close_' char('1') '.mat'];
% 
% data_c= open(file_name_c);

%%
len_z_use= size(simlation_cl.S);
z1= floor(len_z_use(1)*conset.z1/conset.z_len);
z2= ceil(len_z_use(1)*conset.z2/conset.z_len);
dz= conset.dz;
save_name = 'C:\maxwell-LLG\videos\inverst farday\compare_Hx_open_close';

time_plot_2_fields(squeeze(simlation_op.H(1,:,:))',squeeze(simlation_cl.H(1,:,:))',dz,z1, z2,3,"H_x open system","H_x close system",save_name);

%% plot momnt
dt= conset.dt*conset.resultion;
plot_moments(simlation1.ms,dt,5)

Ms= squeeze(sum(simlation1.ms,2));
T= length(Ms);

T_Global= linspace(0,dt*(T-1),T);

figure
plot(T_Global,squeeze(Ms(1,:)),'linewidth',2,'DisplayName','M_x');
hold on
plot(T_Global,squeeze(Ms(2,:)),'linewidth',2,'DisplayName','M_y');
plot(T_Global,squeeze(Ms(3,:)),'linewidth',4,'DisplayName','M_z')
plot(T_Global,squeeze(((Ms(1,:)).^2 + (Ms(2,:)).^2 + (Ms(3,:)).^2).^0.5),'linewidth',1.5,'DisplayName','Magnetization Absolute Value');

title(['all Ms , Ms size chaeng in: 10^{' num2str(floor(log10(simlation1.ms_ratio)+1)) '} of origanl Value' ])
xlabel('time [sec]'), ylabel('magentizasion value [A/m]')

%%

% plot_dir(squeeze(H(1,z1-10,Time_us)),squeeze(H(2,z1-10,Time_us)))

%%
lamda= squeeze(conset.evr_lamda(2,:,:,1,:));
lamda_open= squeeze(conset.evr_lamda(1,:,:,1,:));

ind_Hext= 1 ;
ind_alpha= 3 ;
ind_amp= 4;
hold on
grid on
plot(conset.ms_run,lamda_open(ind_alpha,:,ind_amp))
plot(conset.ms_run,lamda(ind_alpha,:,ind_amp))









