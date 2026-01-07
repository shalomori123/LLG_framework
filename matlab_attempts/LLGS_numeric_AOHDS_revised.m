%% Simulation - Spin Current Amplification & CP pulses in the LLG equation
% Benjamin Assouline 
%
clear all; close all force; clc;
%
close_fig=0; %this closes figures immediatly after they open.
parameter_name='B_Pulse_amp'; 
%
%% General instructions
% This code simulates the LLGS equation under various excitations.
% The parameters that can be swept are: 
% The AC and DC current (resulting in AC-STT and DC-STT),the inital phase of the RF excitation,
% the amplitude and polarization of the H_pump field,  
% the amplitude, duration and number of Gaussian pulses.
%
% The code is divided into the following sections:
% Physical constants - *** do not change ***
% General settings
% Simulation settings
% Type of sweep
% Simulation loop
% Plots
%
%% Physical constants
q=1.60217646e-19;                        %[cb]
eps=8.854e-12;                           %[cb/V*m]  = [F/m]
miu=4*pi*1e-7;                           %[H/m]
eta0=sqrt(miu/eps);                      %[Oham]
c=1/sqrt(eps*miu);                       %[m/sec]
hbar=6.62606885e-34/2/pi;                %[E*sec]
g=2;                                     %[Landau factor]
me=9.1093821545*1e-31;                   %[gram] electron mass
gma_factor=1;
gma=gma_factor*g*q/(2*me);               %[rad/sec*T]
n_id = 1;                                %refractive index of air
%
%% General settings
% Determines wether the simulation will plot anything
do_plot = 0*1;
% This parameter determines if only the last run's plots will be
% displayed, or if all plots from all the runs will be displayed
disp_only_last_plots = 1;
% To save figures at the end of the run, set save_fig = 1;
save_fig = 0;
% To save the figures in their original size, set save_original_size = 1;
save_original_size = 0;
% To save the figures in full screen size, set save_full_screen = 1;
save_full_screen = 1;
% The figures will be saved to the location (in string format)
adress_str = 'C:\Users\owner\Desktop\Spin Current Amplifier\Simulation Results\Dmag\alpha = 0.023 f = 10G Hz\Nz = -1 PMA\singles';
%
% Define time display - either ns, ps or fs (only 1 of the following can be 1).
disp_time_femto = 1;
disp_time_pico = 0;
disp_time_nano = 0;
%
%% Simulation settings
% Gilbert damping 
alpha= 4*0.0125;%0.035; % %0.0075; 
%
% CW frequency [Hz]. In STT and rf driven simulations it is in the GHz range, and
% in AOS simulations in the optical range (~374THz).
Freq_CW= 1.25*3.7474*(10^14);% 1e9;    
%
% Angular frequency
w= 2*pi*Freq_CW; 
%
% Spin Hall angle
Theta_sh=0.148;  
%
% Ferromagnet thickness [m]
Thick=11.5e-10;
%
% Magnetization in cgs is in [emu/cm3]. 1 [emu/cm3] --> 1*10^3 [A/m]
M0= 0.3e6;%0.3e6; %[A/m]  
%
% Perturbation in Magnitization for intitial conditions
dM0= 0e-06; %[A/m] 1000*1.5e-06;
%
% The anisotropy field. Multiply by 1e-9* to reach quickly steady state
B0_anis= 0*0.02; %[T]
%
% Anisotropy field constant[J/m^3]
Ku= 0.5*B0_anis*M0;%1e-9*0.6*0.5*10000;  
K0=Ku/miu;
%
% Low temporal resolution [sec]
dt_low_res= 4*4*2.5e-18;%4*0.25*1e-18;%0.005*1e-18;% 
%
% High temporal resolution [sec]
dt_high_res= 1e-18; %1e-15% 
%
% Determine the DC current with respect to the STT current
STT_DC_ratio = 0;
% 
% Insert an array from outside of the simulation for the AC-STT excitation
insert_S_ac_manually = 0;
%
% Insert an array from outside of the simulation for the H_pump field
insert_H_pump_manually = 0;
%
% Generate optical Gaussian pulses in the simulation
generate_Gaussian_H_pump = 0*1;
%
% Temporal spacing between adjacent pulses, in terms of pulse duration T_p
T_btw_Gs_pls_ratio = 10; 
%
% Generate a random phase between adjacent pulses
random_phase_between_pulses = 1;
%
% Optical frequency of the Gaussian pulses - 3.7474*(10^14) corresponds to 800nm
freq_aos_gs_pls = 1*3.7474*(10^14); % 3.7474*(10^14)
%
% polarization of the Gaussian pulses: +90deg for RCP and -90deg for LCP, 0 and 180 deg for LP
pol_aos_gs_pls_phase  = 90;
%
% Peak time of the first Gaussian pulse [sec]
miu_gs = 10*10^-14;%10^-12;%%7.5*10^-12;%10^-11; %%%%%%10*10^-12;%0.25*10^-11;%% %5*10^-12;%50*10^-15;%(12.5)*10^-12; %
%
% Calculate the parameters of the analytical model of the AOS process: 
% Gamma slope, t_half and H_half
caclulate_Gamma_slope_t_half_H_half = 0;
%
% polarization of the CW signal: +90deg for RCP and -90deg for LCP, 0 and 180 deg for LP
pol_CW_phase  = 90;
%
H_CW_circ_polarization = 1; % set 0 for linear polarization and 1 for circular polarization of the CW signal
%
H_CW_circ_polarization_peak_normalization = 1; % Noramalize the peak amplitude of the CW field for all polarizations        
%
B_CW_amp = 10*10^-3; % Amplitude of the CW signal is scaled with respect to the critical eta=1 field.
B_amp_CW = B_CW_amp*Freq_CW*(1+alpha^2)/(alpha*gma);%0*1e-3; %50*1e-4; %[T] B magnetic CW induction field amplitude.
%
%% Type of sweep
% Only one of the following may be 1
J_dc_sweep = 0*1;
J_ac_sweep = 0;
RF_phase_sweep = 0;
H_pump_sweep = 0;
H_pump_phase_sweep = 0;
H_Gaussian_sweep = 0*1;
CW_polarization_sweep = 0*1;
CW_amplitude_sweep = 1;
%
if J_dc_sweep
% For STT_DC_ratio sweep
sweeping_values = linspace(0,0,1); %(0.7,1.3,800); 
end
%
if J_ac_sweep
% For J_ac sweep
  sweeping_values = logspace(4,10,400); 
end
%
if RF_phase_sweep
% For RF phase sweep
  sweeping_values = linspace(0,0,1);
end
%
if H_pump_sweep
% For H_pump field amplitude strength sweep
  sweeping_values = logspace(6,8,400);
end
%
if CW_polarization_sweep
% For pol_CW_phase sweep
sweeping_values = linspace(90,90,1);%linspace(0,360,721); 
end
%
%
if CW_amplitude_sweep
% Sweep over the amplitude of the CW , scaled with respect to the critical amplitude of w_opt*(1+alpha^2)/(alpha*miu*gma).
   H_CW_min = 10^-4;%((1+0.005^2)/(0.005))*(((1+alpha^2)/(alpha))^-1)*0.1;%0.025; %(15/91.1345)*88*(1/11)*(1/8.4)*((1+0.035^2)/(0.035))*(((1+alpha^2)/(alpha))^-1)*0.00157; %0.01; 
   H_CW_max = 1;%((1+0.005^2)/(0.005))*(((1+alpha^2)/(alpha))^-1)*0.1;%0.025; %(15/91.1345)*88*(1/11)*(1/8.4)*((1+0.035^2)/(0.035))*(((1+alpha^2)/(alpha))^-1)*0.00157; %1;
   num_amplitude_sweeps = 504; %21; %9; "5+4X"
   sweeping_values = logspace(log10(H_CW_min),log10(H_CW_max),num_amplitude_sweeps); %linspace((0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),(0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),1);% 6.7863e+10;
end
%
if H_Gaussian_sweep
   %
   % Sweep over the amplitude of the Gaussian pulse, scaled with respect to the critical amplitude of w_opt*(1+alpha^2)/(alpha*miu*gma).
   H_Gaussian_min = 0.025;%10^-4;%((1+0.005^2)/(0.005))*(((1+alpha^2)/(alpha))^-1)*0.1;%0.025; %(15/91.1345)*88*(1/11)*(1/8.4)*((1+0.035^2)/(0.035))*(((1+alpha^2)/(alpha))^-1)*0.00157; %0.01; 
   H_Gaussian_max = 0.025;%1;%((1+0.005^2)/(0.005))*(((1+alpha^2)/(alpha))^-1)*0.1;%0.025; %(15/91.1345)*88*(1/11)*(1/8.4)*((1+0.035^2)/(0.035))*(((1+alpha^2)/(alpha))^-1)*0.00157; %1;
   num_amplitude_sweeps = 1; %21;
   sweeping_values = logspace(log10(H_Gaussian_min*freq_aos_gs_pls*(1+alpha^2)/(alpha*miu*gma)),log10(H_Gaussian_max*freq_aos_gs_pls*(1+alpha^2)/(alpha*miu*gma)),num_amplitude_sweeps); %linspace((0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),(0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),1);% 6.7863e+10;
   %
   % Sweep over the the pulse duration T_p, scaled with respect to the optical cycle 2*pi/w_opt.
   T_p_min = 3.7474;%1;   %3.7474;%(3/2.7)*1000;%250;%%%1.5*37.037; %230;%9.26;%1; %
   T_p_max = 3.7474;%1000;  %(3/2.7)*1000;%250;%%%1.5*37.037; %230;%3.7474;%9.26;%1000; %1.25*3.7474
   num_T_p_sweeps = 1; %19;
   sweeping_values2 = logspace(log10(T_p_min*freq_aos_gs_pls^-1),log10(T_p_max*freq_aos_gs_pls^-1),num_T_p_sweeps);
   %
   % Sweep over the number of the Gaussian pulses:
   num_Gaussian_pulses_min = 200; 
   num_Gaussian_pulses_max = 200; 
   num_pulses_sweeps = 9; %1; %21;
   sweeping_values3 = 200;%[1 2 4 8 16 32 64 128 256];%[1 2 5 10 20 50 100 150 200];%[200]; %linspace(num_Gaussian_pulses_min,num_Gaussian_pulses_max,num_pulses_sweeps); % logspace(log10(num_Gaussian_pulses_min),log10(num_Gaussian_pulses_max),num_pulses_sweeps) %linspace((0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),(0.1*(3.7474*(10^14))*(1+alpha^2)/(alpha*miu*gma)),1);% 6.7863e+10;
   %
   % Create empty arrays
   % For a sweep over the number of pulses and field amplitude
   H_peak_num_pulses_matrix = ones(length(sweeping_values),length(sweeping_values3)); % lines correpond to varying the number of pulses, columns correspond to varying H_peak
   %
   % For a sweep over the pulse duration and field amplitude
   H_peak_Tp_matrix = ones(length(sweeping_values),length(sweeping_values2)); % lines correpond to varying T_p, columns correspond to varying H_peak
   ZZ_H_crit_matrix = ones(length(sweeping_values),length(sweeping_values2)); 
   ZZ_H_half_matrix = ones(length(sweeping_values),length(sweeping_values2));
   ZZ_Gamma_slope_matrix = ones(length(sweeping_values),length(sweeping_values2));
   %
end
%
% Create empty arrays to fill during the simulation, depending on the type
% of sweep.
%
% %%%% Use initials ZZ to be able to stack all those variables at the 
% bottom of the workspace %%%%%
%
if J_ac_sweep||J_dc_sweep||H_pump_sweep||H_pump_phase_sweep||H_Gaussian_sweep||CW_polarization_sweep||CW_amplitude_sweep
    %
    % Calculate phase between J_ac and J_sp before pulses, driven by J_ac only
    ZZ_phase_shift_in_deg_aft_pls_sz = ones(1,length(sweeping_values));
    % Calculate phase between J_ac and J_sp during pulses
    ZZ_phase_shift_in_deg_bfr_pls_sz = ones(1,length(sweeping_values));
    % Calculate J_SP amplitude before pulses
    ZZ_J_SP_amp_bfr = ones(1,length(sweeping_values));
    % Calculate J_SP amplitude during pulses
    ZZ_J_SP_amp_aft = ones(1,length(sweeping_values));
    % Calculate the Precession magnetization with J_ac infulence
    ZZ_M_z_size = ones(1,length(sweeping_values));
    % Calculate the opening angle
    ZZ_theta_array_deg = ones(1,length(sweeping_values));
    % Calculate dB gain before pulses
    ZZ_J_SP_and_Js_bfr = ones(1,length(sweeping_values));
    % Calculate dB gain during pulses
    ZZ_J_SP_and_Js_aft = ones(1,length(sweeping_values));
    %
    if J_ac_sweep
    % Save the sweeping values of J_ac
    ZZ_J_ac_amp_array = ones(1,length(sweeping_values));
    end
    %
    if J_dc_sweep
    % Save the sweeping values of J_dc
    ZZ_J_dc_amp_array = ones(1,length(sweeping_values));
    ZZ_DC_to_STT_ratio_array = ones(1,length(sweeping_values));
    end
    %
    if H_pump_sweep
    % Save the sweeping values of H_pump
    ZZ_H_pump_array = ones(1,length(sweeping_values));
    end
    %
    if H_pump_phase_sweep
    % Save the sweeping values of the phase of H_pump
    ZZ_H_pump_phase_array = ones(1,length(sweeping_values));
    end
    %
    if CW_polarization_sweep
    % Save the sweeping values of the polarization of the CW
    ZZ_CW_polarization_phase_array = ones(1,length(sweeping_values));
    end
    %
    if CW_amplitude_sweep
    % Save the sweeping values of the amplitude of the CW
    ZZ_CW_amplitude_array = ones(1,length(sweeping_values));  
    end
    %
    if H_Gaussian_sweep
    % Save the sweeping values of the phase of H_pump
    ZZ_num_pulse_array = ones(1,length(sweeping_values3));
    ZZ_H_peak_array = ones(1,length(sweeping_values));
    ZZ_T_p_array = ones(1,length(sweeping_values2));
    ZZ_M_z_size_T_p = ones(1,length(sweeping_values2));
    ZZ_M_z_size_num_pulses = ones(1,length(sweeping_values3));
    ZZ_H_crit_array = ones(1,length(sweeping_values2)); 
    ZZ_H_half_array = ones(1,length(sweeping_values2));
    Gamma_slope_array = ones(1,length(sweeping_values2));
    end
    %
end
%
if RF_phase_sweep
    %
    % Calculate phase between J_ac and J_sp before pulses, driven by J_ac only
    ZZ_phase_shift_in_deg_bfr_pls_sz = ones(1,length(sweeping_values));
    % Calculate phase between J_ac and J_sp during pulses
    ZZ_phase_shift_in_deg_aft_pls_sz = ones(1,length(sweeping_values));
    % Calculate J_SP amplitude before pulses
    ZZ_J_SP_amp_bfr = ones(1,length(sweeping_values));
    % Calculate J_SP amplitude during pulses
    ZZ_J_SP_amp_aft = ones(1,length(sweeping_values));
    % Calculate the Precession magnetization with J_ac infulence
    ZZ_M_z_size = ones(1,length(sweeping_values));
    % Calculate the opening angle
    ZZ_theta_array_deg = ones(1,length(sweeping_values));
    % Calculate dB gain before pulses
    ZZ_J_SP_and_Js_bfr = ones(1,length(sweeping_values));
    % Calculate dB gain during pulses
    ZZ_J_SP_and_Js_aft = ones(1,length(sweeping_values));
    % Save the sweeping values of RF phase
    ZZ_RF_phase_array = ones(1,length(sweeping_values));
    %
    ZZ_Mz_time_array = ones(1,length(sweeping_values));
    %
    ZZ_SP_time_array = ones(1,length(sweeping_values));
    %
    ZZ_log10_JSPJS_JS_amp_time_array = ones(1,length(sweeping_values));
    %
    ZZ_M_z_size_T_p = ones(1,length(sweeping_values2));
end
%
h_param = waitbar(0,'Sweeping parameter, please wait...','position',[335  200  270   56.25]);
% 
% This 'for' loop sweeps over EITHER: J_AC, J_DC, RF_phase, H_pump
% amplitude and polarization, Gaussian pulse amplitude.
%
for param_ind=1:length(sweeping_values)
    %
    if J_ac_sweep
       Jc_ac = sweeping_values(param_ind); %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = 0*6.4035e+07; 
       H_pump_phase = 90;
    end 
    %
    if RF_phase_sweep
       Jc_ac = 1*(10^1); %[A/m^2]
       S_ac_Phase_rf = sweeping_values(param_ind); %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       % sign of H_pump_0:
       % positive for hole (+Mz) pump, negative for electron (-Mz) pump
       H_pump_0 = 0*4.5269256327*1e+13; %[(A\m)^2]
       H_pump_phase = 90; % polarization deg of the x and y components
    end 
    %
    if J_dc_sweep
       Jc_ac = 10^1;%10^10; %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       STT_DC_ratio = sweeping_values(param_ind);
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = 0*4.5*(1+(alpha)^2)*1e+13;
       H_pump_phase = 0*90;
    end
    %
    if H_pump_sweep
       Jc_ac = 10^8;%10^8; %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = sweeping_values(param_ind);
       H_pump_phase = 90;
    end
    %
    if H_pump_phase_sweep
       Jc_ac = 5*(10^10); %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = (1e+08);
       H_pump_phase = sweeping_values(param_ind);
    end
    %
    if CW_polarization_sweep
       Jc_ac = 1*(10^1); %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = 0*(1e+08);
       H_pump_phase = 90;
       pol_CW_phase = sweeping_values(param_ind); 
    end
    %
    if CW_amplitude_sweep
       Jc_ac = 1*(10^1); %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = 0*(1e+08);
       H_pump_phase = 90;
       B_CW_amp = sweeping_values(param_ind); % Amplitude of the CW signal is scaled with respect to the critical eta=1 field.
       B_amp_CW = B_CW_amp*Freq_CW*(1+alpha^2)/(alpha*gma);%0*1e-3; %50*1e-4; %[T] B magnetic CW induction field amplitude.
    end
    %
    % The sweep of the Gaussian pulse includes the amplitude, duration and
    % number of pulses
    %
    if ~H_Gaussian_sweep
        sweeping_values2 = 1;
        sweeping_values3 = 1;
    end
 %
 % This 'for loop' sweeps over the Gaussian pulses' duration, T_p
 for param2_ind=1:length(sweeping_values2)
 %
 % This 'for loop' sweeps over the number of Gaussian pulses
 for param3_ind=1:length(sweeping_values3)
    if H_Gaussian_sweep
       Num_Gaussian_pulses = sweeping_values3(param3_ind); % Number of Gaussian pulses in the pulse train
       T_p = sweeping_values2(param2_ind);
       H_peak = sweeping_values(param_ind);
       Jc_ac = 10^1;%10^10; %[A/m^2]
       S_ac_Phase_rf = 0; %[deg]
       Jc_dc = STT_DC_ratio*((alpha*w/gma)*(2*q*M0*Thick))/((1+alpha^2)*hbar*Theta_sh); %[A/m^2]
       H_pump_0 = 0*4.5*(1+(1*0.025)^2)*1e+13;
       H_pump_phase = 90;
    end
    %
    H_pump_phase_rad = (2*pi/360)*H_pump_phase;
    B_Pulse_amp = 0*0.15;%0.01;
    param_value = sweeping_values(param_ind);
    waitbar((param_ind)/length(sweeping_values),h_param)

    %% Magnetic field settings
    
    B0 = 0*1*(w/gma+B0_anis/2); %+0.5*miu*H0_anis;
    half_delta_B = w*alpha/gma;
    B_samples = 1; % Choose an odd B_sample number so that one B will be resonant
    B_width = 4;   %in fwhm units
    Resonanse_Field = w/gma;
    res_field_over_pls_amp = 0;
    if abs(B_Pulse_amp) > 0
        res_field_over_pls_amp = Resonanse_Field/B_Pulse_amp;
    end
    % For a single B external value, in resonance (also, manually set B_samples = 1):
    B_Xtrnl_vect = B0;
    % For a single B external value' off resonance
    %B_Xtrnl_vect = B0-B_width*half_delta_B;
    if B_samples > 1
    % For multiple B external field values:
        B_Xtrnl_vect = linspace(B0-B_width*half_delta_B,B0+B_width*half_delta_B,B_samples);  %[T] External Magnetic flux density;
    end
    %
    H_Xtrnl_vect= B_Xtrnl_vect/miu;
    
    %% External field Direction
    Field_tilt_angle= 0*(0.1/360)*2*pi;% (0.1/360)*2*pi; 
    H_ext_dirctn=[1 0 tan(Field_tilt_angle)]';   %external field direction only
    H_ext_unit=H_ext_dirctn/norm(H_ext_dirctn);   %field unit vector
    
    %% Demagnetizing Field
    % Nx+Ny+Nz = 1 in SI units
    % Corresponds to shape anisotropy
    Nx = 0;
    Ny = 1;
    Nz = 0;
    N0 = 0*1;
    N_Dmag = N0*[Nx Ny Nz]';
    
    %% STT DC Components
    
    Sx=0;     
    Sy=0;
    Sz=1;
    S_vector=[Sx Sy Sz]';  %spin polarization vector

    %% Generate time sections
    % The three time sections are designed for a sharp pulse perturbation
    % in the second section, in which dt is smaller than the other two.
    T_section{1}=0:dt_low_res:100*dt_low_res;%2e-16;%1e-13%1e-9;
%    T_section{1}=0:dt_low_res:3.5e-07;
    dT_section{1}=dt_low_res;
%    T_section{2}=T_section{1}(end):dt_low_res:T_section{1}(end)+100*300e-12; % time sections are overlaping in start_end points 
%     T_section{2}=T_section{1}(end):dt_high_res:T_section{1}(end)+ 510e-9; % time sections are overlaping in start_end points 
%     dT_section{2}=dt_high_res;
    T_section{2}=T_section{1}(end):dt_low_res:T_section{1}(end)+ 100*dt_low_res;%2e-16;%5e-13%5e-9; % time sections are overlaping in start_end points 
    dT_section{2}=dt_low_res;
    T_section{3}=T_section{2}(end):dt_low_res:2*10^-11;%2*10^-12;%26*10^-12;%15*10^-12;%20*10^-12;%0.25*2e-11;%%100*10^-15;%1e-12;%2.5e-11;%1e-13;%
    dT_section{3}=dt_low_res;
    %
    T_Global=[T_section{1}(1:end-1) T_section{2}(1:end-1) T_section{3}];  % Index starts at 2 since time sections are overlaping in start_end points 
    %% Find The Steady State Teta Angle With Z Axis
    psi=acos(H_ext_unit(3)) ;       % angle of external field with the Z axis
    for ind=1:length(H_Xtrnl_vect)
        %options = optimset('TolFun',1e-15);
        %teta0(ind) = fminbnd(@(x) abs((sin(2*x)/sin(psi-x) - miu*M0*H_Xtrnl_vect(ind)/Ku)) ,-pi/2,pi/2,options);   %teta0 is the steady state angle
        teta_v=0:0.00001:pi/2;
        teta_func=abs(Ku*sin(2*teta_v)-miu*M0^2*sin(teta_v).*cos(teta_v)-miu*M0*H_Xtrnl_vect(ind)*sin(psi-teta_v));
        %     plot(teta_v/pi*180,teta_func,[teta_v(1)/pi*180 teta_v(end)/pi*180] ,[0 0])
        %     return
        [a_tmp,ind_tmp]=min(teta_func);
        teta0(ind)=teta_v(ind_tmp);
        phi_abs(ind)=atan(abs(H_ext_unit(2)/H_ext_unit(1)));      %this is the angle in the xy plane relative to x axis (absolute value ).
        
        Mx=sign(H_ext_unit(1))*sin(teta0(ind))*cos(phi_abs(ind));
        My=sign(H_ext_unit(2))*sin(teta0(ind))*sin(phi_abs(ind));
        Mz=cos(teta0(ind));
        
%     M_SS_dirctn_unit(ind,1:3)=[Mx My Mz]';                  we turned off the S.S  direction and replace it with External  Field direction (Anisotropies turned off)  
    
        M_SS_dirctn_unit(ind,1:3)=H_ext_unit';  
        
        dM_SS_dirctn_unit(ind,1:3)=cross(M_SS_dirctn_unit(ind,1:3),[0 0 1]);  % I chose the direction of the initial perturbation to be orthogonal to Z and M_ss
        dM_SS_dirctn_unit(ind,1:3)=dM_SS_dirctn_unit(ind,1:3)/norm(dM_SS_dirctn_unit(ind,1:3));
        
        M0_0_dirctn(ind,1:3)=M_SS_dirctn_unit(ind,1:3)+dM0*dM_SS_dirctn_unit(ind,1:3);
        M0_0_dirctn(ind,1:3)=M0_0_dirctn(ind,1:3)/norm(M0_0_dirctn(ind,1:3));
        
    end
    %
    %% Spin Pumping
    % Spin Pumping effect occurs at FM NM heterostructure.
    % Spin Pumping enlarges the Gilbert losses, by an additional 
    % lose channel through spin transfer to a nearby NM. 

    % Spin Pumping conductivity
    g_SP = (10^19);
    %
    % Formula from Saitoh:

    % Efficiency of Spin Pumping PHENOMENOLOGICAL only (value between 0 to 1)
    % To nullify the Spin Pumping, set eff_SP = 0.
    % In order to nullify the SP effect in general, also set g_SP = 0
    %eff_SP = 0*0.5*1;

    % alpha_SP = hbar*gma*g_SP*eff_SP*((4*pi*Thick*M0)^-1);

    % Take the sign with respect to the direction between Js-Jsp
    sign_J_s_J_sp = -1;
    %
    % Describes the ratio between the emmited to the backscattered SP (closed loop) -
    % which enlarges the Gilbert damping, by the SP that exits the FM.
    SP_efficiency_ratio = 0.1;
    %
    % Formula from the simulation
    alpha_SP = (-sign_J_s_J_sp)*hbar*gma*g_SP*SP_efficiency_ratio*((4*pi)^-1);
    
    % Enlarge the DC STT current with the new damping:
    Jc_dc = ((alpha + alpha_SP)/alpha)*((1+alpha^2)/(1+(alpha + alpha_SP)^2))*Jc_dc;
  
    %% STT DC parameter
    Hshe_dc=(hbar*Theta_sh*Jc_dc)/(2*q*M0*Thick);  %spin Hall effect field parameter 
    %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% H Field iteration loop
    h = waitbar(0,'Please wait...','position',[650  200  270   56.25]); 
 
    for L=1:length(H_Xtrnl_vect)
        phase_shift_in_radians = [];
        phase_shift_in_radians_bfr_pls = [];
        phase_shift_in_radians_aft_pls = [];
        S_ac_array_initial{L} = [];
        S_ac_array_global{L} = [];
        M_array_global{L} = [];
        S_sp_array_global{L} = [];
        H_CW_array_global{L} = [];
        H_pump_array_global{L} = [];
        H_IFE_global{L} = [];
        gamma_v_global{L} = [];
        gamma_c_global{L} = [];
        ReV12_global{L} = [];
        ImV12_global{L} = [];
        Lambda_h_global{L} = [];
        Lambda_e_global{L} = [];
        w_TLS_global{L} = [];
        gamma_inh_global{L} = [];
        primary_torque_array_global{L} = [];
        damping_torque_array_global{L} = [];
        
        waitbar(L/length(H_Xtrnl_vect),h)
        for Section_ind=1:1:3                        % Iterate on each time section for loop
            clear H_CW H_CW_an H_CW_bn H_CW_cn ...
                  S_ac S_ac_an S_ac_bn S_ac_cn...
                M_tg gma_LL_0 gma_LL_an gma_LL_bn gma_LL_cn H_Anis H_Tot K0_0 K0_an K0_bn K0_cn LL_lada_0 LL_lada_an LL_lada_bn LL_lada_cn...
                t t_xtrct M M_mat M_mat_nrm M_mat_xtrct M_mat_xtrct_nrm S ...
                sig_rf sig_rf_an sig_rf_bn sig_rf_cn...
                S_ac_sig_rf S_ac_sig_rf_an S_ac_sig_rf_bn S_ac_sig_rf_cn ...
                A_0 A_an A_bn A_cn dA_0 dA_an dA_bn dA_cn...
                gma_LLG_0 gma_LLG_an gma_LLG_bn gma_LLG_cn alpha_0 alpha_an alpha_bn alpha_cn H_Pulsed_0 H_Pulsed_an H_Pulsed_bn H_Pulsed_cn
            
            t=T_section{Section_ind};
            dt=dT_section{Section_ind};
            
            Len_t=length(t);            
            M=zeros(3,Len_t);
                        
            %% CW Field settings
            H_CW=ones(3,Len_t);
            H_CW_dirctn=[1 1 0]';     %RF field direction only
            %H_CW_unit=H_CW_dirctn/norm(H_CW_dirctn);   %Rf field unit vector
            H_CW_unit=H_CW_dirctn;   % For the unnormalized field, in CP xy polarization      
            %
            pol_CW_phase_rad  = pol_CW_phase*(2*pi/360); % Convert into radian units
            H_amp_CW=B_amp_CW/miu;   %[A/m]   H magnetic RF field amplitude.
            if H_CW_circ_polarization_peak_normalization
               H_amp_CW=H_amp_CW/sqrt(1+abs(cos(pol_CW_phase_rad))); % Noramalize the peak amplitude of the CW field for all polarizations
            end
            Phase_rf=0;              %[deg]  -NOTE: This is in degrees.
            Phase_rf_rad=(Phase_rf/360)*(2*pi); %[rad]
            %
            sig_rf=H_amp_CW*cos(2*pi*Freq_CW*t+Phase_rf_rad);
            sig_rf_an=H_amp_CW*cos(2*pi*Freq_CW*(t+dt/2)+Phase_rf_rad);   %this adds with the bn eq.
            sig_rf_bn=H_amp_CW*cos(2*pi*Freq_CW*(t+dt/2)+Phase_rf_rad);   %this field adds with the cn eq.
            sig_rf_cn=H_amp_CW*cos(2*pi*Freq_CW*(t+dt)+Phase_rf_rad);     %this field adds with the dn eq.
            
            % Linearly polarized H_CW:
            if ~H_CW_circ_polarization
                for i = 1:1:3
                H_CW(i,:)=H_CW_unit(i)*sig_rf;
                H_CW_an(i,:)=H_CW_unit(i)*sig_rf_an;
                H_CW_bn(i,:)=H_CW_unit(i)*sig_rf_bn;
                H_CW_cn(i,:)=H_CW_unit(i)*sig_rf_cn;
                end
            end
            
            % Circularly polarized H_CW:
            % polarized in the x_y plane 
            if H_CW_circ_polarization                        
                sig_rf_pol=H_amp_CW*cos(2*pi*Freq_CW*t+Phase_rf_rad-pol_CW_phase_rad);
                sig_rf_an_pol=H_amp_CW*cos(2*pi*Freq_CW*(t+dt/2)+Phase_rf_rad-pol_CW_phase_rad);   %this adds with the bn eq.
                sig_rf_bn_pol=H_amp_CW*cos(2*pi*Freq_CW*(t+dt/2)+Phase_rf_rad-pol_CW_phase_rad);   %this field adds with the cn eq.
                sig_rf_cn_pol=H_amp_CW*cos(2*pi*Freq_CW*(t+dt)+Phase_rf_rad-pol_CW_phase_rad);     %this field adds with the dn eq.
                %
                H_CW(1,:)=H_CW_unit(1)*sig_rf;
                H_CW(2,:)=H_CW_unit(2)*sig_rf_pol;
                H_CW(3,:)=H_CW_unit(3)*0;
                H_CW_an(1,:)=H_CW_unit(1)*sig_rf_an;
                H_CW_an(2,:)=H_CW_unit(2)*sig_rf_an_pol;
                H_CW_an(3,:)=H_CW_unit(3)*0;
                H_CW_bn(1,:)=H_CW_unit(1)*sig_rf_bn;
                H_CW_bn(2,:)=H_CW_unit(2)*sig_rf_bn_pol;
                H_CW_bn(3,:)=H_CW_unit(3)*0;
                H_CW_cn(1,:)=H_CW_unit(1)*sig_rf_cn;
                H_CW_cn(2,:)=H_CW_unit(2)*sig_rf_cn_pol;
                H_CW_cn(3,:)=H_CW_unit(3)*0;  
            end
            %
            %% Time Dependent Spin Polarization
            %
            Hshe_ac=(hbar*Theta_sh*Jc_ac)/(2*q*M0*Thick);  %spin Hall effect field parameter 

            S_ac=ones(3,Len_t);
            
            S_ac_dirctn=[1 0 0]';                        %RF field direction only
            S_ac_unit=S_ac_dirctn/norm(S_ac_dirctn);   %Rf field unit vector
            
            S_ac_Phase_rf_rad=(S_ac_Phase_rf/360)*(2*pi); %[rad]
            %
%             T_start_S_ac= (0*0.5+2500*450.5)*(10^-9); % Time in which S_ac signal starts
%             T_end_S_ac = 2500*1.3*10^-6; % Time in which S_ac signal stops
            T_start_S_ac= 0*0.5*10^-6;%4*10^-7; %150*10^-9; % Time in which S_ac signal starts
            T_end_S_ac = T_section{3}(end);%2.85*(10^-7)-T_section{2}(end);%% Time in which S_ac signal stops
            %
            S_ac_sig_rf=(cos(2*pi*Freq_CW*t+S_ac_Phase_rf_rad).*((T_start_S_ac < t)&(t < T_end_S_ac)));
            S_ac_sig_rf_an=(cos(2*pi*Freq_CW*(t+dt/2)+S_ac_Phase_rf_rad).*((T_start_S_ac < t)&(t < T_end_S_ac)));   %this adds with the bn eq.
            S_ac_sig_rf_bn=(cos(2*pi*Freq_CW*(t+dt/2)+S_ac_Phase_rf_rad).*((T_start_S_ac < t)&(t < T_end_S_ac)));   %this field adds with the cn eq.
            S_ac_sig_rf_cn=(cos(2*pi*Freq_CW*(t+dt)+S_ac_Phase_rf_rad).*((T_start_S_ac < t)&(t < T_end_S_ac)));     %this field adds with the dn eq.
            %
            for i = 1:1:3
                S_ac(i,:)=S_ac_unit(i)*S_ac_sig_rf;
                S_ac_an(i,:)=S_ac_unit(i)*S_ac_sig_rf_an;
                S_ac_bn(i,:)=S_ac_unit(i)*S_ac_sig_rf_bn;
                S_ac_cn(i,:)=S_ac_unit(i)*S_ac_sig_rf_cn;
            end
            % 
            % Recover the last S_ac from the previous section
            if Section_ind > 1
                S_ac(:,1) = S_ac_between_sections;
                S_ac_an(:,1) = S_ac_between_sections_an;
                S_ac_bn(:,1) = S_ac_between_sections_bn;
                S_ac_cn(:,1) = S_ac_between_sections_cn;
            end
            %
            %% Anisotropy Recovery Profile
            %
            % Pulse at steady state
            t0 = 1e-6;%130.5e-09;%150.5e-09;
            %t0=30*2e-09;
            % Pulse at non-steady state
            %t0=10*1e-09;
            %
            tau=0.4e-12;
            t1=t0+1e-12;
            tau1=4e-9;
            prtrb_depth=0;%0.5;
            
            K0_0=K0*(1-prtrb_depth*Gen_Pls_Fermi_Fermi(t,t0,tau,t1,tau1));
            K0_an=K0*(1-prtrb_depth*Gen_Pls_Fermi_Fermi(t+dt/2,t0,tau,t1,tau1));
            K0_bn=K0*(1-prtrb_depth*Gen_Pls_Fermi_Fermi(t+dt/2,t0,tau,t1,tau1));
            K0_cn=K0*(1-prtrb_depth*Gen_Pls_Fermi_Fermi(t+dt,t0,tau,t1,tau1));
            
%             sig_rcrd=[sig_rcrd(1:end-1) K0_0];
%             if Section_ind==3
%               plot(T_Global,sig_rcrd,T_section{1}(1:end-1),sig_rcrd(1:length(T_section{1})-1),T_section{3},sig_rcrd(end-length(T_section{3})+1:end))
%             end
                       
            %% Magnetization Recovery Profile
            %
            t0_A=t0; 
            tau_A=5e-12;%150e-12;
            t1_A=t0+30e-12;
            tau1_A=0.05e-9;%tau;
            prtrb_depth_A=0*0.02;%0.026;
            t2_depth_A=0*0.05;
            tau_2_A=0.5e-9; 
            %
            A_0=1-prtrb_depth_A*Gen_Pls_Fermi_Fermi_3T(t,t0_A,tau_A,t1_A,tau1_A,t2_depth_A,tau_2_A,T_Global);
            A_an=1-prtrb_depth_A*Gen_Pls_Fermi_Fermi_3T(t+dt/2,t0_A,tau_A,t1_A,tau1_A,t2_depth_A,tau_2_A,T_Global);
            A_bn=1-prtrb_depth_A*Gen_Pls_Fermi_Fermi_3T(t+dt/2,t0_A,tau_A,t1_A,tau1_A,t2_depth_A,tau_2_A,T_Global);
            A_cn=1-prtrb_depth_A*Gen_Pls_Fermi_Fermi_3T(t+dt,t0_A,tau_A,t1_A,tau1_A,t2_depth_A,tau_2_A,T_Global);
            
            %             sig_rcrd=[sig_rcrd(1:end-1) A_0];
            %             if Section_ind==3
            %                 hold on
            %                 plot(T_Global,sig_rcrd,T_section{1}(1:end-1),sig_rcrd(1:length(T_section{1})-1),T_section{3},sig_rcrd(end-length(T_section{3})+1:end))
            %                 %figure;plot(T_section{2}(1:end),sig_rcrd(length(T_section{1}):end-length(T_section{3})+1))
            %                 return
            %             end
            
            
            warning('what to do with +dt/2')
            warning('differentiation not continuous')            
            warning('one differentiation might be missing')
            
            dA_0=diff(A_0)/dt;   %dA_0=[dA_0 dA_0(end)];          % last point is artifitially duplicated to get equal length vectors
            dA_an=diff(A_an)/dt; %dA_an=[dA_an dA_an(end)];       % last point is artifitially duplicated to get equal length vectors
            dA_bn=diff(A_bn)/dt; %dA_bn=[dA_bn dA_bn(end)];       % last point is artifitially duplicated to get equal length vectors
            dA_cn=diff(A_cn)/dt; %dA_cn=[dA_cn dA_cn(end)];       % last point is artifitially duplicated to get equal length vectors
            
            %             sig_rcrd=[sig_rcrd(1:end-1) dA_0];
            %             if Section_ind==3
            %                 plot([T_section{1}(1:end-2) T_section{2}(1:end-2) T_section{3}(1:end-1)],sig_rcrd,...
            %                     T_section{1}(1:end-2),sig_rcrd(1:length(T_section{1})-2),T_section{3}(1:end-1),sig_rcrd(end-length(T_section{3})+2:end))
            %                 return
            %             end
            
            
            %% Damping recovery profile
            % alpha is the Gilbert dumping coefficient
            t0_alpha=t0; 
            tau_alpha=3e-12;
            t1_alpha=t0+30e-12;
            tau1_alpha=0.5e-9;%tau;
            prtrb_depth_alpha=0;%0.5;
                      
            alpha_0=alpha*(1-prtrb_depth_alpha*Gen_Pls_Fermi_Fermi(t,t0_alpha,tau_alpha,t1_alpha,tau1_alpha));
            alpha_an=alpha*(1-prtrb_depth_alpha*Gen_Pls_Fermi_Fermi(t+dt/2,t0_alpha,tau_alpha,t1_alpha,tau1_alpha));
            alpha_bn=alpha*(1-prtrb_depth_alpha*Gen_Pls_Fermi_Fermi(t+dt/2,t0_alpha,tau_alpha,t1_alpha,tau1_alpha));
            alpha_cn=alpha*(1-prtrb_depth_alpha*Gen_Pls_Fermi_Fermi(t+dt,t0_alpha,tau_alpha,t1_alpha,tau1_alpha));
            
            %sig_rcrd=[sig_rcrd(1:end-1) alpha_0];
            %if Section_ind==3
            %   plot(T_Global,sig_rcrd,T_section{1}(1:end-1),sig_rcrd(1:length(T_section{1})-1),T_section{3},sig_rcrd(end-length(T_section{3})+1:end))
            %return
            %end
            
            %% Gamma Recovery Profile
            % gma is the gyromagnetic coefficient
            %gma_LLG=gma;
            %LL_lada=gma_LLG*alpha/(1+alpha^2);
            %gma_LL=gma_LLG/(1+alpha^2);
            t0_gma_LLG=t0; 
            tau_gma_LLG=5e-12;
            t1_gma_LLG=t0+30e-12;
            tau1_gma_LLG=0.4e-9;%tau;
            prtrb_depth_gma_LLG=0;%0.05;
            
            gma_LLG_0=gma*(1-prtrb_depth_gma_LLG*Gen_Pls_Fermi_Fermi(t,t0_gma_LLG,tau_gma_LLG,t1_gma_LLG,tau1_gma_LLG));
            gma_LLG_an=gma*(1-prtrb_depth_gma_LLG*Gen_Pls_Fermi_Fermi(t+dt/2,t0_gma_LLG,tau_gma_LLG,t1_gma_LLG,tau1_gma_LLG));
            gma_LLG_bn=gma*(1-prtrb_depth_gma_LLG*Gen_Pls_Fermi_Fermi(t+dt/2,t0_gma_LLG,tau_gma_LLG,t1_gma_LLG,tau1_gma_LLG));
            gma_LLG_cn=gma*(1-prtrb_depth_gma_LLG*Gen_Pls_Fermi_Fermi(t+dt,t0_gma_LLG,tau_gma_LLG,t1_gma_LLG,tau1_gma_LLG));
            
            %             sig_rcrd=[sig_rcrd(1:end-1) gma_LLG_0];
            %             if Section_ind==3
            %               plot(T_Global,sig_rcrd,T_section{1}(1:end-1),sig_rcrd(1:length(T_section{1})-1),T_section{3},sig_rcrd(end-length(T_section{3})+1:end))
            %             return
            %             end
            
            gma_LL_0=gma_LLG_0./((1+alpha_0.^2));
            gma_LL_an=gma_LLG_an./(1+alpha_an.^2);    %this adds with the bn eq.
            gma_LL_bn=gma_LLG_bn./(1+alpha_bn.^2);    %this field adds with the cn eq.
            gma_LL_cn=gma_LLG_cn./(1+alpha_cn.^2);    %this field adds with the dn eq.
                        
            LL_lada_0=gma_LLG_0.*alpha_0./(1+alpha_0.^2);
            LL_lada_an=gma_LLG_an.*alpha_an./(1+alpha_an.^2);    %this adds with the bn eq.
            LL_lada_bn=gma_LLG_bn.*alpha_bn./(1+alpha_bn.^2);    %this field adds with the cn eq.
            LL_lada_cn=gma_LLG_cn.*alpha_cn./(1+alpha_cn.^2);     %this field adds with the dn eq.
                        
            %sig_rcrd=[sig_rcrd(1:end-1) LL_lada_0];
            %if Section_ind==3
            %     plot(T_Global,sig_rcrd,T_section{1}(1:end-1),sig_rcrd(1:length(T_section{1})-1),T_section{3},sig_rcrd(end-length(T_section{3})+1:end))
            %return
            %end
            
            %% External pulsed Field profile stemming from a sharp change in the magnetization.
            %Pulsed_Field_tilt_angle=-89.99/360*2*pi;
            
            %H_Pulsed_dirctn=[-1 0 tan(Pulsed_Field_tilt_angle)]';   %external field direction only
            %H_Pulsed_unit=H_Pulsed_dirctn/norm(H_Pulsed_dirctn);   %field unit vector
            H_Pulsed_unit=[0 1 0];%-M0_0_dirctn(L,1:3);
            
            %B_Pulse_amp=400e-4;%1e1*20000*1e-4;       %50*1e-4; %[T] B magnetic RF induction field amplitude.
            H_Pulse_amp=B_Pulse_amp/miu;     %[A/m]   H magnetic RF field amplitude.
            t0_H_pulsed= 1.6*10^-14;%0.8*T_end_S_ac;%(5*(10^-16)-0*T_section{2}(end));%0.2e-6;%4e-6;
            tau_H_pulsed=0.4e-12;
            t1_H_pulsed=t0_H_pulsed+3e-12;
            tau1_H_Pulsed=0.4e-12;
            % For a pulse train:
            num_pulses = 1;%50;
            time_int_betw_pls = 10^-8;%27*10^-7;%76*10^-7;%
            %
            % when using pulses, set accordingly
            %t0_H_pulsed_end = t0_H_pulsed + time_int_betw_pls*num_pulses;
            %
            %when NOT using pulses, set arbitrarly
            t0_H_pulsed_end  = 20*10^-12 ;%2*10^-14;%T_end_S_ac; %5*10^-15;%
            %
            % Number of time samples for calculations of phase \ amplitude.
            % Takes into account the last 5 pulses
            %T_int = 5*(time_int_betw_pls/dt_low_res);
            % For STT simulation without pulses, to measure the J_s J_sp
            % phase only in presence of J_s
            % short simulation
            %T_int = round(0.2*(time_int_betw_pls/dt_low_res));
            % long simulation
            T_int_fraction = 0.01;
            T_int_fraction_small = 0.001;
            T_int1 = round(T_int_fraction*(t0_H_pulsed/dt));
            T_int11 = round(T_int1*T_int_fraction_small);
            T_int2 = round(T_int_fraction*((t0_H_pulsed_end-t0_H_pulsed)/dt));
            T_int22 = round(T_int2*T_int_fraction_small);
%             %
%             S_ac_indx = uint64((t0_H_pulsed-t(1))/dt);
%             %
            H_sig_Pulsed_0 = zeros(size(S_ac_sig_rf));
            H_sig_Pulsed_an = zeros(size(S_ac_sig_rf_an));
            H_sig_Pulsed_bn = zeros(size(S_ac_sig_rf_bn));
            H_sig_Pulsed_cn = zeros(size(S_ac_sig_rf_cn));

          % For a pulse train - 
          % In each time section, go over all pulses and see of they belong
          % to the time section
          %
            for n = 1:(num_pulses)
                % check if the pulse is in the time section
                if (t(1)<t0_H_pulsed + time_int_betw_pls*(n-1))&(t0_H_pulsed + time_int_betw_pls*(n-1)<t(length(t)))
                  % Phase locked to the Jac 
                    H_sig_Pulsed_0=H_sig_Pulsed_0 + H_Pulse_amp*Gen_Pls_Fermi_Fermi(t,t0_H_pulsed + time_int_betw_pls*(n-1),tau_H_pulsed,t1_H_pulsed + time_int_betw_pls*(n-1),tau1_H_Pulsed);
                    H_sig_Pulsed_an=H_sig_Pulsed_an + H_Pulse_amp*Gen_Pls_Fermi_Fermi(t+dt/2,t0_H_pulsed + time_int_betw_pls*(n-1),tau_H_pulsed,t1_H_pulsed + time_int_betw_pls*(n-1),tau1_H_Pulsed);
                    H_sig_Pulsed_bn=H_sig_Pulsed_bn + H_Pulse_amp*Gen_Pls_Fermi_Fermi(t+dt/2,t0_H_pulsed + time_int_betw_pls*(n-1),tau_H_pulsed,t1_H_pulsed + time_int_betw_pls*(n-1),tau1_H_Pulsed);
                    H_sig_Pulsed_cn=H_sig_Pulsed_cn + H_Pulse_amp*Gen_Pls_Fermi_Fermi(t+dt,t0_H_pulsed + time_int_betw_pls*(n-1),tau_H_pulsed,t1_H_pulsed + time_int_betw_pls*(n-1),tau1_H_Pulsed);
                end
            end  
            %
            % Set the pulse as a vector
            %
            H_Pulsed_0(1,:)=H_Pulsed_unit(1)*H_sig_Pulsed_0;
            H_Pulsed_0(2,:)=H_Pulsed_unit(2)*H_sig_Pulsed_0;
            H_Pulsed_0(3,:)=H_Pulsed_unit(3)*H_sig_Pulsed_0;
            H_Pulsed_an(1,:)=H_Pulsed_unit(1)*H_sig_Pulsed_an;
            H_Pulsed_an(2,:)=H_Pulsed_unit(2)*H_sig_Pulsed_an;
            H_Pulsed_an(3,:)=H_Pulsed_unit(3)*H_sig_Pulsed_an;
            H_Pulsed_bn(1,:)=H_Pulsed_unit(1)*H_sig_Pulsed_bn;
            H_Pulsed_bn(2,:)=H_Pulsed_unit(2)*H_sig_Pulsed_bn;
            H_Pulsed_bn(3,:)=H_Pulsed_unit(3)*H_sig_Pulsed_bn;
            H_Pulsed_cn(1,:)=H_Pulsed_unit(1)*H_sig_Pulsed_cn;
            H_Pulsed_cn(2,:)=H_Pulsed_unit(2)*H_sig_Pulsed_cn;
            H_Pulsed_cn(3,:)=H_Pulsed_unit(3)*H_sig_Pulsed_cn;
            %
            
            %% initialize variables
            %
            M_array{length(H_Xtrnl_vect)}=[];       %prealoocate cell array
            S{length(H_Xtrnl_vect)}=[];
            f_shift{length(H_Xtrnl_vect)}=[];
            freq(length(H_Xtrnl_vect))=0;
            H_Tot=zeros(3,length(t));
            H_Anis=zeros(3,length(t));
                        
            H_Xtrnl=H_Xtrnl_vect(L)*H_ext_unit;      % [A/m] Magnetic field;
            
            if Section_ind==1
                M_tg(:,1)=A_0(1)*M0*M0_0_dirctn(L,1:3);    % Starting point for the magnetization.
            else
                M_tg(:,1)=M_array_global{L}(:,end);
            end
            %
            %% S_ac feedback loop response  
            %
            % To choose either an  open \ closed loop of SP currrent
            % Open loop - the SP current exits the FM and doesn't interact
            % with the magnetization.
            %
            % Closed loop - the SP current interacts with the magnetization
            SP_closed_loop = 1;
            %
            % Feedback coefficients of the dM/dt SP term (usually negligble)
            Im_SP_xy = 0*4.33*10^-17; % "dM/dt J_sp" in the x&y directions
            Im_SP_z = 0*4.33*10^-17; % "dM/dt J_sp" in the z directions
            %
            % Feedback coefficients of the M x dM/dt SP term (usually dominant)
            % For isotropic M x dM/dt SP terms:
            Re_SP = g_SP*hbar/(4*pi);
            % For anisotropic M x dM/dt SP terms:
            Re_SP_xy = 0*g_SP*hbar/(4*pi); % "(M x dM/dt) J_sp" in the x&y directions
            Re_SP_z = 0*g_SP*hbar/(4*pi); % "(M x dM/dt) J_sp" in the z directions
            %
            % Phenomenological term that arises from a SP proportional to
            % dM/dt. Analogous to alpha_SP (which is introduced phenomenologically) vs A_feedback_xy/z_SP ('implemented')
            % To nullify this response, set A_feedback = 0
            A_feedback = 0*50*2.93*10^-17;
            A_Global = 1 + (A_feedback*gma*Hshe_ac*M0); % Total feedback response - divides all LLGS right handside by A_Global
            %
            % A_feedback_xy / _z behave like A_Global for the overdamped regime,
            % and differ for the underdamped regime.
            %
            % Define an array that will retain the inital wave S_ac (no feedback yet)
            S_ac_initial = [];
            S_ac_initial(:,1) = S_ac(:,1);
            % Define an initial Spin Pumping vector
            S_sp = [];
            S_sp(:,1) = [0 0 0]';
            %
            if Section_ind > 1
                S_sp(:,1) = S_sp_temp;
                t0 = t0_between_sections;
                t0_H_pulsed = t0_H_pulsed_between_sections;
                T_int1 = T_int1_between_sections;
                T_int2 = T_int2_between_sections;
            end
            %
            %% TLS terms
            %
            gamma_v = zeros(1,length(t));
            gamma_c = zeros(1,length(t));
            ReV12 = zeros(1,length(t));
            ImV12 = zeros(1,length(t));
            Lambda_h = zeros(1,length(t));
            Lambda_e = zeros(1,length(t));
            w_TLS = zeros(1,length(t));
            gamma_inh = zeros(1,length(t));
            %
            if Section_ind == 1
                % Assume that the pump field is 0 in the beginning, as its
                % envelope diverges
                gamma_v(:,1) = -gma*miu*alpha*((M0*(1+(alpha^2)))^-1)*((H_Xtrnl_vect*(1-STT_DC_ratio)+(M_tg(3,1))*N_Dmag(3))*(M0-(M_tg(3,1)))-0*H_pump_0*cos(H_pump_phase_rad));
                gamma_c(:,1) = gma*miu*alpha*((M0*(1+(alpha^2)))^-1)*((H_Xtrnl_vect*(1-STT_DC_ratio)+(M_tg(3,1))*N_Dmag(3))*(M0+(M_tg(3,1)))+0*H_pump_0*cos(H_pump_phase_rad));
                ReV12(:,1) = gma*miu*((2*(1+(alpha^2)))^-1)*((M_tg(1,1))*N_Dmag(1)-(alpha/M0)*(M_tg(3,1))*((M_tg(2,1))*N_Dmag(2)+sig_rf(1)));
                ImV12(:,1) = -gma*miu*((2*(1+(alpha^2)))^-1)*((alpha/M0)*((M_tg(3,1)))*(M_tg(1,1))*N_Dmag(1)+((M_tg(2,1))*N_Dmag(2)+sig_rf(1)))+(gma/(2*M0))*(M_tg(3,1))*Hshe_ac*S_ac_initial(1,1);
                Lambda_h(:,1) = 0*gma*miu*H_pump_0*((1+(alpha^2))^-1)*(sin(H_pump_phase_rad)-alpha*cos(H_pump_phase_rad));
                lambda_e(:,1) = 0*gma*miu*H_pump_0*((1+(alpha^2))^-1)*(sin(H_pump_phase_rad)+alpha*cos(H_pump_phase_rad));
                w_TLS(:,1) = gma*miu*((1+(alpha^2))^-1)*(H_Xtrnl_vect+(M_tg(3,1))*N_Dmag(3)+(1+(alpha)^2)*(1/(miu*M0))*(M_tg(2,1))*Hshe_ac*S_ac_initial(1,1)-(M_tg(3,1))*0*H_pump_0*cos(H_pump_phase_rad)/(M0^2 - ((M_tg(3,1)))^2)+(alpha/M0)*((((M0)^2)*0*H_pump_0*sin(H_pump_phase_rad)/(M0^2 - ((M_tg(3,1)))^2))-((N_Dmag(1)-N_Dmag(2))*((M_tg(1,1))))*((M_tg(2,1)))+(M_tg(1,1))*sig_rf(1)));
                gamma_inh(:,1) = gma*miu*((1+(alpha^2))^-1)*((M_tg(3,1))/M0)*((0*H_pump_0/(M0^2 - ((M_tg(3,1)))^2))*(alpha*cos(H_pump_phase_rad)*(M_tg(3,1))-M0*sin(H_pump_phase_rad))+alpha*(H_Xtrnl_vect*(1-STT_DC_ratio)+(M_tg(3,1))*N_Dmag(3)));
                %
            end
             % Recover the last TLS terms from the previous section
            if Section_ind > 1
                gamma_v(:,1) = gamma_v_between_sections;
                gamma_c(:,1) = gamma_c_between_sections;
                ReV12(:,1) = ReV12_between_sections;
                ImV12(:,1) = ImV12_between_sections;
                Lambda_h(:,1) = Lambda_h_between_sections;
                Lambda_e(:,1) = Lambda_e_between_sections;
                w_TLS(:,1) = w_TLS_between_sections;
                gamma_inh(:,1) = gamma_inh_between_sections;
                primary_torque(:,1) = primary_torque_between_sections;
                damping_torque(:,1) = damping_torque_between_sections;
            end
            %
            %% Pump field array
            %
            H_pump = zeros(3,length(t));
            H_pump_an = zeros(3,length(t));
            H_pump_bn = zeros(3,length(t));
            H_pump_cn = zeros(3,length(t));
            %
            H_IFE = zeros(3,length(t));
            H_IFE_an = zeros(3,length(t));
            H_IFE_bn = zeros(3,length(t));
            H_IFE_cn = zeros(3,length(t));
            %
            if insert_H_pump_manually
                H_pump_manual = load('C:\Users\owner\Desktop\Magnetization dynamics in the TLS picture\simulation results\alpha = 0.023\RCP\workspace.mat');
                %     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % !!! insert the name of the MATLAB file manually !!!
                %     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                H_pump_manual = (140)*(0.25*(10^-4))*H_pump_manual.ZZ_new_H_pump_manual; %H_pump_manual.ZZ_new_H_pump_manual;
                %
            end
            %
            if generate_Gaussian_H_pump
                %
                T_btw_Gs_pls = T_btw_Gs_pls_ratio*T_p; 
                %
                % The amplitude ratio between the x and y components of the pulse determine the polarization -
                % for x_amp_parameter = y_amp_parameter, the polarization can be circular
                x_amp_parameter = 1;  
                y_amp_parameter = 1;   
                %
                pol_aos_gs_pls_phase_rad  = pol_aos_gs_pls_phase*(2*pi/360); % Convert into radian units
                %
%                 H_pmp_sym_srt_time = miu_gs - T_p;
%                 H_pmp_sym_srt_ind = H_pmp_sym_srt_time/dt_low_res;
%                 H_pmp_sym_end_time = miu_gs + T_p;
%                 H_pmp_sym_end_ind = H_pmp_sym_end_time/dt_low_res;
                sigma_gs =  T_p; 
                % Create the Gaussian envelopes
                Gaussian_array_x = ((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs).^2)./(2*(sigma_gs^2)))));
                Gaussian_array_y = ((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs).^2)./(2*(sigma_gs^2)))));
                % For a single pulse, multiply the envelope by a cos function.
                if Num_Gaussian_pulses == 1
                    % Normalize the Gaussian envelopes
                    Gaussian_array_x = x_amp_parameter*Gaussian_array_x/(max(Gaussian_array_x));
                    Gaussian_array_y = y_amp_parameter*Gaussian_array_y/(max(Gaussian_array_y));
                    % Modulate the Gaussian envelope with the optical frequency
                    Gaussian_array_x = H_peak*cos(2*pi*freq_aos_gs_pls.*T_Global).*Gaussian_array_x;
                    Gaussian_array_y = H_peak*cos(2*pi*freq_aos_gs_pls.*T_Global - pol_aos_gs_pls_phase_rad).*Gaussian_array_y;
                    H_pump_manual = zeros(3,length(T_Global));
                    H_pump_manual(1,:) = Gaussian_array_x;
                    H_pump_manual(2,:) = Gaussian_array_y;
                end
                % For multiple pulses, create multiple Gaussian envelopes
                if ~random_phase_between_pulses
                    if Num_Gaussian_pulses > 1
                        for i=2:Num_Gaussian_pulses
                            Gaussian_array_x = Gaussian_array_x + ((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs-T_btw_Gs_pls*(i-1)).^2)./(2*(sigma_gs^2)))));
                            Gaussian_array_y = Gaussian_array_y + ((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs-T_btw_Gs_pls*(i-1)).^2)./(2*(sigma_gs^2)))));
                        end
                    end
                    % Normalize the Gaussian envelopes
                    Gaussian_array_x = x_amp_parameter*Gaussian_array_x/(max(Gaussian_array_x));
                    Gaussian_array_y = y_amp_parameter*Gaussian_array_y/(max(Gaussian_array_y));
                    % Modulate the Gaussian envelope with the optical frequency
                    Gaussian_array_x = H_peak*cos(2*pi*freq_aos_gs_pls.*T_Global).*Gaussian_array_x;
                    Gaussian_array_y = H_peak*cos(2*pi*freq_aos_gs_pls.*T_Global - pol_aos_gs_pls_phase_rad).*Gaussian_array_y;
                    H_pump_manual = zeros(3,length(T_Global));
                    H_pump_manual(1,:) = Gaussian_array_x;
                    H_pump_manual(2,:) = Gaussian_array_y;
                end
                %
                if random_phase_between_pulses
                    if Num_Gaussian_pulses > 1
                       random_phase = 2*pi*rand; % Generate a random phase for the first pulse.
                       Gaussian_array_x = cos(2*pi*freq_aos_gs_pls.*T_Global+random_phase).*Gaussian_array_x;
                       Gaussian_array_y = cos(2*pi*freq_aos_gs_pls.*T_Global - pol_aos_gs_pls_phase_rad+random_phase).*Gaussian_array_y;
                        for i=2:Num_Gaussian_pulses
                            random_phase = 2*pi*rand; % Generate a random phase for each pulse.
                            Gaussian_array_x = Gaussian_array_x + cos(2*pi*freq_aos_gs_pls.*T_Global+random_phase).*((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs-T_btw_Gs_pls*(i-1)).^2)./(2*(sigma_gs^2)))));
                            Gaussian_array_y = Gaussian_array_y + cos(2*pi*freq_aos_gs_pls.*T_Global - pol_aos_gs_pls_phase_rad+random_phase).*((1/(sqrt(2*pi)*sigma_gs)).*exp(-(((T_Global-miu_gs-T_btw_Gs_pls*(i-1)).^2)./(2*(sigma_gs^2)))));
                        end
                    end
                    % Normalize the Gaussian envelopes
                    Gaussian_array_x = x_amp_parameter*Gaussian_array_x/(max(Gaussian_array_x));
                    Gaussian_array_y = y_amp_parameter*Gaussian_array_y/(max(Gaussian_array_y));
                    % Scale the Gaussian envelope with the peak amplitude
                    Gaussian_array_x = H_peak.*Gaussian_array_x;
                    Gaussian_array_y = H_peak.*Gaussian_array_y;
                    H_pump_manual = zeros(3,length(T_Global));
                    H_pump_manual(1,:) = Gaussian_array_x;
                    H_pump_manual(2,:) = Gaussian_array_y;
                end
            end
            %
            if insert_S_ac_manually
                S_ac_manual = load('C:\Users\owner\Desktop\Magnetization dynamics in the TLS picture\simulation results\alpha=0.025 and lambda=800nm\LCP\workspace.mat');
                %     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % !!! insert the name of the MATLAB file manually !!!
                %     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                S_ac_manual = 0*(10^5)*S_ac_manual.ZZ_new_H_pump_manual; %H_pump_manual.ZZ_new_H_pump_manual;
                %
                S_ac = zeros(3,length(t));
                S_ac_an = zeros(3,length(t));
                S_ac_bn = zeros(3,length(t));
                S_ac_cn = zeros(3,length(t));
            end
            % Non - manual (generation of H_pump)
            % Be aware of time overlap between the three time section
            % Choose the start time such as it is longer than time section 2&3
            %
            H_pump_t_start = 1*(10^-14); %sec
            H_pump_t_end = 13*(10^-14); %sec
            %
            H_pump_t_start = H_pump_t_start-T_section{2}(end);%0.2*(130.5*(10^-9) - T_section{2}(end)); %sec
            H_pump_t_end = H_pump_t_end-T_section{2}(end);%- 2*10^-14; %- T_section{2}(end); % sec
            %
            n_min = round(H_pump_t_start/dt);
            n_max = round(H_pump_t_end/dt);
            H_pump_on_off = 0;
            IFE_field_amplitude = 0*4*(10^-6)*(n_id*c*eps)/2;%0.8*(10^-7); % creates a field according to the IFE
            %
            % # of time samples in each rect period of H_pump
            H_pump_on_off_time_period = 10*(10^-17); % sec
            H_pump_on_off_period = round(H_pump_on_off_time_period/dt);
            %
            % Calculate the maximal amplitude of H_pump (Tesla\miu0) in the
            % run
            calc_max_pump = 1;
            H_pump_max = 0;
            %
            % Give an absolute maximal value that H_pump might attain
            H_pump_chop = 0*1;
            H_pump_chop_val = 0.1*H_Xtrnl(3);         
            %
            % Reverse the sign of H_pump after Mz reaches a certain
            % absolute value
            % This will result in a pump in the opposite direction
            H_pump_change_sign = 1;
            Mz_val_for_H_pump_sgn_change = 0.99*M0;
            %
            % Turn off H_pump after Mz reaches a certain
            % absolute value
            % This will result in a single switching in each time section
            H_pump_turn_off = 0;
            Mz_val_for_H_pump_turn_off = 0*M0;
            % mark the switching point as the 'pulse time'
            % only for a simulation with no pulses
            set_t_pls_as_switching_point = 0;
            if Section_ind == 1
               % internal value
               % DONT CHANGE
               turn_off_pump = 0;
               % inner parameters of the sign cahnge implementation of H_pump
               % DONT CHANGE
               H_pump_inner_sign = 1;
               inner_sign_count = 0; 
               H_pump_count_sign = 0;
               %
               % inner parameters of the on-off implementation of H_pump
               % DONT CHANGE
               H_pump_count = 0;
               calc_H_pump = 1;
            end
            %%
            %% Simulation loop
            for n = 1:(length(t)-1)
                Ht_an=H_Xtrnl+(2*K0_0(n)/(A_0(n)*M0^4))*[-(M_tg(3,n))^2 -(M_tg(3,n))^2 ((M_tg(1,n))^2+(M_tg(2,n))^2)]'.*M_tg(:,n)...
                    -N_Dmag.*M_tg(:,n)+H_CW(:,n)+H_Pulsed_0(:,n)+H_pump(:,n)+H_IFE(:,n);      %Ht is the total field
                an=(1/A_Global)*(-gma_LL_0(n)*miu*cross(M_tg(:,n),Ht_an)-(LL_lada_0(n)*miu/M0)*cross(M_tg(:,n),cross(M_tg(:,n),Ht_an))-(1/A_0(n))*dA_0(n)*M_tg(:,n)...
                    +(gma*Hshe_ac/M0)*cross(M_tg(:,n),cross(M_tg(:,n),S_ac(:,n)))+(gma*Hshe_dc/M0)*cross(M_tg(:,n),cross(M_tg(:,n),S_vector)));
                Ht_bn=H_Xtrnl+...
                    (2*K0_an(n)/(A_an(n)*M0^4))*[-(M_tg(3,n)+(dt/2)*an(3))^2 -(M_tg(3,n)+(dt/2)*an(3))^2 ((M_tg(1,n)+(dt/2)*an(1))^2+(M_tg(2,n)+(dt/2)*an(2))^2)]'.*(M_tg(:,n)+(dt/2)*an)...
                    -N_Dmag.*(M_tg(:,n)+(dt/2)*an)+H_CW_an(:,n)+H_Pulsed_an(:,n)+H_pump_an(:,n)+H_IFE_an(:,n);
                bn=(1/A_Global)*(-gma_LL_an(n)*miu*cross(M_tg(:,n)+(dt/2)*an,Ht_bn)-(LL_lada_an(n)*miu/M0)*cross(M_tg(:,n)+(dt/2)*an,cross(M_tg(:,n)+(dt/2)*an,Ht_bn))-(1/A_an(n))*dA_an(n)*(M_tg(:,n)+(dt/2)*an)...
                    +(gma*Hshe_ac/M0)*cross(M_tg(:,n)+(dt/2)*an,cross(M_tg(:,n)+(dt/2)*an,S_ac_an(:,n)))+(gma*Hshe_dc/M0)*cross(M_tg(:,n)+(dt/2)*an,cross(M_tg(:,n)+(dt/2)*an,S_vector)));
                
                Ht_cn=H_Xtrnl+...
                    (2*K0_bn(n)/(A_bn(n)*M0^4))*[-(M_tg(3,n)+(dt/2)*bn(3))^2 -(M_tg(3,n)+(dt/2)*bn(3))^2 ((M_tg(1,n)+(dt/2)*bn(1))^2+(M_tg(2,n)+(dt/2)*bn(2))^2)]'.*(M_tg(:,n)+(dt/2)*bn)...
                    -N_Dmag.*(M_tg(:,n)+(dt/2)*bn)+H_CW_bn(:,n)+H_Pulsed_bn(:,n)+H_pump_bn(:,n)+H_IFE_bn(:,n);
                cn=(1/A_Global)*(-gma_LL_bn(n)*miu*cross(M_tg(:,n)+(dt/2)*bn,Ht_cn)-(LL_lada_bn(n)*miu/M0)*cross(M_tg(:,n)+(dt/2)*bn,cross(M_tg(:,n)+(dt/2)*bn,Ht_cn))-(1/A_bn(n))*dA_bn(n)*(M_tg(:,n)+(dt/2)*bn)...
                    +(gma*Hshe_ac/M0)*cross(M_tg(:,n)+(dt/2)*bn,cross(M_tg(:,n)+(dt/2)*bn,S_ac_bn(:,n)))+(gma*Hshe_dc/M0)*cross(M_tg(:,n)+(dt/2)*bn,cross(M_tg(:,n)+(dt/2)*bn,S_vector)));
                
                Ht_dn=H_Xtrnl+...
                    (2*K0_cn(n)/(A_cn(n)*M0^4))*[-(M_tg(3,n)+dt*cn(3))^2 -(M_tg(3,n)+dt*cn(3))^2 ((M_tg(1,n)+dt*cn(1))^2+(M_tg(2,n)+dt*cn(2))^2)]'.*(M_tg(:,n)+dt*cn)...
                    -N_Dmag.*(M_tg(:,n)+dt*cn)+H_CW_cn(:,n)+H_Pulsed_cn(:,n)+H_pump_cn(:,n)+H_IFE_cn(:,n);
                dn=(1/A_Global)*(-gma_LL_cn(n)*miu*cross(M_tg(:,n)+dt*cn,Ht_dn)-(LL_lada_cn(n)*miu/M0)*cross(M_tg(:,n)+dt*cn,cross(M_tg(:,n)+dt*cn,Ht_dn))-(1/A_cn(n))*dA_cn(n)*(M_tg(:,n)+dt*cn)...
                    +(gma*Hshe_ac/M0)*cross(M_tg(:,n)+dt*cn,cross(M_tg(:,n)+dt*cn,S_ac_cn(:,n)))+(gma*Hshe_dc/M0)*cross(M_tg(:,n)+dt*cn,cross(M_tg(:,n)+dt*cn,S_vector)));
                
                M_tg(:,n+1) = M_tg(:,n)+(dt/6)*(an+2*bn+2*cn+dn);
                
                %% S_ac response to dM/dt - A_feedback_xy/z
                % Additional response of S_ac to change in magnetization
                % This response is affected only by the x,y components
                % The dipole moment in the 'spin picture' is proportional 
                % to mz, and thus we will multiply the interaction term by
                % mz (normalized).
                % 
               if insert_S_ac_manually
                   S_ac(:,n+1) = S_ac_manual(:,round((T_section{2}(end))/dt)+n+1);
                   S_ac_an(:,n+1) = S_ac_manual(:,round((T_section{2}(end))/dt)+n+1);
                   S_ac_bn(:,n+1) = S_ac_manual(:,round((T_section{2}(end))/dt)+n+1);
                   S_ac_cn(:,n+1) = S_ac_manual(:,round((T_section{2}(end))/dt)+n+1);
               end  
                %
                % Save the initial value of S_ac, consisting only of the 
                % extrenal AC excitation.
                %
                S_ac_initial(:,n+1) = S_ac(:,n+1);
                %
                norm_M = normalize(M_tg(:,n),'norm');
                % To nullify the dipole moment mz, set norm_M(3)=1;
                % norm_M(3) = 1;
                if Im_SP_xy ~= 0
                    S_ac(1,n+1) = S_ac(1,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(1)+2*bn(1)+2*cn(1)+dn(1));
                    S_ac(2,n+1) = S_ac(2,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(2)+2*bn(2)+2*cn(2)+dn(2));
                    S_ac_an(1,n+1) = S_ac_an(1,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(1)+2*bn(1)+2*cn(1)+dn(1));
                    S_ac_an(2,n+1) = S_ac_an(2,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(2)+2*bn(2)+2*cn(2)+dn(2));
                    S_ac_bn(1,n+1) = S_ac_bn(1,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(1)+2*bn(1)+2*cn(1)+dn(1));
                    S_ac_bn(2,n+1) = S_ac_bn(2,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(2)+2*bn(2)+2*cn(2)+dn(2));
                    S_ac_cn(1,n+1) = S_ac_cn(1,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(1)+2*bn(1)+2*cn(1)+dn(1));
                    S_ac_cn(2,n+1) = S_ac_cn(2,n+1) + norm_M(3)*Im_SP_xy*(1/6)*(an(2)+2*bn(2)+2*cn(2)+dn(2));
                end
                %
                % Addition of the z component to the feedback
                % to nullify it put it as a comment
                if Im_SP_z ~= 0
                    S_ac(3,n+1) = S_ac(3,n+1) + norm_M(3)*Im_SP_z*(1/6)*(an(3)+2*bn(3)+2*cn(3)+dn(3));
                    S_ac_an(3,n+1) = S_ac_an(3,n+1) + norm_M(3)*Im_SP_z*(1/6)*(an(3)+2*bn(3)+2*cn(3)+dn(3));
                    S_ac_bn(3,n+1) = S_ac_bn(3,n+1) + norm_M(3)*Im_SP_z*(1/6)*(an(3)+2*bn(3)+2*cn(3)+dn(3));
                    S_ac_cn(3,n+1) = S_ac_cn(3,n+1) + norm_M(3)*Im_SP_z*(1/6)*(an(3)+2*bn(3)+2*cn(3)+dn(3));
                end
                %
                % M(:,n+1)=M(:,n+1)*M0/sqrt(M(:,n+1)'*M(:,n+1));  % to preserve the size of the magnetization - 
                % correct error according to Carlos J. Garcia Cervera (Numerical micromagnetics: A Review)
                %
                %% S_ac response to (M x dM/dt) - A_feedback_xy/z_SP
                % Additional response of S_ac to change in magnetization
                % This response is affected only by the x,y components
                % The dipole moment in the 'spin picture' is proportional 
                % to mz, and thus we will multiply the interaction term by
                % mz (normalized)
                %
                norm_M = normalize(M_tg(:,n),'norm');
                % To nullify the dipole moment mz, set norm_M(3)=1;
                norm_M(3) = 1;
                % Define the (M x dM/dt) vector, by using the dM/dt 
                % just calculated, and a mean value for M, using the 
                % time steps n and n + 1
                S_sp(:,n+1) = (1/M0^2)*(1/6)*cross(0.5*(M_tg(:,n+1) + M_tg(:,n)),(an+2*bn+2*cn+dn));
                S_sp_temp = S_sp(:,n+1);
                %
                if SP_closed_loop
                % Writtern as a vector form, for an equal effect on x,y,z:
                % Multiply by (1/(Hshe_ac)) in order to cancel it in
                % LLGS, taking into account only {A_feedback_SP*gma/M0} and the
                % backscatter ratio
                % Hshe_ac=(hbar*Theta_sh*Jc_ac)/(2*q*M0*Thick);  %spin Hall effect field parameter 

                Tot_Re_SP_AC = sign_J_s_J_sp*SP_efficiency_ratio*(1/(Hshe_ac))*norm_M(3)*Re_SP;
                    if Re_SP ~= 0
                        S_ac(:,n+1) = S_ac(:,n+1) + Tot_Re_SP_AC*S_sp_temp;
                        S_ac_an(:,n+1) = S_ac_an(:,n+1) + Tot_Re_SP_AC*S_sp_temp;
                        S_ac_bn(:,n+1) = S_ac_bn(:,n+1) + Tot_Re_SP_AC*S_sp_temp;
                        S_ac_cn(:,n+1) = S_ac_cn(:,n+1) + Tot_Re_SP_AC*S_sp_temp;
                    end
                    %
                    % Writtern in a component form, for a NON-equal effect on x,y,z.
                    % To nullify, set A_feedback_xy/z_SP = 0
                    %
                   if  Re_SP_xy ~= 0 
                        S_ac(1,n+1) = S_ac(1,n+1) + Tot_Re_SP_AC*S_sp_temp(1);
                        S_ac(2,n+1) = S_ac(2,n+1) + Tot_Re_SP_AC*S_sp_temp(2);
                        S_ac_an(1,n+1) = S_ac_an(1,n+1) + Tot_Re_SP_AC*S_sp_temp(1);
                        S_ac_an(2,n+1) = S_ac_an(2,n+1) + Tot_Re_SP_AC*S_sp_temp(2);
                        S_ac_bn(1,n+1) = S_ac_bn(1,n+1) + Tot_Re_SP_AC*S_sp_temp(1);
                        S_ac_bn(2,n+1) = S_ac_bn(2,n+1) + Tot_Re_SP_AC*S_sp_temp(2);
                        S_ac_cn(1,n+1) = S_ac_cn(1,n+1) + Tot_Re_SP_AC*S_sp_temp(1);
                        S_ac_cn(2,n+1) = S_ac_cn(2,n+1) + Tot_Re_SP_AC*S_sp_temp(2);            
                   end
                      %
                      % Addition of the z component to the feedback
                    if  Re_SP_z ~= 0 
                        S_ac(3,n+1) = S_ac(3,n+1) + Tot_Re_SP_AC*S_sp_temp(3);
                        S_ac_an(3,n+1) = S_ac_an(3,n+1) + Tot_Re_SP_AC*S_sp_temp(3);
                        S_ac_bn(3,n+1) = S_ac_bn(3,n+1) + Tot_Re_SP_AC*S_sp_temp(3);
                        S_ac_cn(3,n+1) = S_ac_cn(3,n+1) + Tot_Re_SP_AC*S_sp_temp(3);
                    end
                end
                %
                %M(:,n+1)=M(:,n+1)*M0/sqrt(M(:,n+1)'*M(:,n+1));  % to preserve the size of the magnetization - 
                %correct error according to Carlos J. Garcia Cervera (Numerical micromagnetics: A Review)
                %
                % H_pump field
                %
                if ~ insert_H_pump_manually && ~generate_Gaussian_H_pump && ~turn_off_pump && n>n_min && n<n_max
                    %
                    if H_pump_turn_off && (0.5*(M_tg(3,n+1)+M_tg(3,n))) < 0 && abs(0.5*(M_tg(3,n+1)+M_tg(3,n))) > abs(Mz_val_for_H_pump_turn_off) 
                        n_min = n_max + 1;
                        if set_t_pls_as_switching_point
                            t0_H_pulsed = t(n);
%                             t0_H_pulsed = t(n) + (Section_ind > 1)*T_section{Section_ind - 1}(end);
                            T_int1 = round(T_int_fraction*(t0_H_pulsed/dt));
                            T_int2 = abs(round(T_int_fraction*((t0_H_pulsed_end-t0_H_pulsed)/dt)));
                        end
                        turn_off_pump = ~turn_off_pump;
                    end
                    %
                    if H_pump_on_off == 1
                        H_pump_count = H_pump_count + 1;
                        if H_pump_count == H_pump_on_off_period
                           H_pump_count = 0;
                           calc_H_pump = ~ calc_H_pump;
                        end
                        if calc_H_pump
                            H_pump(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n))*cos(H_pump_phase_rad)-(M_tg(2,n))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n)))^2);
                            H_pump(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n))*sin(H_pump_phase_rad)+(M_tg(2,n))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n)))^2);
                            H_pump_an(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*an(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+(dt/2)*an(2))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*an(3))^2);
                            H_pump_an(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*an(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+(dt/2)*an(2))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*an(3))^2);
                            H_pump_bn(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*bn(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+(dt/2)*bn(2))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*bn(3))^2);
                            H_pump_bn(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*bn(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+(dt/2)*bn(2))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*bn(3))^2);
                            H_pump_cn(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+dt*cn(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+dt*cn(1))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+dt*cn(1))^2);
                            H_pump_cn(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+dt*cn(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+dt*cn(1))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+dt*cn(1))^2);
                        end
                    end
                    %
                    if H_pump_on_off == 0
                            H_pump(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n))*cos(H_pump_phase_rad)-(M_tg(2,n))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n)))^2);
                            H_pump(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n))*sin(H_pump_phase_rad)+(M_tg(2,n))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n)))^2);
                            H_pump_an(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*an(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+(dt/2)*an(2))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*an(3))^2);
                            H_pump_an(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*an(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+(dt/2)*an(2))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*an(3))^2);
                            H_pump_bn(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*bn(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+(dt/2)*bn(2))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*bn(3))^2);
                            H_pump_bn(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+(dt/2)*bn(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+(dt/2)*bn(2))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+(dt/2)*bn(3))^2);
                            H_pump_cn(1,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+dt*cn(1))*cos(H_pump_phase_rad)-(M_tg(2,n)+dt*cn(2))*sin(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+dt*cn(3))^2);
                            H_pump_cn(2,n+1) = H_pump_inner_sign*H_pump_0*((M_tg(1,n)+dt*cn(1))*sin(H_pump_phase_rad)+(M_tg(2,n)+dt*cn(2))*cos(H_pump_phase_rad))/(M0^2 - ((M_tg(3,n))+dt*cn(3))^2);
                            %
                            H_IFE(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*((H_pump_0)^2)/(M0^2 - ((M_tg(3,n)))^2);
                            H_IFE_an(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*((H_pump_0)^2)/(M0^2 - ((M_tg(3,n))+(dt/2)*an(3))^2);
                            H_IFE_bn(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*((H_pump_0)^2)/(M0^2 - ((M_tg(3,n))+(dt/2)*bn(3))^2);
                            H_IFE_cn(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*((H_pump_0)^2)/(M0^2 - ((M_tg(3,n))+dt*cn(3))^2);
                    end
                    %
                    if H_pump_chop
                        if  abs(H_pump(2,n+1)) > H_pump_chop_val
                            sgn_H_pump_y = sign (H_pump(2,n+1));
                            H_pump(2,n+1) = sgn_H_pump_y*H_pump_chop_val;
                        end
                        %
                        if  abs(H_pump(1,n+1)) > H_pump_chop_val
                            sgn_H_pump_x = sign (H_pump(1,n+1));
                            H_pump(1,n+1) = sgn_H_pump_x*H_pump_chop_val;
                        end
                    end
                    if calc_max_pump
                    %
                        if abs(H_pump(2,n+1)) > H_pump_max
                            H_pump_max = abs(H_pump(2,n+1));
                        end
                        %
                        if abs(H_pump(1,n+1)) > H_pump_max
                            H_pump_max = abs(H_pump(1,n+1));
                        end
                    end
                    %
                    if inner_sign_count 
                        H_pump_count_sign = H_pump_count_sign + 1;
                        if H_pump_count_sign == H_pump_on_off_period
                           H_pump_count_sign = 0;
                           inner_sign_count = ~ inner_sign_count;
                        end
                    end
                    %
                    if ~inner_sign_count && H_pump_change_sign && abs(0.5*(M_tg(3,n+1)+M_tg(3,n))) > abs(Mz_val_for_H_pump_sgn_change) 
                       H_pump_inner_sign = (-1)*H_pump_inner_sign;
                       inner_sign_count = 1;
                    end
                    %
                end
                %
                if insert_H_pump_manually
                    H_pump(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_an(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_bn(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_cn(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_IFE(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*((H_pump(1,n+1)^2)+(H_pump(2,n+1)^2));
                    H_IFE_an(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*(((H_pump(1,n+1)+(dt/2)*an(1))^2)+((H_pump(2,n+1)+(dt/2)*an(2))^2));
                    H_IFE_bn(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*(((H_pump(1,n+1)+(dt/2)*bn(1))^2)+((H_pump(2,n+1)+(dt/2)*bn(2))^2));
                    H_IFE_cn(3,n+1) = H_pump_inner_sign*IFE_field_amplitude*(((H_pump(1,n+1)+(dt)*cn(1))^2)+((H_pump(2,n+1)+(dt)*cn(2))^2));
                end
                % 
                if generate_Gaussian_H_pump
                    H_pump(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_an(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_bn(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                    H_pump_cn(:,n+1) = H_pump_manual(:,round((T_section{2}(end))/dt)+n+1);
                end
                % Calculate TLS terms during run
                gamma_v(:,n+1) = -gma*miu*alpha*((M0*(1+(alpha^2)))^-1)*((H_Xtrnl_vect*(1-STT_DC_ratio)+(0.5*(M_tg(3,n+1)+M_tg(3,n)))*N_Dmag(3))*(M0-0.5*(M_tg(3,n+1)+M_tg(3,n)))-H_pump_0*cos(H_pump_phase_rad)*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max));
                gamma_c(:,n+1) = gma*miu*alpha*((M0*(1+(alpha^2)))^-1)*((H_Xtrnl_vect*(1-STT_DC_ratio)+(0.5*(M_tg(3,n+1)+M_tg(3,n)))*N_Dmag(3))*(M0+0.5*(M_tg(3,n+1)+M_tg(3,n)))+H_pump_0*cos(H_pump_phase_rad)*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max));
                ReV12(:,n+1) = gma*miu*((2*(1+(alpha^2)))^-1)*((0.5*(M_tg(1,n+1)+M_tg(1,n)))*N_Dmag(1)-(alpha/M0)*((0.5*(M_tg(3,n+1)+M_tg(3,n))))*((0.5*(M_tg(2,n+1)+M_tg(2,n)))*N_Dmag(2)+sig_rf(n)));
                ImV12(:,n+1) = -gma*miu*((2*(1+(alpha^2)))^-1)*((alpha/M0)*((0.5*(M_tg(3,n+1)+M_tg(3,n))))*(0.5*(M_tg(1,n+1)+M_tg(1,n)))*N_Dmag(1)+(0.5*(M_tg(2,n+1)+M_tg(2,n)))*N_Dmag(2)+sig_rf(n))+(gma/(2*M0))*0.5*(M_tg(3,n+1)+M_tg(3,n))*Hshe_ac*S_ac_initial(1,n+1);
                Lambda_h(:,n+1) = gma*miu*H_pump_0*((1+(alpha^2))^-1)*(sin(H_pump_phase_rad)-alpha*cos(H_pump_phase_rad))*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max);
                Lambda_e(:,n+1) = gma*miu*H_pump_0*((1+(alpha^2))^-1)*(sin(H_pump_phase_rad)+alpha*cos(H_pump_phase_rad))*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max);
                w_TLS(:,n+1) = gma*miu*((1+(alpha^2))^-1)*(H_Xtrnl_vect+0.5*(M_tg(3,n+1)+M_tg(3,n))*N_Dmag(3)+(1+(alpha)^2)*(0.5/(M0*miu))*(M_tg(2,n+1)+M_tg(2,n))*Hshe_ac*S_ac_initial(1,n+1)-((0.5*(M_tg(3,n+1)+M_tg(3,n)))*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max)*H_pump_0*cos(H_pump_phase_rad)/(M0^2 - (0.5*(M_tg(3,n+1)+M_tg(3,n)))^2))+(alpha/M0)*((((M0)^2)*(~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max)*H_pump_0*sin(H_pump_phase_rad)/(M0^2 - (0.5*(M_tg(3,n+1)+M_tg(3,n)))^2))-((N_Dmag(1)-N_Dmag(2))*(0.5*(M_tg(1,n+1)+M_tg(1,n))))*0.5*((M_tg(2,n+1)+M_tg(2,n)))+(0.5*(M_tg(1,n+1)+M_tg(1,n)))*sig_rf(n)));
                gamma_inh(:,n+1) = gma*miu*((1+(alpha^2))^-1)*((0.5*(M_tg(3,n+1)+M_tg(3,n)))/M0)*(((~ insert_H_pump_manually && ~turn_off_pump && n>n_min && n<n_max)*H_pump_0/(M0^2 - (0.5*(M_tg(3,n+1)+M_tg(3,n)))^2))*(alpha*cos(H_pump_phase_rad)*0.5*(M_tg(3,n+1)+M_tg(3,n))-M0*sin(H_pump_phase_rad))+alpha*(H_Xtrnl_vect*(1-STT_DC_ratio)+(0.5*(M_tg(3,n+1)+M_tg(3,n)))*N_Dmag(3)));
                %
                H_Tot(:,n+1)=Ht_an; % Record the total field from which Anis us recorded as well by subtraction of H_anis
                H_Anis(:,n+1)=Ht_an-H_Xtrnl-H_CW(:,n); % There is a very slight misaccuracy of these recorded vectors - thy are shifted in time by one index
                %

                primary_torque(:,n+1) = (1/A_Global)*(-gma_LL_0(n)*miu*cross(M_tg(:,n),Ht_an));
                damping_torque(:,n+1) = (1/A_Global)*(-LL_lada_0(n)*miu/M0)*cross(M_tg(:,n),cross(M_tg(:,n),Ht_an));
                %
            end
            % 
            t0_between_sections = t0;
            t0_H_pulsed_between_sections = t0_H_pulsed;
            T_int1_between_sections = T_int1;
            T_int2_between_sections = T_int2;
            %
            % Save the last value of S_ac
            S_ac_between_sections = S_ac(:,n+1);
            S_ac_between_sections_an = S_ac_an(:,n+1);
            S_ac_between_sections_bn = S_ac_bn(:,n+1);
            S_ac_between_sections_cn = S_ac_cn(:,n+1);
            % Save the last value of the TLS terms
            gamma_v_between_sections = gamma_v(:,n+1);
            gamma_c_between_sections = gamma_c(:,n+1);
            ReV12_between_sections = ReV12(:,n+1);
            ImV12_between_sections = ImV12(:,n+1);
            Lambda_h_between_sections = Lambda_h(:,n+1);
            Lambda_e_between_sections = Lambda_e(:,n+1);
            w_TLS_between_sections = w_TLS(:,n+1);
            gamma_inh_between_sections = gamma_inh(:,n+1);
            primary_torque_between_sections = primary_torque(:,n+1);
            damping_torque_between_sections = damping_torque(:,n+1);
            %
            % Enlarge the arrays by the calculated time section
            M_array_global{L} = [M_array_global{L}(:,1:end-1) M_tg];
            S_ac_array_global{L} = [S_ac_array_global{L}(:,1:end-1) S_ac];
            S_ac_array_initial{L} = [S_ac_array_initial{L}(:,1:end-1) S_ac_initial];
            S_sp_array_global{L} = [S_sp_array_global{L}(:,1:end-1) Re_SP*S_sp];
            H_CW_array_global{L} = [H_CW_array_global{L}(:,1:end-1) H_CW];
            H_pump_array_global{L} = [H_pump_array_global{L}(:,1:end-1) H_pump];
            H_IFE_global{L} = [H_IFE_global{L}(:,1:end-1) H_IFE];
            gamma_v_global{L} = [gamma_v_global{L}(:,1:end-1) gamma_v];
            gamma_c_global{L} = [gamma_c_global{L}(:,1:end-1) gamma_c];
            ReV12_global{L} = [ReV12_global{L}(:,1:end-1) ReV12];
            ImV12_global{L} = [ImV12_global{L}(:,1:end-1) ImV12];
            Lambda_h_global{L} = [Lambda_h_global{L}(:,1:end-1) Lambda_h];
            Lambda_e_global{L} = [Lambda_e_global{L}(:,1:end-1) Lambda_e];
            w_TLS_global{L} = [w_TLS_global{L}(:,1:end-1) w_TLS];
            gamma_inh_global{L} = [gamma_inh_global{L}(:,1:end-1) gamma_inh];
            primary_torque_array_global{L} = [primary_torque_array_global{L}(:,1:end-1) primary_torque];
            damping_torque_array_global{L} = [damping_torque_array_global{L}(:,1:end-1) damping_torque];
            %
        end                    % End of Time Loops
        %
        % Update the arrays at the end of the loop
        M_array{L} = M_array_global{L};
        S_ac_array{L} = S_ac_array_global{L};
        S_ac_initial_array{L} = S_ac_array_initial{L};
        S_sp_array{L} = S_sp_array_global{L};
        H_pump_array{L} = H_pump_array_global{L};
        H_IFE_array{L} = H_IFE_global{L};
        gamma_v_array{L} = gamma_v_global{L};
        gamma_c_array{L} = gamma_c_global{L};
        ReV12_array{L} = ReV12_global{L};
        ImV12_array{L} = ImV12_global{L};
        Lambda_h_array{L} = Lambda_h_global{L};
        Lambda_e_array{L} = Lambda_e_global{L};
        w_TLS_array{L} = w_TLS_global{L};
        gamma_inh_array{L} = gamma_inh_global{L};
        primary_torque_array{L} = primary_torque_array_global{L};
        damping_torque_array{L} = damping_torque_array_global{L};
        %
        theta_01(L)=(acos(M_array_global{L}(3,1)/M0))*180/pi;
        phi_01(L)=(atan(M_array_global{L}(2,1)/M_array_global{L}(1,1)))*180/pi;

        theta_0_end(L)=(acos(M_array_global{L}(3,end)/M0))*180/pi;
        phi_0_end(L)=(atan(M_array_global{L}(2,end)/M_array_global{L}(1,end)))*180/pi;
       
        a=cos(w*T_Global);
        
        b=sin(w*T_Global);
        
        % Components of Rotated Magnetization 
        % The 'rotated' part (desired to be removed) is the trivial 
        % sinusoidal oscillation of the transverse (xy) compenents
            
        M_array_Rot_x{L}(1,:)= M_array_global{L}(1,:).*a+M_array_global{L}(2,:).*b;   
        M_array_Rot_y{L}(1,:)= -M_array_global{L}(1,:).*b+M_array_global{L}(2,:).*a;
        M_array_Rot_z{L}(1,:) = M_array_global{L}(3,:);    

        % Starting and Ending Points of Rotated Magnetization
        
        theta_01_Rot(L)=(acos(M_array_Rot_z{L}(1,1)/M0))*180/pi;
        
        [azimuth_01(L),elevation_01(L),r_01(L)]=cart2sph(M_array_Rot_x{L}(1,1),M_array_Rot_y{L}(1,1),M_array_Rot_z{L}(1,1));
        
        phi_01_Rot(L)=(azimuth_01(L))*180/pi;
        
%         phi_01_Rot(L)=(atan(M_array_Rot_y{L}(1,1)/M_array_Rot_x{L}(1,1)))*180/pi;

        theta_0_end_Rot(L)=(acos(M_array_Rot_z{L}(1,end)/M0))*180/pi;
        
        [azimuth_end(L),elevation_end(L),r_end(L)]=cart2sph(M_array_Rot_x{L}(1,end),M_array_Rot_y{L}(1,end),M_array_Rot_z{L}(1,end));
        
        phi_0_end_Rot(L)=(azimuth_end(L))*180/pi;
        
%         phi_0_end_Rot(L)=(atan(M_array_Rot_y{L}(1,end)/M_array_Rot_x{L}(1,end)))*180/pi;


        if 0                   %plot last calculated magnetization
            figure;
            plot(T_Global,M_array_global{L}(2,:),...
                T_section{1}(1:end-1),M_array_global{L}(2,1:length(T_section{1})-1),...
                T_section{3},M_array_global{L}(2,end-length(T_section{3})+1:end));
        end
        
    end                       %end of L Loops
    close(h)
    
    %% Plots 
    % Define font for figures
    ax_num_font_size = 30; % font of numers in x\y axes
    axis_font = 35; % font of x\y axes
    title_font = 40; % font of titles
    smaller_font_diff = 9;
    lgd_font = 20;
    %
    if disp_time_femto
       t_correction = 10^15;
       f_correction = 10^-15;
    end
    %
    if disp_time_pico
       t_correction = 10^12;
       f_correction = 10^-12;
    end
    %
    if disp_time_nano
       t_correction = 10^9;
       f_correction = 10^-9;
    end
    T_Global = T_Global*t_correction;
    %
    % Fourier Transform variables:
    %
    % Define row vector of zeros for padding before fourier transform in
    % order to increase the frequency resolution
    zero_pad_array_for_fourier = zeros(1,5*length(T_Global));
    %
    F_range_for_fourier_plots = 0.5; % in [PHz] unit, used for the anlysis of optical pulses in LLG
    enlarge_f_range = 15;
    % Define frequency range around the central CW frequency to be zoomed
    % in fourier plots
    Freq_CW_fraction_for_fourier_plot = 0.04*Freq_CW;
    %
    % (0) Plots all the relevant parameters of the simulation
    plot_param = 0*1;
    plot_param = plot_param*do_plot;
    %
    % (1) preferable for a single B value analysis
    % To plot the normalized Mx,My,Mz in different figures, displaying
    % both the envelope and oscillations, set plot_M_env_osc = 1
    plot_M_env_osc = 0;
    plot_M_env_osc = plot_M_env_osc*do_plot;
    %
    % (2) preferable for a single B value analysis
    % To plot the normalized  Sx_ac,Sy_ac,Sz_ac in different figures,
    % displaying both the envelope and oscillations, set plot_S_ac_env_osc = 1
    plot_S_ac_env_osc = 0;
    plot_S_ac_env_osc = plot_S_ac_env_osc*do_plot;
    %
    % (3)
    % To plot S_ac_x,S_ac_y,S_ac_z in the same figure,
    % set plot_S_ac_legend_not_norm = 1
    % Not normalized S_ac is useful in comparing the different 
    % components of S_ac in the same plot
    plot_S_ac_legend_not_norm = 0*1;
    plot_S_ac_legend_not_norm = plot_S_ac_legend_not_norm*do_plot;
    % To include the initial and final S_ac_X,
    % set view_initial_and_final_S_ac = 1
    view_initial_and_final_S_ac = 1;
    %
    % (4)
    % In the OPEN loop of spin pumping (i.e. the SP current only exits the 
    % FM and doesn't interact), plot the SP current.
    % In the CLOSED loop, a fraction of the spin pumping will be
    % backscattered, enlarging the Gilbert damping.
    plot_S_sp_legend_not_norm = 0*1;
    plot_S_sp_legend_not_norm = plot_S_sp_legend_not_norm*do_plot;
    %
    %(5)
    % To plot the sum of S_sp and S_initial, which is supposed to be
    % measured at the out, set plot_total_S_sp_initial = 1;
    plot_total_S_sp_initial = 0*1;
    plot_total_S_sp_initial = plot_total_S_sp_initial*do_plot;
    %
    % In order to view only the absorption\amplification, only the x
    % component is required, because the initial spin current is
    % x-polarized, such that S_sp y,z are generated, not
    % absorbed\amplified.
    tot_sp_int_only_x = 1;
    %
    % Plot the phase between J_sp and J_ac.
    % In case of pulses, set plot_phase_entire_time = 0,
    % dividing the phase plots to before and during the pulses.
    % NOTE: use only for dt equal for all time sections!
    plot_phase_entire_time = 0;
    %
    % (6)
    % To plot the normalized Mx,My,Mz in the same figure, displaying only 
    % the envelope (using the rotated magnetization) set plot_M_env = 1
    plot_M_env = 0;
    plot_M_env = plot_M_env*do_plot;
    %
    % (7)
    % To plot Mx,My,Mz in the same figure, set plot_M_env = 1
    % Not normalized magnetization is useful in comparing the different 
    % components of the magnetization in the same plot
    %
    % displaying only the envelope(using the rotated magnetization):
    plot_M_legend_not_norm_rot = 0;
    plot_M_legend_not_norm_rot = plot_M_legend_not_norm_rot*do_plot;
    %
    % displaying the not-rotated magnetization:
    plot_M_legend_not_norm = 1;
    plot_M_legend_not_norm = plot_M_legend_not_norm*do_plot;
    %
    % (8) only for multiple B value analysis
    % To plot all Mz as function of time with their respective B,
    % set plot_Mz_legend_offset = 1
    % This is usful in order to see the resonace nature of B in LLGS
    plot_Mz_legend_offset = 0;
    plot_Mz_legend_offset = plot_Mz_legend_offset*do_plot;
    %
    % (9)
    % To plot three figures of Mx,My,Mz and their respective power
    % spectrums, set plot_M_with_Power_Spectrum = 1
    plot_M_with_Power_Spectrum = 1;
    plot_Mx_My_filtered = 0*1; % Additionally, it is possible to filter out the optical carrier from Mx and My
    plot_M_with_Power_Spectrum = plot_M_with_Power_Spectrum*do_plot;
    plot_Mx_My_filtered = plot_Mx_My_filtered*do_plot;
    %
    % (10)
    % To plot two figures of H_pump_x,H_pump_y and their respective power
    % spectrums, set plot_H_pump = 1
    plot_H_pump = 1;
    plot_H_pump = plot_H_pump*do_plot;
    plot_zoomed_H_pump_spectrum = 0;
    %
    plot_H_pump_and_IFE = 0*1;
    plot_H_pump_and_IFE = plot_H_pump_and_IFE*do_plot;
    %
    % (11)
    % To plot the primary torque induced by a field (no anis and stt), and 
    % the LL damping torque, set plot_surf_norm_extr_M = 1;
    plot_toruqes = 0*1;
    plot_toruqes = plot_toruqes*do_plot;
    %
    % (12)
    % only for multiple B value analysis
    % To plot  M as function of time and external DC magnetic 
    % field B, set plot_surf_norm_extr_M = 1;
    plot_surf_M = 0;
    plot_surf_M = plot_surf_M*do_plot;
    %
    % (13)
    % only for multiple B value analysis
    % To plot normalized M as function of time and external DC magnetic 
    % field B in an extracted time to be chosen by Ts,Tf,
    % set plot_surf_norm_extr_M = 1;
    plot_surf_norm_extr_M = 0;
    plot_surf_norm_extr_M = plot_surf_norm_extr_M*do_plot;
    % Define start and end time of extraction
    % Be aware of the units of time of T_Global!
    Ts = 4*10^-7;
    Tf = 8*10^-7;
    %
    Ts = Ts*t_correction;
    Tf = Tf*t_correction;
    %
    % (14)
    % only for multiple B value analysis
    % To plot the amplitude of S_initial+S_sp as function of both time and
    % external B, set plot_surf_S_tot = 1
    plot_surf_S_tot = 0;
    plot_surf_S_tot = plot_surf_S_tot*do_plot;
    %
    % (15)
    % only for multiple B value analysis
    % To plot the amplitude of S_sp as function of both time and
    % external B, set plot_surf_S_sp_tot = 1
    plot_surf_S_sp_tot = 0;
    plot_surf_S_sp_tot = plot_surf_S_sp_tot*do_plot;
    %
    % For 3D plots, choose if to display all pulses in legend or just the
    % first and last
    all_pls_x_lin = 0;
    %
    %% plot traces   
    %
    if disp_only_last_plots
        close all
    end
    %
    if param_ind ==1
        close all
    % define a figure index to be run on when poltting distinct figured
        fig_indx = 1;
        fig_indx_initial = 1;
    end
    %
    % The offset is used to distinguish between plots of different B
    % it simply puts the plot one above the other
    offset_L = 0;
    %
    M_mat=zeros(length(H_Xtrnl_vect),length(T_Global));
    %
    if plot_param
        % Empty plot to show annotations on
        figure(fig_indx)
        plot(T_Global,0*(M_array{L}(1,:)),'linewidth',0.5,'color','w'); 
        title('Simulation Parameters','FontSize',24)
        %   
        % Create two textboxes containing the parameters with which the simulation
        % ran - A_feedback_xy,A_feedback_z,alpha_SP,H_Xtrnl,Jc_dc,Jc_ac,N_Dmag,H_amp_CW,H_Pulse_amp
        % The annotation appear in a single empty plot.
        %
        dim = [.7 .33 .18 .59];
        str = {'- SP FEEDBACKS -','SP closed loop ',num2str(SP_closed_loop),'g_S_P conductivity [1/m^2] ',num2str(g_SP),'{\alpha}_S_P ',num2str(alpha_SP),'{\alpha}_S_P / {\alpha}_G_i_l_b_e_r_t ',num2str(alpha_SP/alpha),'SP efficiency ratio ',num2str(SP_efficiency_ratio),'dM/dt xy SP ', num2str(Im_SP_xy),'dM/dt z SP ', num2str(Im_SP_z),'dM/dt feedback ',num2str(A_feedback)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        dim = [.45 .33 .18 .59];
        str = {'{Js}_d_c / {Js}_S_T_T ',num2str(STT_DC_ratio),'{Jc}_d_c [A/m^2] ',['10^',num2str(log10(Jc_dc))],'{Hshe}_d_c ',num2str(Hshe_dc),'{Jc}_a_c [A/m^2] ',['10^',num2str(log10(Jc_ac))],'{Hshe}_a_c ',num2str(Hshe_ac),'Start time {Jc}_a_c [nsec] ',num2str((t_correction)*T_start_S_ac),'End time of {Jc}_a_c [nsec] ',num2str((t_correction)*T_end_S_ac),'{Jc}_a_c starting phase (deg) ',num2str(S_ac_Phase_rf),'Demag parameter ',num2str(N0)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        dim = [.2 .33 .18 .59];
        str = {'- GENERAL -','{\alpha}_G_i_l_b_e_r_t ',num2str(alpha),'RF Frequency [GHz] ',num2str(Freq_CW*(f_correction)),'B_r_e_s_o_n_a_n_c_e [T] ',num2str(w/gma),'Spacing Between B [T] ',num2str(B_width*half_delta_B/B_samples),'# B values ',num2str(B_samples),'N_D_m_a_g ', num2str(N_Dmag'),'Anisotropy Field Constant ',num2str(K0),'M_z starting value ',num2str(M_tg(3,1)),'RF field amplitude',num2str(H_amp_CW)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        fig_indx = fig_indx + 1;
        %
        figure(fig_indx)
        plot(T_Global,0*(M_array{L}(1,:)),'linewidth',0.5,'color','w'); 
        title('Simulation Parameters 2','FontSize',24)
        %   
        % Create two textboxes containing the parameters with which the simulation
        % ran - A_feedback_xy,A_feedback_z,alpha_SP,H_Xtrnl,Jc_dc,Jc_ac,N_Dmag,H_amp_CW,H_Pulse_amp
        % The annotation appear in a single empty plot.
        %
        dim = [.7 .33 .18 .59];
        str = {'- PUMP FIELD 2 -','B_p_u_m_p_ sgn change ',num2str(H_pump_change_sign),'|Mz| val for sgn change ',num2str(Mz_val_for_H_pump_sgn_change),'B_p_u_m_p turn off ',num2str(H_pump_turn_off),'|Mz| val for turn off ',num2str(Mz_val_for_H_pump_turn_off),'set pulse start time ','as switching point ',num2str(set_t_pls_as_switching_point),'IFE coefficient ',num2str(IFE_field_amplitude)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        dim = [.45 .33 .18 .59];
        str = {'- PUMP FIELD -','B_p_u_m_p_ amplitude ',num2str(H_pump_0),'B_p_u_m_p_ start time [nsec] ',num2str((t_correction)*H_pump_t_start),'B_p_u_m_p_ end time [nsec] ',num2str((t_correction)*H_pump_t_end),'B_p_u_m_p_ on/off mode ',num2str(H_pump_on_off),'B_p_u_m_p_ on/off time period ',num2str(H_pump_on_off_time_period),'B_p_u_m_p_ chop ',num2str(H_pump_chop),'B_p_u_m_p_ chop value ',num2str(H_pump_chop_val)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        dim = [.2 .33 .18 .59];
        str = {'- PULSES -','B-Pulse-amp ', num2str(B_Pulse_amp),'set pulse start time ','as switching point ',num2str(set_t_pls_as_switching_point),'Pulse start time [nsec] ',num2str((t_correction)*t0_H_pulsed),'Pulse end time [nsec] ',num2str((t_correction)*t0_H_pulsed_end), '# pulses ',num2str(num_pulses),'Time interval btw pls [nsec] ',num2str((t_correction)*time_int_betw_pls),'B_r_e_s / Pulse amp ',num2str(res_field_over_pls_amp)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20)
        %
        fig_indx = fig_indx + 1;
        %
    end
    %
   T_end_S_ac = T_end_S_ac*t_correction;
   t0_H_pulsed = t0_H_pulsed*t_correction;
   time_int_betw_pls = time_int_betw_pls*t_correction;
    %
    %% plot time trafigureces of magnetization
    
%    figure;
    colors = colormap(lines(length(H_Xtrnl_vect)));
    if plot_M_env_osc
        for L=1:length(H_Xtrnl_vect)
            for j = 1:3
                Mag_compnnt=j;    %chooses which component to analyze x,y,or z  (1,2,or 3)
                for l=1:length(H_Xtrnl_vect)
                    max_L(l)=max(abs(M_array{l}(Mag_compnnt,:)));   % This extracts a global maximum of the specific magnetization component
                end
                M_max=max(max_L);
                figure(fig_indx)
                plot(T_Global,offset_L*2*(L-1)+M_array{L}(Mag_compnnt,:)/M_max,'linewidth',1.5,'color',colors(L,:)); 
                %
                title({['Normalized Magnetization component ',num2str(Mag_compnnt),' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
                %
                if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
                end
                if disp_time_pico
                    xlabel('Time (psec)','fontsize',axis_font)
                end
                if disp_time_nano
                    xlabel('Time (nsec)','fontsize',axis_font)
                end
                %
                ylabel(['Norm. M ',num2str(Mag_compnnt),' component (a.u)'],'fontsize',axis_font)
                grid
                xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
                fig_indx = fig_indx +1;
            end
        end 
    end
    %

    %% plot time traces of spin current
    if plot_S_ac_env_osc
        for L=1:length(H_Xtrnl_vect)
            for i = 1:3
                S_ac_compnnt = i;
                for l=1:length(H_Xtrnl_vect)
                    max_L(l)=max(abs(S_ac_array{l}(S_ac_compnnt,:)));   % This extracts a global maximum of the specific magnetization comonent
                end
                M_max_S_ac=max(max_L);
                figure(fig_indx)
                plot(T_Global,offset_L*2*(L-1)+S_ac_array{L}(S_ac_compnnt,:)/M_max_S_ac,'linewidth',1.5,'color',colors(L,:)); 
                title({['Normalized J_S ac component ',num2str(S_ac_compnnt),' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
                %
                if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
                end
                if disp_time_pico
                    xlabel('Time (psec)','fontsize',axis_font)
                end
                if disp_time_nano
                    xlabel('Time (nsec)','fontsize',axis_font)
                end
                %
                ylabel(['Norm. J_S ',num2str(S_ac_compnnt),'component','(a.u)'],'fontsize',axis_font)
                grid 
                hold on
                xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
                fig_indx = fig_indx + 1;
            end
        end
    end
    %
    if plot_S_ac_legend_not_norm
        for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            if view_initial_and_final_S_ac
                % For S_ac initially polarized in x direction, set S_ac_initial(1,:)
                plot(T_Global,offset_L*2*(L-1)+ S_ac_initial_array{L}(1,:),'linewidth',1.5,'DisplayName','J_S_{-in} X component');
                hold on 
                plot(T_Global,(offset_L*2*(L-1)+ S_ac_array{L}(1,:)),'linewidth',1.5,'DisplayName','J_S X component');
                hold on
            end
            plot(T_Global,offset_L*2*(L-1)+ S_ac_array{L}(1,:)-S_ac_initial_array{L}(1,:),'linewidth',1.5,'DisplayName','J_S X component initial - final');
            hold on 
            plot(T_Global,(offset_L*2*(L-1)+ S_ac_array{L}(2,:)),'linewidth',1.5,'DisplayName','J_S Y component');
            hold on
            plot(T_Global,offset_L*2*(L-1)+ S_ac_array{L}(3,:),'linewidth',1.5,'DisplayName','J_S Z component');
            hold off
            title({['J_S components ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
            %
            if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
            end            
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('J_S (a.u)','fontsize',axis_font)
            grid           
            lgd = legend;
            lgd.FontSize = lgd_font;
            lgd.Title.String = 'J_S components';
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
           fig_indx = fig_indx +1;
        end 
    end
    %
    if plot_S_sp_legend_not_norm
        for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            plot(T_Global,offset_L*2*(L-1)+ S_sp_array{L}(1,:),'linewidth',1.5,'DisplayName','J_S_P X component');
            hold on 
            plot(T_Global,offset_L*2*(L-1)+ S_sp_array{L}(2,:),'linewidth',1.5,'DisplayName','J_S_P Y component');
            hold on
            plot(T_Global,offset_L*2*(L-1)+ S_sp_array{L}(3,:),'linewidth',1.5,'DisplayName','J_S_P Z component');
            hold off
            title({['J_S_P components ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('J_S_P (a.u)','fontsize',axis_font)
            grid            
            lgd = legend;
            lgd.FontSize = lgd_font;
            lgd.Title.String = 'J_S_P components';
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
           ax = gca;
           ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
           ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
           fig_indx = fig_indx +1;
        end 
    end
    %
    if plot_total_S_sp_initial
        for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            plot(T_Global,offset_L*2*(L-1) + Hshe_ac*S_ac_initial_array{L}(1,:) + sign_J_s_J_sp*S_sp_array{L}(1,:),'linewidth',1.5,'DisplayName','J_S_P + J_S_{-in} X component');
            yline(max(Hshe_ac*S_ac_initial_array{L}(1,:)),'-g',{'Upper envelope','of J_S_{-in}'},'linewidth',2,'DisplayName','Upper envelope of J_S_{-in}','fontsize',lgd_font);
            yline(min(Hshe_ac*S_ac_initial_array{L}(1,:)),'-y',{'Lower envelope','of J_S_{-in}'},'linewidth',2,'DisplayName','Lower envelope of J_S_{-in}','fontsize',lgd_font);
            hold on 
            if  ~tot_sp_int_only_x 
                plot(T_Global,offset_L*2*(L-1) + Hshe_ac*S_ac_initial_array{L}(2,:) + sign_J_s_J_sp*S_sp_array{L}(2,:),'linewidth',1.5,'DisplayName','J_S_P + J_S_{-in} Y component');
                hold on
                plot(T_Global,offset_L*2*(L-1) + Hshe_ac*S_ac_initial_array{L}(3,:) + sign_J_s_J_sp*S_sp_array{L}(3,:),'linewidth',1.5,'DisplayName','J_S_P + J_S_{-in} Z component');
                hold off
            end
            title({['J_S_P + J_S_{-in} components ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('J_S_P + J_S_{-in} (a.u)','fontsize',axis_font)
            grid         
            lgd = legend;
            lgd.FontSize = lgd_font;
            lgd.Title.String = 'J_S_P + J_S_{-in} components';
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            %
            if plot_phase_entire_time
                % Calculate the phase between J_ac and J_sp for the entire time 
                dot_product = dot(Hshe_ac*S_ac_initial_array{L}(1,:),sign_J_s_J_sp*S_sp_array{L}(1,:));
                norm_product = (norm(Hshe_ac*S_ac_initial_array{L}(1,:)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,:)));
                phase_shift_in_radians(L) = acos(dot_product/norm_product);
                fig_indx = fig_indx +1;
            end
            %
            if ~plot_phase_entire_time
                t0_H_pulsed = t0_H_pulsed*(f_correction);
                % Calculate phase between J_ac and J_sp before pulses, driven by J_ac only
                dot_product_bfr_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11),sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11));
                norm_product_bfr_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11))));
                phase_shift_in_radians_bfr_pls(L) = acos(dot_product_bfr_pls/norm_product_bfr_pls);
                % Calculate phase between J_ac and J_sp during pulses
                dot_product_aft_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt) - T_int2):round(t0_H_pulsed_end/dt)-T_int22),sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt) - T_int2):round(t0_H_pulsed_end/dt)-T_int22));
                norm_product_aft_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt) - T_int2):round(t0_H_pulsed_end/dt)-T_int22)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt) - T_int2):round(t0_H_pulsed_end/dt)-T_int22)));
                phase_shift_in_radians_aft_pls(L) = acos(dot_product_aft_pls/norm_product_aft_pls);
                t0_H_pulsed = t0_H_pulsed*(t_correction);
                fig_indx = fig_indx +1;
            end
            %
        end
            %
        % Plot the phase as a function of the external field
        if plot_phase_entire_time
            figure(fig_indx)
            plot(B_Xtrnl_vect,phase_shift_in_radians*(360/(2*pi)),'linewidth',4,'DisplayName','Phase response J_S_P to J_S_{-in}')
            hold on
            scatter(B_Xtrnl_vect,phase_shift_in_radians*(360/(2*pi)),'filled','LineWidth',3)
            title('{\phi}{_{J_S_{-in}}_{J_S_P}} ','FontSize',title_font + smaller_font_diff)
            xlabel('B (T)','fontsize',axis_font)
            ylabel('Phase (deg)','fontsize',axis_font)
            grid
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            fig_indx = fig_indx +1;
        end
        %
        if ~plot_phase_entire_time
            figure(fig_indx)
            plot(B_Xtrnl_vect,phase_shift_in_radians_bfr_pls*(360/(2*pi)),'linewidth',4,'DisplayName','Phase response J_S_P to J_S_{-in} before pulses')
            hold on
            scatter(B_Xtrnl_vect,phase_shift_in_radians_bfr_pls*(360/(2*pi)),'filled','LineWidth',3)
            title('{\phi}{_{J_S_{-in}}_{J_S_P}} before pulses ','FontSize',title_font + smaller_font_diff)
            xlabel('B (T)','fontsize',axis_font)
            ylabel('Phase (deg)','fontsize',axis_font)
            grid
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            fig_indx = fig_indx +1;
            %
            figure(fig_indx)
            plot(B_Xtrnl_vect,phase_shift_in_radians_aft_pls*(360/(2*pi)),'linewidth',4,'DisplayName','Phase response J_S_P to J_S_{-in} during pulses')
            hold on
            scatter(B_Xtrnl_vect,phase_shift_in_radians_aft_pls*(360/(2*pi)),'filled','LineWidth',3)
            title('{\phi}{_{J_S_{-in}}_{J_S_P}} during pulses ','FontSize',title_font + smaller_font_diff)
            xlabel('B (T)','fontsize',axis_font)
            ylabel('Phase (deg)','fontsize',axis_font)
            grid
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            fig_indx = fig_indx +1;
        end
        %
    end    
    %
    %% Plot the primary and damping torque as a function of time
    if plot_toruqes
            for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            subplot(2,1,1);
            plot(T_Global,primary_torque_array{L}(1,:),'linewidth',2,'color','#77AC30','DisplayName','primary torque x');
            hold on
            plot(T_Global,primary_torque_array{L}(2,:),'linewidth',2,'color','#4DBEEE','DisplayName','primary torque y');
            hold on
            plot(T_Global,primary_torque_array{L}(3,:),'linewidth',2,'color','m','DisplayName','primary torque z');
            hold on
            title({['Primary Torque ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',round(0.5*title_font))
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('dM/dt [A/(m*sec)]','fontsize',axis_font)
            lgd = legend;
            lgd.FontSize =lgd_font;
            lgd.Title.String = 'primary torque components';
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = round(0.25*ax_num_font_size); %for x-axis 
            ax.YAxis.FontSize = round(0.25*ax_num_font_size); %for y-axis 
            fig_indx = fig_indx +1;
            %
            subplot(2,1,2);
            plot(T_Global,damping_torque_array{L}(1,:),'linewidth',2,'color','#77AC30','DisplayName','damping torque x');
            hold on
            plot(T_Global,damping_torque_array{L}(2,:),'linewidth',2,'color','#4DBEEE','DisplayName','damping torque y');
            hold on
            plot(T_Global,damping_torque_array{L}(3,:),'linewidth',2,'color','#A2142F','DisplayName','damping torque z');
            hold on
            title({['Damping Torque ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',round(0.5*title_font))
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('dM/dt [A/(m*sec)]','fontsize',axis_font)
            lgd = legend;
            lgd.FontSize =lgd_font;
            lgd.Title.String = 'damping torque components';
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = round(0.25*ax_num_font_size); %for x-axis 
            ax.YAxis.FontSize = round(0.25*ax_num_font_size); %for y-axis 
            fig_indx = fig_indx +1;
            % #D95319 #A2142F
            figure(fig_indx)
            subplot(3,1,1);
            plot(T_Global(1:end-1),damping_torque_array{L}(3,1:end-1)+primary_torque_array{L}(3,1:end-1),'linewidth',2,'color','#D95319','DisplayName','total torque z');
            hold on
            plot(T_Global(1:end-1),diff(M_array_Rot_z{L}(1,1:end))/dt_low_res,'linewidth',2,'color','b','DisplayName','primary torque z');
            title({['Total z Torque ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',round(0.5*title_font))
%             ylim([-5*10^19 5*10^19]);
%             xlim([60 75]);
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('dM/dt [A/(m*sec)]','fontsize',axis_font)
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = round(0.25*ax_num_font_size); %for x-axis 
            ax.YAxis.FontSize = round(0.25*ax_num_font_size); %for y-axis 
            subplot(3,1,2);
            plot(T_Global,((damping_torque_array{L}(3,:)+primary_torque_array{L}(3,:))./((( (damping_torque_array{L}(1,:)+primary_torque_array{L}(1,:)).^2 + (damping_torque_array{L}(2,:)+primary_torque_array{L}(2,:)).^2 + (damping_torque_array{L}(3,:)+primary_torque_array{L}(3,:)).^2  )).^0.5)),'linewidth',2,'color','k','DisplayName','total torque z fraction from total torque');
            title({['Total z Torque Fraction ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',round(0.5*title_font))
%             ylim([-5*10^19 5*10^19]);
%             xlim([60 75]);
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('dM/dt [A/(m*sec)]','fontsize',axis_font)
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = round(0.25*ax_num_font_size); %for x-axis 
            ax.YAxis.FontSize = round(0.25*ax_num_font_size); %for y-axis 
            subplot(3,1,3);
            plot(T_Global,damping_torque_array{L}(1,:)+primary_torque_array{L}(1,:),'linewidth',2,'color','#77AC30','DisplayName','total torque x');
            hold on
            plot(T_Global,damping_torque_array{L}(2,:)+primary_torque_array{L}(2,:),'linewidth',2,'color','#4DBEEE','DisplayName','total torque y');
            title({['Total x and y Torque ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',round(0.5*title_font))
%             ylim([-5*10^19 5*10^19]);
%             xlim([60 75]);
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('dM/dt [A/(m*sec)]','fontsize',axis_font)
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = round(0.25*ax_num_font_size); %for x-axis 
            ax.YAxis.FontSize = round(0.25*ax_num_font_size); %for y-axis 
            fig_indx = fig_indx +1;
        end  
    end
    %
     %% plot traces of Rotated Magnetization

    for L=1:length(H_Xtrnl_vect)
        % This extracts a global maximum of all magnetization components
        max_L_Rot_x(L)=max(abs(M_array_Rot_x{L}(1,:)));
        max_L_Rot_y(L)=max(abs(M_array_Rot_y{L}(1,:)));
        max_L_Rot_z(L)=max(abs(M_array_Rot_z{L}(1,:)));
    end
    M_max_Rot_x=max(max_L_Rot_x);
    M_max_Rot_y=max(max_L_Rot_y);
    M_max_Rot_z=max(max_L_Rot_z);
    %
    M_mat_Rot=zeros(length(H_Xtrnl_vect),length(T_Global));
      
    %% plot time traces of normalized Rotated Magnetization    
        
%    figure;
    colors=colormap(lines(length(H_Xtrnl_vect)));
    if plot_M_env
        % At  overdamping, one transverse component might look as if it is
        % oscillating ('unrotated') - the oscillations are really weak,
        % which can be seen in the plot_M_legend_not_norm figure.
        for L=1:length(H_Xtrnl_vect)
                % plot Mx
                figure(fig_indx)
                subplot(3,1,1)
                plot(T_Global,abs(offset_L*2*(L-1)+M_array_Rot_x{L}(1,:)/M_max_Rot_x)); 
                title({['Normalized Rotated M_x Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'fontsize',title_font)
                %
                if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
                end
                if disp_time_pico
                    xlabel('Time (psec)','fontsize',axis_font)
                end
                if disp_time_nano
                    xlabel('Time (nsec)','fontsize',axis_font)
                end
                %
                ylabel('Norm. Rotated Mx (a.u)','fontsize',axis_font)
                xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
                hold on
                % plot My
                subplot(3,1,2)
                plot(T_Global,abs(offset_L*2*(L-1)+M_array_Rot_y{L}(1,:)/M_max_Rot_y)); 
                title({['Normalized Rotated M_y Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'fontsize',title_font)
                %
                if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
                end
                if disp_time_pico
                    xlabel('Time (psec)','fontsize',axis_font)
                end
                if disp_time_nano
                    xlabel('Time (nsec)','fontsize',axis_font)
                end
                %
                ylabel('Norm. Rotated My (a.u)','fontsize',axis_font)
                xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
                hold on
                % plot Mz
                subplot(3,1,3)
                plot(T_Global,offset_L*2*(L-1)+M_array_Rot_z{L}(1,:)/M_max_Rot_z); 
                title({['Normalized Rotated M_z Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'fontsize',title_font)
                %
                if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
                end
                if disp_time_pico
                    xlabel('Time (psec)','fontsize',axis_font)
                end
                if disp_time_nano
                    xlabel('Time (nsec)','fontsize',axis_font)
                end
                %
                ylabel('Norm. Rotated Mz (a.u)','fontsize',axis_font)
                hold on
                xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
                fig_indx = fig_indx + 1;
                %
        end
    end
    %% plot time traces of NOT normalized Rotated Magnetization    
        
    if plot_M_legend_not_norm
        for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            plot(T_Global,(offset_L*2*(L-1)+ M_array_global{L}(1,:)),'linewidth',2,'DisplayName','M_x');
            hold on
            plot(T_Global,(offset_L*2*(L-1)+ M_array_global{L}(2,:)),'linewidth',2,'DisplayName','M_y');
            hold on
            plot(T_Global,offset_L*2*(L-1)+ M_array_global{L}(3,:),'linewidth',4,'DisplayName','M_z');
            hold on
            plot(T_Global,offset_L*2*(L-1)+((M_array_global{L}(1,:)).^2 + (M_array_global{L}(2,:)).^2 + (M_array_global{L}(3,:)).^2).^0.5,'linewidth',1.5,'DisplayName','Magnetization Absolute Value');
            hold off
            title({['Magnetization components ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
            %
            if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('M [A/m]','fontsize',axis_font)
%             grid
%             lgd = legend;
%             lgd.FontSize =lgd_font;
%             lgd.Title.String = 'Magnetization components';
%             xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            fig_indx = fig_indx +1;
        end 
    end
    %
    if plot_M_legend_not_norm_rot
        for L=1:length(H_Xtrnl_vect)
            figure(fig_indx)
            plot(T_Global,abs(offset_L*2*(L-1)+ M_array_Rot_x{L}(1,:)),'linewidth',2,'DisplayName','Rotated M_x');
            hold on
            plot(T_Global,abs(offset_L*2*(L-1)+ M_array_Rot_y{L}(1,:)),'linewidth',2,'DisplayName','Rotated M_y');
            hold on
            plot(T_Global,offset_L*2*(L-1)+ M_array_Rot_z{L}(1,:),'linewidth',3,'DisplayName','Rotated M_z');
            hold on
            plot(T_Global,offset_L*2*(L-1)+((M_array_Rot_z{L}(1,:)).^2 + (M_array_Rot_y{L}(1,:)).^2 + (M_array_Rot_x{L}(1,:)).^2).^0.5,'linewidth',1.5,'DisplayName','Magnetization Absolute Value');
            hold off
            title({['Rotated Magnetization components ',' Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font)
            %
            if disp_time_femto
                    xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font)
            end
            %
            ylabel('Rotated M [A/m]','fontsize',axis_font)
            grid
            lgd = legend;
            lgd.FontSize =lgd_font;
            lgd.Title.String = 'Rotated Magnetization components';
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
            fig_indx = fig_indx +1;
        end 
    end
    %
%  %   NOT FOR FIGURE
%     if plot_Mz_legend_offset
%         figure(fig_indx)
%         for L=1:length(H_Xtrnl_vect)
%             plot(T_Global,offset_L*2*(L-1)+ M_array_Rot_z{L}(1,:),'linewidth',1.5,'DisplayName',['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']);
%             hold on
%         end 
%         title('NOT Normalized Mz for all B Time Traces','fontsize',title_font)
                %
%                 if disp_time_femto
%                     xlabel('Time (fsec)','fontsize',axis_font)
%                 end
%             if disp_time_pico
%                 xlabel('Time (psec)','fontsize',axis_font)
%             end
%             if disp_time_nano
%                 xlabel('Time (nsec)','fontsize',axis_font)
%             end
%         ylabel('NOT Norm & offset. Mz (a.u)','fontsize',axis_font)
%         grid
% %       h = text(450,650,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
% %       set(h,'color','g','fontsize',24)
%         lgd = legend;
%         lgd.FontSize = lgd_font;
%         lgd.Title.String = 'NOT Normalized Mz for B values';
%         xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation','fontsize',lgd_font);
%         for z = 1:num_pulses
%             xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%         end
%         ax = gca;
%         ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
%         ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
%         fig_indx = fig_indx +1;
%     end
    % FOR FIGURE
    if plot_Mz_legend_offset
        figure(fig_indx)
        for L=1:length(H_Xtrnl_vect)
            plot(T_Global,offset_L*2*(L-1)+ M_array_Rot_z{L}(1,:),'linewidth',3.5,'DisplayName',['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']);
            hold on
        end 
        %title('M_z(t) on resonance','fontsize',title_font)
        title('M_z(t)','fontsize',title_font)
        %
        if disp_time_femto
            xlabel('Time (fsec)','fontsize',axis_font)
        end
        if disp_time_pico
            xlabel('Time (psec)','fontsize',axis_font)
        end
        if disp_time_nano
            xlabel('Time (nsec)','fontsize',axis_font)
        end
        %
        ylabel('Mz ([A/m])','fontsize',axis_font)
        lgd = legend;
        lgd.FontSize = lgd_font;
        lgd.Title.String = 'M_z(t) for B values';
        xline(T_end_S_ac,'-k',{'End time','of excitation'},'linewidth',6,'DisplayName','End time of excitation','fontsize',45);
        % FOR INSET
%          for z = 1:num_pulses
%              xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%          end
        % NOT FOR INSET
        xline(t0_H_pulsed ,'-k',{'pulse start'},'linewidth',6,'DisplayName','pulse start','FontSize',title_font,'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left'); 
        xline(t0_H_pulsed + (num_pulses-1)*time_int_betw_pls,'-k',{'pulse end'},'linewidth',6,'DisplayName','pulse end','FontSize',title_font,'LabelVerticalAlignment','bottom');
        %
        ax = gca;
        ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
        ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
        %
        fig_indx = fig_indx +1;
    end
    %
    if plot_M_with_Power_Spectrum
        for L=1:length(H_Xtrnl_vect)
            %
            % Mark the moment an optical excitation distrubes the steady precession. used in the FMR regime, not in the optical regime.
            disp_H_pulse = 0;
            %
            f_correction_optical = 10^-15; % Defined for the Fourier transform visibility.
            % Mx
            figure(fig_indx)
            subplot(2,1,1)
            plot(T_Global,(M_array{L}(1,:)/(max(abs(M_array{L}(1,:))))),'b','linewidth',0.5); 
            title({['Normalized M_x Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %            
            ylabel('Norm. Mx (a.u)','fontsize',axis_font - smaller_font_diff)
            set(gcf,'outerposition',[10 65 1370 500])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            if disp_H_pulse
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            % 
            % F{Mx}
            M_zero_pad=[M_array{L}(1,:)/max(abs(M_array{L}(1,:))) zero_pad_array_for_fourier] ;
            [S_1{L},f_shift_1{L}]=trans_fourier(M_zero_pad,dt);
            subplot(2,1,2)
            plot((f_correction_optical)*f_shift_1{L},abs(S_1{L}.^2)/max(abs(S_1{L}.^2)),'linewidth',2.5)
            title({['Normalized M_x Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{Mx}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_1{L});
            freq(L)=abs(f_shift_1{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            % My
            figure(fig_indx)
            subplot(2,1,1)
            plot(T_Global,(M_array{L}(2,:)/(max(abs(M_array{L}(2,:))))),'b','linewidth',0.5); 
            title({['Normalized M_y Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %           
            ylabel('Norm. My (a.u)','fontsize',axis_font - smaller_font_diff)
            set(gcf,'outerposition',[10 65 1370 500])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            if disp_H_pulse
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            % 
            % F{My}
            M_zero_pad=[M_array{L}(2,:)/max(abs(M_array{L}(2,:))) zero_pad_array_for_fourier] ;
            [S_2{L},f_shift_2{L}]=trans_fourier(M_zero_pad,dt);
            subplot(2,1,2)
            plot((f_correction_optical)*f_shift_2{L},abs(S_2{L}.^2)/max(abs(S_2{L}.^2)),'linewidth',2.5)
            title({['Normalized M_y Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{My}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_2{L});
            freq(L)=abs(f_shift_2{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            % Mz
            figure(fig_indx)
            subplot(2,1,1)
            plot(T_Global,(M_array{L}(3,:)/(max(abs(M_array{L}(3,:))))),'b','linewidth',1); 
            title({['Normalized M_z Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %           
            ylabel('Norm. Mz (a.u)','fontsize',axis_font - smaller_font_diff)
            set(gcf,'outerposition',[10 65 1370 500])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            if disp_H_pulse
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            %
            % F{Mz}
            M_zero_pad=[M_array{L}(3,:)/max(abs(M_array{L}(3,:))) zero_pad_array_for_fourier] ;
            [S_3{L},f_shift_3{L}]=trans_fourier(M_zero_pad,dt);
            subplot(2,1,2)
            plot((f_correction_optical)*f_shift_3{L},abs(S_3{L}.^2)/max(abs(S_3{L}.^2)),'linewidth',2.5)
            title({['Normalized M_z Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{Mz}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_3{L});
            freq(L)=abs(f_shift_3{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            figure(fig_indx)
            subplot(2,1,1)
            plot((f_correction_optical)*f_shift_1{L},abs(S_1{L}.^2)/max(abs(S_1{L}.^2)),'linewidth',2.5)
            title({['Normalized Mx zoomed Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{Mx}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([(f_correction)*(Freq_CW-Freq_CW_fraction_for_fourier_plot) (f_correction)*(Freq_CW+Freq_CW_fraction_for_fourier_plot)])
            xline((f_correction)*Freq_CW,'-m',{'RF frequency','of ac excitation'},'linewidth',3,'DisplayName','frequency of RF excitation','fontsize',lgd_font); 
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            subplot(2,1,2)
            plot((f_correction)*f_shift_2{L},abs(S_2{L}.^2)/max(abs(S_2{L}.^2)),'linewidth',2.5)
            title({['Normalized My zommed Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{My}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([(f_correction)*(Freq_CW-Freq_CW_fraction_for_fourier_plot) (f_correction)*(Freq_CW+Freq_CW_fraction_for_fourier_plot)])
            xline((f_correction)*Freq_CW,'-m',{'RF frequency','of ac excitation'},'linewidth',3,'DisplayName','frequency of RF excitation','fontsize',lgd_font);
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            %
            fig_indx = fig_indx +1;
        end
        if plot_Mx_My_filtered
            figure(fig_indx)
            y_limit_in_T_Hz = 2000;
            % Define filter specifications
            fs = 1/dt_low_res; % sample frequency is the inverse of the temporal resolution
            cutoffFrequency = 100e12;  % Cutoff frequency (100 THz, 100000 THz in Hz)
            filterOrder = 500;   % Filter order (adjust as needed)

            % Normalize the cutoff frequency
            normalizedCutoff = cutoffFrequency / (fs/2);

            % Design the low-pass filter using fir1
            filterCoefficients = fir1(filterOrder, normalizedCutoff, 'low');

            % Frequency response of the filter
            freqResponse = fftshift(abs(fft(filterCoefficients, 1000)));

            % Frequency vector for plotting
            frequencies = linspace(-fs/2, fs/2, 1000);
            % 
            % Plot the frequency response
            plot(frequencies, freqResponse);
            xlabel('Frequency (Hz)');
            ylabel('Magnitude');
            title('Low-Pass Filter Frequency Response');
            xlim([-500 500]*10^12)
            % 
            fig_indx = fig_indx +1;
            figure(fig_indx)
            filteredMy = filter(filterCoefficients, 1, M_array_global{L}(2,:));
            subplot(2,1,1)
            %filteredMx = filter(filterOrder, normalizedCutoff, M_array_global{L}(1,:));
            plot(T_Global,(filteredMy(1,:)/(max(abs(filteredMy(1,:))))),'b','linewidth',0.5); 
            title({['Normalized filtered M_y Time Traces'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %            
            ylabel('Norm. filtered Mx (a.u)','fontsize',axis_font - smaller_font_diff)
            set(gcf,'outerposition',[10 65 1370 500])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            if disp_H_pulse
                for z = 1:num_pulses
                    xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
                end
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            subplot(2,1,2)
            M_zero_pad=[filteredMy(1,:)/max(abs(filteredMy(1,:))) zero_pad_array_for_fourier] ;
            [S_1{L},f_shift_1{L}]=trans_fourier(M_zero_pad,dt);
            plot((f_correction_optical)*f_shift_1{L},abs(S_1{L}.^2)/max(abs(S_1{L}.^2)),'linewidth',2.5)
            title({['Normalized filtered M_y Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{Mx}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_1{L});
            freq(L)=abs(f_shift_1{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            figure(fig_indx)
            instfreq(filteredMy/M0,fs)
            ylim([0 y_limit_in_T_Hz]);
            title('Instantaneous Frequency of filtered Mx and My');
            fig_indx = fig_indx +1;
            %
        end
    end
    %    
    if plot_H_pump            
            % H_pump_x
            figure(fig_indx)
            subplot(2,1,1)
            plot(T_Global,(H_pump_array{L}(1,:)),'b','linewidth',1); 
            title({['H_p_u_m_p x component Time Trace'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %           
            ylabel('H_p_u_m_p x (A/m)','fontsize',axis_font - smaller_font_diff)
%            set(gcf,'outerposition',[10 65 1370 500])
%            axis([0 5e+04 -1.5 1.5])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            %text(1.3e+04,0.6,'H_0 = ',num2str(H_Xtrnl_vect(L)),' T','FontSize',35)
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            %
            % F{Bx}
            B_zero_pad=[(H_pump_array{L}(1,:)) zero_pad_array_for_fourier] ;
            [S_1{L},f_shift_1{L}]=trans_fourier(B_zero_pad,dt);
            subplot(2,1,2)
            plot((f_correction)*f_shift_1{L},abs(S_1{L}.^2)/max(abs(S_1{L}.^2)))
            title({['H_p_u_m_p x normalized Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{B_p_u_m_p x}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_1{L});
            freq(L)=abs(f_shift_1{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            % H_pump_y
            figure(fig_indx)
            subplot(2,1,1)
            plot(T_Global,(H_pump_array{L}(2,:)),'b','linewidth',1); 
            title({['H_p_u_m_p y component Time Trace'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %           
            ylabel('H_p_u_m_p y (A/m)','fontsize',axis_font - smaller_font_diff)
%             set(gca,'fontsize',20)
%             set(gca,'linewidth',2)
%            set(gcf,'outerposition',[10 65 1370 500])
%            axis([0 5e+04 -1.5 1.5])
            xline(T_end_S_ac,'-m',{'End time','of excitation'},'linewidth',1.5,'DisplayName','End time of excitation','fontsize',lgd_font);
            %text(1.3e+04,0.6,'H_0 = ',num2str(H_Xtrnl_vect(L)),' T','FontSize',35)
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
            end
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            %
            % F{By}
            B_zero_pad=[(H_pump_array{L}(2,:)) zero_pad_array_for_fourier] ;
            [S_2{L},f_shift_2{L}]=trans_fourier(B_zero_pad,dt);
            subplot(2,1,2)
            plot((f_correction)*f_shift_2{L},abs(S_2{L}.^2)/max(abs(S_2{L}.^2)))
            title({['H_p_u_m_p y normalized Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            if disp_time_femto
                xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_pico
                xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
            end
            ylabel('Norm. |F{B_p_u_m_p y}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
    %             h = text(900,100,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
    %             set(h,'color',colors(L,:),'fontsize',18)
            xlim([-F_range_for_fourier_plots F_range_for_fourier_plots])
            [val,indx]=max(S_2{L});
            freq(L)=abs(f_shift_2{L}(indx));
            %
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            fig_indx = fig_indx +1;
            %
            if plot_zoomed_H_pump_spectrum
                % plot F{Bx}, F{By} near the RF spectrum
                figure(fig_indx)
                subplot(2,1,1)
                plot((f_correction)*f_shift_1{L},abs(S_1{L}.^2)/max(abs(S_1{L}.^2)))
                title({['H_p_u_m_p x normalized zoomed Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
                if disp_time_femto
                    xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_pico
                    xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_nano
                    xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
                end
                ylabel('Norm. |F{H_p_u_m_p x}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
                xlim([(f_correction)*(Freq_CW-Freq_CW_fraction_for_fourier_plot) (f_correction)*(Freq_CW+Freq_CW_fraction_for_fourier_plot)])
                xline((f_correction)*Freq_CW,'-m',{'RF frequency','of ac excitation'},'linewidth',3,'DisplayName','frequency of RF excitation','fontsize',lgd_font);
                %
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
                %
                subplot(2,1,2)
                plot((f_correction)*f_shift_2{L},abs(S_2{L}.^2)/max(abs(S_2{L}.^2)))
                title({['H_p_u_m_p y normalized zoomed Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
                if disp_time_femto
                    xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_pico
                    xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_nano
                    xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
                end
                ylabel('Norm. |F{H_p_u_m_p y}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
                xlim([(f_correction)*(Freq_CW-Freq_CW_fraction_for_fourier_plot) (f_correction)*(Freq_CW+Freq_CW_fraction_for_fourier_plot)])
                xline((f_correction)*Freq_CW,'-m',{'RF frequency','of ac excitation'},'linewidth',3,'DisplayName','frequency of RF excitation','fontsize',lgd_font);
                %
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
                %
                fig_indx = fig_indx +1;
                %
                % Plot the power spectrum between H_pump and J_S_in, which
                % corresponds to the difference between an RF signal that is
                % capable of switching and another one that isn't.
                %
                % BE SURE THAT H_pump and J_S_in HAVE THE SAME TIME DURATION!
                %
                figure(fig_indx)
                S_ac_initial_zero_pad=[S_ac_initial_array{L}(1,:) zero_pad_array_for_fourier] ;
                [S_11{L},f_shift_11{L}]=trans_fourier(S_ac_initial_zero_pad,dt);
                S_x_diff{L} = abs( abs(S_11{L}.^2)/max(abs(S_11{L}.^2)) - abs(S_1{L}.^2)/max(abs(S_1{L}.^2)) );
                % S_x_diff is already defined squared
                plot((f_correction)*f_shift_11{L},(S_x_diff{L})/max(abs(S_x_diff{L})))
                title({['[J_S_{-in} - H_p_u_m_p] x normalized Power Spectrum'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
                if disp_time_femto
                    xlabel('Frequency (PHz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_pico
                    xlabel('Frequency (THz)','fontsize',axis_font - smaller_font_diff)
                end
                if disp_time_nano
                    xlabel('Frequency (GHz)','fontsize',axis_font - smaller_font_diff)
                end
                ylabel('Norm. |F{J_S_{-in} - H_p_u_m_p x}|^2 (a.u)','fontsize',axis_font - smaller_font_diff)
                xline((f_correction)*Freq_CW,'-m',{'RF frequency','of ac excitation'},'linewidth',3,'DisplayName','frequency of RF excitation','fontsize',lgd_font);
                % xlim([-F_range_for_fourier_plots*enlarge_f_range F_range_for_fourier_plots*enlarge_f_range])
                xlim([(f_correction)*(Freq_CW-Freq_CW_fraction_for_fourier_plot*enlarge_f_range) (f_correction)*(Freq_CW+Freq_CW_fraction_for_fourier_plot*enlarge_f_range)])
                ax = gca;
                ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
                ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
                %
                fig_indx = fig_indx +1;
                %
            end
    end
    %
    if plot_H_pump_and_IFE
            figure(fig_indx)
%            yyaxis left
            plot(T_Global,(H_pump_array{L}(1,:)),'color','#D95319','linewidth',8,'DisplayName','H_p_u_m_p_x');
            hold on
            plot(T_Global,(H_pump_array{L}(2,:)),'color','b','linewidth',8,'DisplayName','H_p_u_m_p_y'); 
            hold on 
%            yyaxis right
%            plot(T_Global,H_IFE_array{L}(3,:),'r','linewidth',3,'DisplayName','H_I_F_E');
%            title({['H_p_u_m_p and IFE field Time Trace'],['B = {\omega}/{\gamma} + ',num2str(L -( 1 + (length(H_Xtrnl_vect) - 1)/2)),'*{\Delta}B']},'FontSize',title_font - smaller_font_diff)
            %
            if disp_time_femto
                xlabel('Time (fsec)','fontsize',axis_font)
            end
            if disp_time_pico
                xlabel('Time (psec)','fontsize',axis_font - smaller_font_diff)
            end
            if disp_time_nano
                xlabel('Time (nsec)','fontsize',axis_font - smaller_font_diff)
            end
            %           
%             ylabel(ax(1),'H_p_u_m_p (A/m)','fontsize',axis_font - smaller_font_diff)
%             ylabel(ax(2),'IFE field (A/m)','fontsize',axis_font - smaller_font_diff)
%             for z = 1:num_pulses
%                 xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-c',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse','fontsize',lgd_font); 
%             end
            lgd = legend;
            lgd.FontSize = lgd_font;
            lgd.Title.String = 'H_p_u_m_p and IFE field';
            ax = gca;
            ax.XAxis.FontSize = ax_num_font_size - smaller_font_diff; %for x-axis 
%            ax.YAxis.FontSize = ax_num_font_size - smaller_font_diff; %for y-axis 
            hold on
            % 
            fig_indx = fig_indx +1;
    end
%     if close_fig==1; close; end;
%   
    %% plot surface
    %
    for L=1:length(H_Xtrnl_vect)
        % To analyze a non-rotated component, choose {1,2,3}
        % default - 3 , the z component
        M_mat(L,:) = M_array{L}(3,:);
        %
        % To analyze a rotated component, choose {x,y,z}
        % default - M_array_Rot_z{L}(1,:) , the z component
        M_mat_Rot(L,:) = M_array_Rot_z{L}(1,:);
        %
        % To analyze the sum of S_sp & S_initial, look only at the first
        % (i.e. x) component
        % Scaling: 
        % the initial array is multiplied by Hshe_ac coversion factor,
        % S_sp is already multiplied by g_SP*(hbar/(4pi))
        S_tot_mat(L,:) = Hshe_ac*S_ac_initial_array{L}(1,:) + sign_J_s_J_sp*S_sp_array{L}(1,:);
    end 
    %
    colormap default 
    %
    if plot_surf_M
    figure(fig_indx)
    surf(T_Global,B_Xtrnl_vect',M_mat,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    %figure;h=pcolor(T_global(1:end)*1e12,Field',X_global(:,1:end));set(h,'edgecolor','none');
    %figure;contour(T_global(1:end)*1e12,Field',X_global(:,1:end))
    hold on ;  
    xline(T_end_S_ac,'-r',{'End time','of excitation'},'linewidth',3,'DisplayName','End time of excitation','FontSize',26);
    %
    if all_pls_x_lin
        for z = 1:num_pulses
            xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-r',{'H - pulse'},'linewidth',0.5,'DisplayName','H - pulse'); 
        end
    end 
    %
    if ~all_pls_x_lin
        xline(t0_H_pulsed ,'-r',{'H - pulse start'},'linewidth',3,'DisplayName','H - pulse start','FontSize',26); 
        xline(t0_H_pulsed + (num_pulses-1)*time_int_betw_pls,'-r',{'H - pulse end'},'linewidth',3,'DisplayName','H - pulse end','FontSize',26);
    end 
    %
      % Mark the original resonant field value
%         p2 = plot3([T_Global(1) T_Global(end)],[Resonanse_Field Resonanse_Field],[max(max(M_mat)) max(max(M_mat))],'m');  
%         p2.LineWidth = 2;
    %yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',2,'DisplayName','End time of excitation','FontSize',18);
    %plot3([t0_H_pulsed t0_H_pulsed],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'r')
    yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',3,'DisplayName','Original Resonanc B Field','FontSize',22,'LabelHorizontalAlignment','left');
    view(0,90)
    colorbar('FontSize',22)
    title(['NOT Normalized Magnetization Z component 3D plot VS B'],'FontSize',title_font)
    %
    if disp_time_femto
        xlabel('Time (fsec)','fontsize',axis_font)
    end
    if disp_time_pico
        xlabel('Time (psec)','fontsize',axis_font)
    end
    if disp_time_nano
        xlabel('Time (nsec)','fontsize',axis_font)
    end
    %    
    ylabel('Field (T)','FontSize',axis_font)
    zlabel('X (a.u.)')
    ax = gca;
    ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
    ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
    %plot3([t0_H_pulsed t0_H_pulsed],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'r')
%     h = text(900,100,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
%     set(h,'color',colors(L,:),'fontsize',18)
    fig_indx = fig_indx +1;
    %if close_fig==1; close; end;
    end  
    %
    if plot_surf_S_tot
        figure(fig_indx)
        up_env_tot = [];
        for L=1:length(H_Xtrnl_vect)
            [up_env, low_env] = envelope(S_tot_mat(L,:));
            up_env_tot(L,:) = up_env;
        end 
        surf(T_Global,B_Xtrnl_vect',up_env_tot,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
        %figure;h=pcolor(T_global(1:end)*1e12,Field',X_global(:,1:end));set(h,'edgecolor','none');
        %figure;contour(T_global(1:end)*1e12,Field',X_global(:,1:end))
        hold on ;  
        % Mark the end time of the J_s excitation
%         p1 = plot3([T_end_S_ac T_end_S_ac],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'m');
%         p1.LineWidth = 2;
        xline(T_end_S_ac,'-r',{'End time','of excitation'},'linewidth',3,'DisplayName','End time of excitation','FontSize',26);
        %
        if all_pls_x_lin
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-r',{'H - pulse'},'linewidth',0.5,'DisplayName','H - pulse'); 
            end
        end 
        %
        if ~all_pls_x_lin
            xline(t0_H_pulsed ,'-r',{'H - pulse start'},'linewidth',3,'DisplayName','H - pulse start','FontSize',26); 
            xline(t0_H_pulsed + (num_pulses-1)*time_int_betw_pls,'-r',{'H - pulse end'},'linewidth',3,'DisplayName','H - pulse end','FontSize',26);
        end 
        %
          % Mark the original resonant field value
    %         p2 = plot3([T_Global(1) T_Global(end)],[Resonanse_Field Resonanse_Field],[max(max(M_mat)) max(max(M_mat))],'m');  
    %         p2.LineWidth = 2;
        %yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',2,'DisplayName','End time of excitation','FontSize',18);
        %plot3([t0_H_pulsed t0_H_pulsed*],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'r')
        yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',3,'DisplayName','Original Resonanc B Field','FontSize',22,'LabelHorizontalAlignment','left');
        view(0,90)
        colorbar('FontSize',22)
        title(['S-sp + S-initial envelope'],'FontSize',title_font)
        %
        if disp_time_femto
            xlabel('Time (fsec)','fontsize',axis_font)
        end
        if disp_time_pico
            xlabel('Time (psec)','fontsize',axis_font)
        end
        if disp_time_nano
            xlabel('Time (nsec)','fontsize',axis_font)
        end
        %       
        ylabel('Field (T)','FontSize',axis_font)
        zlabel('X (a.u.)')
        ax = gca;
        ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
        ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
    %     h = text(900,100,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
    %     set(h,'color',colors(L,:),'fontsize',18)
        fig_indx = fig_indx +1;
        %if close_fig==1; close; end;
    end  
    %
    if plot_surf_S_sp_tot
        figure(fig_indx)
        up_env_tot_sp = [];
        for L=1:length(H_Xtrnl_vect)
            [up_env_sp, low_env_sp] = envelope(S_sp_array{L}(1,:));
            up_env_tot_sp(L,:) = up_env_sp;
        end 
        surf(T_Global,B_Xtrnl_vect',up_env_tot_sp,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
        %figure;h=pcolor(T_global(1:end)*1e12,Field',X_global(:,1:end));set(h,'edgecolor','none');
        %figure;contour(T_global(1:end)*1e12,Field',X_global(:,1:end))
        hold on ;  
        % Mark the end time of the J_s excitation
%         p1 = plot3([T_end_S_ac T_end_S_ac],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'m');
%         p1.LineWidth = 2;
        xline(T_end_S_ac,'-r',{'End time','of excitation'},'linewidth',3,'DisplayName','End time of excitation','FontSize',26);
        %
        if all_pls_x_lin
            for z = 1:num_pulses
                xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-r',{'H - pulse'},'linewidth',0.5,'DisplayName','H - pulse'); 
            end
        end 
        %
        if ~all_pls_x_lin
            xline(t0_H_pulsed ,'-r',{'H - pulse start'},'linewidth',3,'DisplayName','H - pulse start','FontSize',26); 
            xline(t0_H_pulsed + (num_pulses-1)*time_int_betw_pls,'-r',{'H - pulse end'},'linewidth',3,'DisplayName','H - pulse end','FontSize',26);
        end 
        %
        %yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',2,'DisplayName','End time of excitation','FontSize',18);
        %plot3([t0_H_pulsed t0_H_pulsed],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat)) max(max(M_mat))],'r')
        yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',3,'DisplayName','Original Resonanc B Field','FontSize',22,'LabelHorizontalAlignment','left');
        view(0,90)
        colorbar('FontSize',22)
        title(['S-sp envelope'],'FontSize',title_font)
        %
        if disp_time_femto
            xlabel('Time (fsec)','fontsize',axis_font)
        end
        if disp_time_pico
            xlabel('Time (psec)','fontsize',axis_font)
        end
        if disp_time_nano
            xlabel('Time (nsec)','fontsize',axis_font)
        end
        %        
        ylabel('Field (T)','FontSize',axis_font)
        zlabel('X (a.u.)')
        ax = gca;
        ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
        ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
    %     h = text(900,100,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
    %     set(h,'color',colors(L,:),'fontsize',18)
        fig_indx = fig_indx +1;
        %if close_fig==1; close; end;
    end  
    %
    %% plot extracted surface
    %
    t_xtrct = T_Global(T_Global>Ts & T_Global<Tf);
    M_mat_xtrct = M_mat(:,T_Global>Ts & T_Global<Tf) ;
    for L = 1:length(H_Xtrnl_vect)
        M_mat_xtrct_nrm(L,:) = (M_mat_xtrct(L,:)-mean(M_mat_xtrct(L,:)))/max(M_mat_xtrct(L,:)-mean(M_mat_xtrct(L,:)));
        %M_mat_xtrct_nrm(L,:)=(2*((M_mat_xtrct(L,:))-min(M_mat_xtrct(L,:)))/(max(M_mat_xtrct(L,:))-min(M_mat_xtrct(L,:))))-1;
        M_mat_nrm(L,:) = M_mat(L,:)/max(M_mat(L,:));
    end
    
%     figure;
%     surf(T_Global,B_Xtrnl_vect',M_mat_nrm,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%     hold on ;
%     plot3([t0_H_pulsed t0_H_pulsed],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat_nrm)) max(max(M_mat_nrm))],'r')
%     view(0,90)
%     colorbar
% %     %
%                 if disp_time_femto
%                     xlabel('Time (fsec)','fontsize',axis_font)
%                 end
%     if disp_time_pico
%         xlabel('Time (psec)','fontsize',axis_font)
%     end
%     if disp_time_nano
%         xlabel('Time (nsec)','fontsize',axis_font)
%     end
%     %
%     ylabel('Field (T)')
%     zlabel('X (a.u.)')
%     title('Extracted trace')
    
%     figure;
%     surf(T_Global,B_Xtrnl_vect',abs(M_mat_nrm),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%     hold on ;
%     plot3([t0_H_pulsed t0_H_pulsed],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat_nrm)) max(max(M_mat_nrm))],'r')
%     set(gca,'fontsize',14)
%     view(0,90)
%     colorbar
%     %
% if disp_time_femto
%                     xlabel('Time (fsec)','fontsize',axis_font)
%                 end
%     if disp_time_pico
%         xlabel('Time (psec)','fontsize',axis_font)
%     end
%     if disp_time_nano
%         xlabel('Time (nsec)','fontsize',axis_font)
%     end
%     %
%     ylabel('Field (T)')
%     zlabel('X (a.u.)')
%     title('Extracted trace')
%
    if plot_surf_norm_extr_M
        figure(fig_indx)
        surf(t_xtrct*1e12,B_Xtrnl_vect',((M_mat_xtrct_nrm)),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
        hold on ;
        % Mark the end time of the J_s excitation
%         p3 = plot3([T_end_S_ac T_end_S_ac],[B_Xtrnl_vect(1) B_Xtrnl_vect(end)],[max(max(M_mat_xtrct_nrm)) max(max(M_mat_xtrct_nrm))],'m');
%         p3.LineWidth = 2;
        xline(T_end_S_ac,'-r',{'End time','of excitation'},'linewidth',2,'DisplayName','End time of excitation');
        for z = 1:num_pulses
            xline(t0_H_pulsed + (z-1)*time_int_betw_pls,'-r',{'H - pulse'},'linewidth',2,'DisplayName','H - pulse'); 
        end
          % Mark the original resonant field value
%         p4 = plot3([1e12*t_xtrct(1) 1e12*t_xtrct(end)],[Resonanse_Field Resonanse_Field],[max(max(M_mat_xtrct_nrm)) max(max(M_mat_xtrct_nrm))],'m');       
%         p4.LineWidth = 2;
        yline(Resonanse_Field,'-r',{'Original Resonance','B Field'},'linewidth',2,'DisplayName','End time of excitation');
        set(gca,'fontsize',20)
        pbaspect([1.5 1 1])
        view(0,90)
        colorbar
        title({['Normalized Extracted Magnetization Z component 3D plot VS B'],['extracted in times: ',num2str(Ts*t_correction),' [nsec] to ',num2str(Tf*t_correction),' [nsec]']},'FontSize',14)
        %
        if disp_time_femto
            xlabel('Time (fsec)','fontsize',axis_font)
        end
        if disp_time_pico
            xlabel('Time (psec)','fontsize',axis_font)
        end
        if disp_time_nano
            xlabel('Time (nsec)','fontsize',axis_font)
        end
        % 
        ylabel('Field (T)')
        zlabel('X (a.u.)')
        set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
        axis([Ts*(t_correction) Tf*(t_correction) B0-B_width*half_delta_B-0.005 0.005+B0+B_width*half_delta_B]);
    %     h = text(900,100,['B = ',num2str(B_Xtrnl_vect(L)) ' T']);
    %     set(h,'color',colors(L,:),'fontsize',18)
        set(gca, 'Layer', 'top');
        grid off
        box on
        ax = gca;
        ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
        ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
        fig_indx = fig_indx +1;
    end
    % 
    % Save all figures to a desired adress
    if save_fig
        for j = fig_indx_initial:fig_indx - 1
            if save_full_screen 
                temp_fig = figure(j);
                set(temp_fig, 'Position', get(0, 'Screensize'));
                saveas(temp_fig,[adress_str,'\figure',num2str(j),'.jpg']);
            end
            %
            if save_original_size
                saveas(figure(j),[adress_str,'\figure',num2str(j),'.jpg']); 
            end
        end
    end
    fig_indx_initial = fig_indx;
    %%
    %    
    if J_ac_sweep||J_dc_sweep||H_pump_sweep||H_pump_phase_sweep||H_Gaussian_sweep||CW_polarization_sweep||CW_amplitude_sweep
        % Calculate phase between J_ac and J_sp before pulses, driven by J_ac only
        t0_H_pulsed = t0_H_pulsed*(f_correction);
        dot_product_bfr_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11),sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11));
        norm_product_bfr_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)));
        ZZ_phase_shift_in_deg_bfr_pls_sz(param_ind) = 360/(2*pi)*acos(dot_product_bfr_pls/norm_product_bfr_pls);
        % Calculate phase between J_ac and J_sp during pulses
        dot_product_aft_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22),sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22));
        norm_product_aft_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)));
        ZZ_phase_shift_in_deg_aft_pls_sz(param_ind) = 360/(2*pi)*acos(dot_product_aft_pls/norm_product_aft_pls);
        % Calculate J_SP amplitude before pulses
        ZZ_J_SP_amp_bfr(param_ind) = max(S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11));
        % Calculate J_SP amplitude during pulses
        ZZ_J_SP_amp_aft(param_ind) = max(S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22));
        % Calculate the Precession magnetization with J_ac AND PULSES infulence
        ZZ_M_z_size_min(param_ind) = min(M_array_Rot_z{L}(1,round(t0_H_pulsed_end/dt)-T_int2:round(t0_H_pulsed_end/dt)-T_int22));
        ZZ_M_z_size_max(param_ind) = max(M_array_Rot_z{L}(1,round(t0_H_pulsed_end/dt)-T_int2:round(t0_H_pulsed_end/dt)-T_int22));
        % Calculate the opening angle
        ZZ_theta_array_deg_max(param_ind) = 360/(2*pi)*acos(ZZ_M_z_size_min(param_ind)/M0);
        ZZ_theta_array_deg_min(param_ind) = 360/(2*pi)*acos(ZZ_M_z_size_max(param_ind)/M0);
        % Calculate before pulses
        ZZ_J_SP_and_Js_bfr(param_ind) = log10(max(abs(S_tot_mat(L,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)/Hshe_ac)));
        % Calculate during pulses
        ZZ_J_SP_and_Js_aft(param_ind) = log10(max(abs(S_tot_mat(L,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)/Hshe_ac)));
         %
        if CW_amplitude_sweep
            ZZ_M_z_average_size(param_ind) = 0.5*(ZZ_M_z_size_min(param_ind)+ZZ_M_z_size_max(param_ind));
            ZZ_CW_amplitude_array(param_ind) = sweeping_values(param_ind);
        end
        %
        if CW_polarization_sweep
            ZZ_M_z_average_size(param_ind) = 0.5*(ZZ_M_z_size_min(param_ind)+ZZ_M_z_size_max(param_ind));
            ZZ_CW_polarization_phase_array(param_ind) = pol_CW_phase;
        end
        %
        if J_ac_sweep
            % The amplitude vector of J_ac
            ZZ_J_ac_amp_array(param_ind) = Jc_ac;
        end
        %
        if J_dc_sweep
             % The amplitude vector of J_dc
            ZZ_J_dc_amp_array(param_ind) = Jc_dc;
            ZZ_DC_to_STT_ratio_array(param_ind) = sweeping_values(param_ind);
        end
        %   
        if H_pump_sweep
             % The amplitude vector of H_pump_0
            ZZ_B_pump_array(param_ind) = sweeping_values(param_ind);
        end
        %
        if H_pump_phase_sweep
             % The amplitude vector of H_pump_0
            ZZ_B_pump_phase_array(param_ind) = sweeping_values(param_ind);
        end
        %
        if H_Gaussian_sweep           
%             Z_torque_abs_min=min(abs(damping_torque_array{L}(3,(miu_gs/dt_low_res):85000)+primary_torque_array{L}(3,(miu_gs/dt_low_res):85000)));
%             ZZ_t_crit_array(param2_ind) = find((abs(damping_torque_array{L}(3,(miu_gs/dt_low_res):85000)+primary_torque_array{L}(3,(miu_gs/dt_low_res):85000)))==min(Z_torque_abs_min));  
            Z_torque_ratio = ((damping_torque_array{L}(3,:)+primary_torque_array{L}(3,:))./((( (damping_torque_array{L}(1,:)+primary_torque_array{L}(1,:)).^2 + (damping_torque_array{L}(2,:)+primary_torque_array{L}(2,:)).^2 + (damping_torque_array{L}(3,:)+primary_torque_array{L}(3,:)).^2  )).^0.5));
            ZZ_t_crit_array(param2_ind) = min(find(Z_torque_ratio((miu_gs/dt_low_res):end)>0.02));
            ZZ_H_crit_array(param2_ind) = H_peak*exp(-(((dt_low_res*ZZ_t_crit_array(param2_ind)).^2)./(2*(sigma_gs^2))));
            if max(abs( M_array_Rot_z{L}(1,(miu_gs/dt_low_res):T_section{3}(end)/dt_low_res)))>0.5*M0
            Mz_half = find(abs( M_array_Rot_z{L}(1,(miu_gs/dt_low_res):T_section{3}(end)/dt_low_res))>0.5*M0);
            end
            %
            if max(abs( M_array_Rot_z{L}(1,(miu_gs/dt_low_res):T_section{3}(end)/dt_low_res)))<0.5*M0
            Mz_half = 1;    
            end
            %
            if caclulate_Gamma_slope_t_half_H_half
                ZZ_t_half_array(param2_ind) = Mz_half(1);
                ZZ_H_half_array(param2_ind) = H_peak*exp(-(((dt_low_res*ZZ_t_half_array(param2_ind)).^2)./(2*(sigma_gs^2))));
                Gamma_slope_array(param2_ind) = (1/(M0*200*dt_low_res))*(M_array_Rot_z{L}(1,(miu_gs/dt_low_res)+ZZ_t_half_array(param2_ind)+100)-M_array_Rot_z{L}(1,(miu_gs/dt_low_res)+ZZ_t_half_array(param2_ind)-100));
            end
            % The amplitude vector of H_peak
            ZZ_H_peak_array(param_ind) = sweeping_values(param_ind);
             % The vector of Gaussian pulse standard deviation T_p 
            ZZ_T_p_array(param2_ind) = sweeping_values2(param2_ind);
            % The number of Gaussian pulses
            ZZ_num_pulse_array(param3_ind) = sweeping_values3(param3_ind);
            %
            ZZ_M_z_size_num_pulses(param3_ind) = M_array_Rot_z{L}(1,round(T_section{3}(end)/dt));
            if abs(ZZ_M_z_size_num_pulses(param3_ind)) > 0.9999*M0 && param3_ind<length(sweeping_values3) && num_T_p_sweeps~=1
                ZZ_M_z_size_num_pulses(param3_ind+1:end) = sign(ZZ_M_z_size_num_pulses(param3_ind))*M0;
                for j = param3_ind+1:length(sweeping_values3) 
                    ZZ_M_z_size_num_pulses(j) = sweeping_values3(j);
                end
                break
            end
            %
            ZZ_M_z_size_T_p(param2_ind) = M_array_Rot_z{L}(1,round(T_section{3}(end)/dt));
            if abs(ZZ_M_z_size_T_p(param2_ind)) > 0.9999*M0 && param2_ind<length(sweeping_values2) && num_pulses_sweeps~=1
                ZZ_M_z_size_T_p(param2_ind+1:end) = sign(ZZ_M_z_size_T_p(param2_ind))*M0;
                for j = param2_ind+1:length(sweeping_values2) 
                    ZZ_T_p_array(j) = sweeping_values2(j);
                end
                break
            end
            H_peak_Tp_matrix(param_ind,:) = ZZ_M_z_size_T_p(:);
            H_peak_num_pulses_matrix(param_ind,:) = ZZ_M_z_size_num_pulses(:);
        end
        %
        % Calculate the Precession magnetization with J_ac infulence
        ZZ_M_z_size(param_ind) = M_array_Rot_z{L}(1,round(t0_H_pulsed/dt)-T_int1);
        % Calculate the opening angle
        ZZ_theta_array_deg(param_ind) = 360/(2*pi)*acos(ZZ_M_z_size(param_ind)/M0);
        %
        ZZ_T_end_S_ac = T_end_S_ac;
        %
        ZZ_num_pulses = num_pulses;
        %
        ZZ_time_int_betw_pls = time_int_betw_pls;
        %
        ZZ_t0_H_pulsed = t0_H_pulsed;
        % Calculate the Precession magnetization with only J_ac infulence
        ZZ_M_z_size_bfr_pls = M_array_Rot_z{L}(1,round(t0_H_pulsed/dt)-T_int1);
        % Calculate the opening angle before pls
        ZZ_theta_array_deg_bfr_pls = 360/(2*pi)*acos(ZZ_M_z_size_bfr_pls/M0);
        %
        ZZ_T_end_S_ac = T_end_S_ac;
        % # pulses
        %num_pulses = ZZ_num_pulses;
        ZZ_num_pulses = num_pulses;
        % time interval between adjacent pulses
        %time_int_betw_pls = ZZ_time_int_betw_pls;
        ZZ_time_int_betw_pls = time_int_betw_pls;
        %t0 = ZZ_t0_H_pulsed;
        ZZ_t0 = t0;
        %
    end
    %
    if RF_phase_sweep
        t0_H_pulsed = t0_H_pulsed*(f_correction);
        % Calculate phase between J_ac and J_sp before pulses, driven by J_ac only
        dot_product_bfr_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11),sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11));
        norm_product_bfr_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)));
        ZZ_phase_shift_in_deg_bfr_pls_sz(param_ind) = 360/(2*pi)*acos(dot_product_bfr_pls/norm_product_bfr_pls);
        % Calculate phase between J_ac and J_sp during pulses
        dot_product_aft_pls = dot(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22),sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22));
        norm_product_aft_pls = (norm(Hshe_ac*S_ac_initial_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)))*(norm(sign_J_s_J_sp*S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)));
        ZZ_phase_shift_in_deg_aft_pls_sz(param_ind) = 360/(2*pi)*acos(dot_product_aft_pls/norm_product_aft_pls);
        % Calculate J_SP amplitude before pulses
        ZZ_J_SP_amp_bfr(param_ind) = max(S_sp_array{L}(1,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11));
        %
        % Calculate J_SP amplitude during pulses
        ZZ_J_SP_amp_aft(param_ind) = max(S_sp_array{L}(1,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22));
        %
        % Calculate the Precession magnetization with J_ac AND PULSES infulence
        ZZ_M_z_size_min(param_ind) = min(M_array_Rot_z{L}(1,round(t0_H_pulsed_end/dt)-T_int2:round(t0_H_pulsed_end/dt)-T_int22));
        ZZ_M_z_size_max(param_ind) = max(M_array_Rot_z{L}(1,round(t0_H_pulsed_end/dt)-T_int2:round(t0_H_pulsed_end/dt)-T_int22));
        %
        % Calculate the opening angle
        ZZ_theta_array_deg_max(param_ind) = 360/(2*pi)*acos(ZZ_M_z_size_min(param_ind)/M0);
        ZZ_theta_array_deg_min(param_ind) = 360/(2*pi)*acos(ZZ_M_z_size_max(param_ind)/M0);
        %
        % Calculate before pulses
        ZZ_J_SP_and_Js_bfr(param_ind) = log10(max(abs(S_tot_mat(L,round(t0_H_pulsed/dt)-T_int1:round(t0_H_pulsed/dt)-T_int11)/Hshe_ac)));
        %
        % Calculate during pulses
        ZZ_J_SP_and_Js_aft(param_ind) = log10(max(abs(S_tot_mat(L,(round(t0_H_pulsed_end/dt)-T_int2):round(t0_H_pulsed_end/dt)-T_int22)/Hshe_ac)));
        %
        % The amplitude vector of J_ac
        ZZ_RF_phase_array(param_ind) = sweeping_values(param_ind);
        %
        ZZ_T_end_S_ac = T_end_S_ac;
        %
        ZZ_num_pulses = num_pulses;
        %
        ZZ_time_int_betw_pls = time_int_betw_pls;
        %
        ZZ_t0_H_pulsed = t0_H_pulsed;
        % Calculate the Precession magnetization with only J_ac infulence
        ZZ_M_z_size_bfr_pls = M_array_Rot_z{L}(1,round(t0_H_pulsed/dt)-T_int1);
        % Calculate the opening angle before pls
        ZZ_theta_array_deg_bfr_pls = 360/(2*pi)*acos(ZZ_M_z_size_bfr_pls/M0);
        %
        ZZ_T_end_S_ac = T_end_S_ac;
        % # pulses
        %num_pulses = ZZ_num_pulses;
        ZZ_num_pulses = num_pulses;
        % time interval between adjacent pulses
        %time_int_betw_pls = ZZ_time_int_betw_pls;
        ZZ_time_int_betw_pls = time_int_betw_pls;
        %t0 = ZZ_t0_H_pulsed;
        ZZ_t0 = t0;
        %
        % LARGER ARRAYS - might be ommited for long simulations
    %     %
    %     ZZ_Mz_time_array(param_ind,:)  = M_mat(L,:);
    %     % Only envelope
    %     [up_env_sp, low_env_sp] = envelope(S_sp_array{L}(1,:));
    %     ZZ_SP_time_array(param_ind,:)  =  up_env_sp;
    %     % Only envelope
    %     [up_env_s_sp, low_env_s_sp] = envelope(S_tot_mat(L,:));
    %     ZZ_S_tot_mat_env_time_array(param_ind,:) = up_env_s_sp;
    %     % on logarithmic scale
    %     ZZ_log10_JSPJS_JS_amp_time_array(param_ind,:) = log10(abs(up_env_s_sp)/Hshe_ac);
    %     %
    %     ZZ_T_Global = T_Global;
    %     %
    end
    ZZ_T_Global = T_Global;
    ZZ_gamma_v_array = gamma_v_array{L};
    ZZ_gamma_c_array = gamma_c_array{L};
    ZZ_ReV12_array = ReV12_array{L};
    ZZ_ImV12_array = ImV12_array{L};
    ZZ_Lambda_h_array = Lambda_h_array{L};
    ZZ_Lambda_e_array = Lambda_e_array{L};
    ZZ_w_TLS_array = w_TLS_array{L};
    ZZ_gamma_inh_array = gamma_inh_array{L};
 end
end
 %
 if H_Gaussian_sweep && caclulate_Gamma_slope_t_half_H_half
    ZZ_H_crit_matrix(param_ind,:) = ZZ_H_crit_array(:); 
    ZZ_H_half_matrix(param_ind,:) = ZZ_H_half_array(:);
    ZZ_Gamma_slope_matrix(param_ind,:) = Gamma_slope_array(:);
    %
    H_crit = (1+alpha^2)*freq_aos_gs_pls/(gma*miu*alpha); 
    T_crit = ((2)^1.5)/freq_aos_gs_pls;
    t_crit = T_p*(2*log(H_peak/(0.707*H_crit)))^(0.5);
    t_half = T_p*(2*log(H_peak/(0.27*H_crit)))^(0.5);
    %
 end
 %
end

   % return
