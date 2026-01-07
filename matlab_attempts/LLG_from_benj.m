clc
clear all
% Low temporal resolution [sec]
dt_low_res= 4*4*2.5e-18;%4*0.25*1e-18;%0.005*1e-18;% 
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

% define H_optical(:,n)
% Gilbert damping 
alpha= 0.5;% 10^-4 light metals --> ~10^-2 heacy metals; 
%
gma_LL=gma./((1+alpha.^2));
LL_lamda=gma.*alpha./(1+alpha.^2);
%
% Magnetization in cgs is in [emu/cm3]. 1 [emu/cm3] --> 1*10^3 [A/m]
M0= 0.3e6;%0.3e6; %[A/m]  
%
% Create inital magnetization M_tg
%

lamda0 = 8e-7;
time_sycel = lamda0/c;
time_cut = 100;


max_feald = 6e9;
bwhp= 675e-12;
bwhp_to_us= 20*time_sycel;
[H_optical ,dt, lamda_0] = gassian_gen_new(max_feald, bwhp_to_us,pi/2, time_cut=400);
M_tg = zeros(3,length(H_optical));
M_tg(1,1) = M0;



for n = 1:(length(H_optical)-1)
                Ht_an = H_optical(:,n);      %Ht is the total field
                an= (-gma_LL*miu*cross(M_tg(:,n),Ht_an)-(LL_lamda*miu/M0)*cross(M_tg(:,n),cross(M_tg(:,n),Ht_an)));
                    
                
                Ht_bn = H_optical(:,n);
                bn=(-gma_LL*miu*cross(M_tg(:,n)+(dt/2)*an,Ht_bn)-(LL_lamda*miu/M0)*cross(M_tg(:,n)+(dt/2)*an,cross(M_tg(:,n)+(dt/2)*an,Ht_bn)));
                    
                
                Ht_cn = H_optical(:,n);
                cn=(-gma_LL*miu*cross(M_tg(:,n)+(dt/2)*bn,Ht_cn)-(LL_lamda*miu/M0)*cross(M_tg(:,n)+(dt/2)*bn,cross(M_tg(:,n)+(dt/2)*bn,Ht_cn)));
                   
                
                Ht_dn = H_optical(:,n);
                dn=(-gma_LL*miu*cross(M_tg(:,n)+dt*cn,Ht_dn)-(LL_lamda*miu/M0)*cross(M_tg(:,n)+dt*cn,cross(M_tg(:,n)+dt*cn,Ht_dn)));
                    
                
                M_tg(:,n+1) = M_tg(:,n)+(dt/6)*(an+2*bn+2*cn+dn);
end            

% T_Global is the time array
% ax_num_font_size = 20;
% 
%          figure(fig_indx)
%             plot(T_Global,(M_array_global(1,:)),'linewidth',2,'DisplayName','M_x');
%             hold on
%             plot(T_Global,(M_array_global(2,:)),'linewidth',2,'DisplayName','M_y');
%             hold on
%             plot(T_Global,M_array_global(3,:),'linewidth',4,'DisplayName','M_z');
%             hold on
%             plot(T_Global,((M_array_global(1,:)).^2 + (M_array_global(2,:)).^2 + (M_array_global(3,:)).^2).^0.5,'linewidth',1.5,'DisplayName','Magnetization Absolute Value');
%             hold off
%             title({['Magnetization components ',' Time Traces']},'FontSize',title_font)
%             %
%             if disp_time_femto
%                     xlabel('Time (fsec)','fontsize',axis_font)
%             end
%             if disp_time_pico
%                 xlabel('Time (psec)','fontsize',axis_font)
%             end
%             if disp_time_nano
%                 xlabel('Time (nsec)','fontsize',axis_font)
%             end
%             %
%             ylabel('M [A/m]','fontsize',axis_font)
%             ax = gca;
%             ax.XAxis.FontSize = ax_num_font_size; %for x-axis 
%             ax.YAxis.FontSize = ax_num_font_size; %for y-axis 
%             fig_indx = fig_indx +1;
%%
t = dt*(0:(length(H_optical)-1));
M_tg = real(M_tg);
m_tg_size = sqrt(M_tg(1,:).^2 + M_tg(2,:).^2 + M_tg(3,:).^2);
hold on
plot(t, squeeze(real(M_tg(1,:))), DisplayName="x")
plot(t, squeeze(real(M_tg(2,:))), DisplayName="y")
plot(t, squeeze(real(M_tg(3,:))), DisplayName="z")
plot(t, squeeze(m_tg_size), "k", DisplayName="size")

legend()

