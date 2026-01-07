function [magnetization] = runge_kutte_LLG_1(M, H, dt, alpha)
%This function runs 1 step of the runge-kutte algorithm, and returns the
%magnetization 

    q=1.60217646e-19;                        %[cb]
    miu=4*pi*1e-7;                           %[H/m]
    g=2;                                     %[Landau factor]
    me=9.1093821545*1e-31;                   %[gram] electron mass
    gma_factor=1;
    gma = gma_factor*g*q/(2*me);               %[rad/sec*T]

    
    gma_LL=gma./((1+alpha.^2));
    LL_lamda=gma.*alpha./(1+alpha.^2);
    %
    % Magnetization in cgs is in [emu/cm3]. 1 [emu/cm3] --> 1*10^3 [A/m]
    M0= sqrt(sum(M.^2)) ;% 0.3e6;%0.3e6; %[A/m]  
    %
    % Create inital magnetization M_tg
    

   an = -gma_LL*miu*cross(M,H)-(LL_lamda*miu./M0).*cross(M,cross(M,H));
        
   bn = -gma_LL*miu*cross(M+(dt/2)*an,H)-(LL_lamda*miu./M0).*cross(M+(dt/2)*an,cross(M+(dt/2)*an,H));
        
   cn= -gma_LL*miu*cross(M+(dt/2)*bn,H)-(LL_lamda*miu./M0).*cross(M+(dt/2)*bn,cross(M+(dt/2)*bn,H));
       
   dn=-gma_LL*miu*cross(M+dt*cn,H)-(LL_lamda*miu./M0).*cross(M + dt*cn,cross(M +dt*cn,H));
        
    
   magnetization = M + (dt/6)*(an+2*bn+2*cn+dn);




end