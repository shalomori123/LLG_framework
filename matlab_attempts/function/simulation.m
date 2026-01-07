function [E,H,ms_i,z1,z2,lenz,dz, Tsteps]  = simulation(e_r,sigma0,ms0, alpha, H_ext, Pulse, dt, close_system, edge_reflectance ,NameValueArgs)

    arguments
            e_r (2,1) double % 2 dimensional epsilon r of the material 
            sigma0 (2,1) double % 2 dimensional conductivity of the material [S/m]
            ms0 (1,3) double % initial state of magnetization [A/m]
            alpha (1,1) double % parameter for the LLG equation, between 0-1
            H_ext (1,3) double % 
            Pulse double 
            dt (1,1) double %dt  [sec]
            close_system  (1, 1) logical % True for including integration, false otherwize
            edge_reflectance (1, 1) logical % True for including reflectance at edge, false otherwize
            NameValueArgs.len_simulation= 5*800e-9 % length of space [m]
            NameValueArgs.start_matiral= 1.5*800e-9 % start of material [m]
            NameValueArgs.len_matiral= 2*800e-9 % length of material [m]
    end
    % phisical constants;
    eps=8.854e-12; %[F/m]
    miu=4*pi*1e-7; %[H/m]
    c=1/sqrt(eps*miu);
    dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c

    Tsteps = length(Pulse);
      
    % spece defintion
    L= NameValueArgs.len_simulation;
    lenz=ceil(L/dz+1);

    %field start point
    H=zeros(3,lenz,Tsteps);
    E=zeros(3,lenz,Tsteps);
  

   % the material defniton
    z1=ceil(NameValueArgs.start_matiral/dz + 1); %location of material start
    z2= min([ceil(NameValueArgs.len_matiral/dz +z1), lenz]); %material end
    if edge_reflectance
        mat_vec= ([zeros(1,z1) ones(1,(z2- z1)) zeros(1,lenz-z2)]);
    else
        finsh_len = ceil(0.8 * (z2- z1));
        finsh_vec = 0.5 *(1 - erf(linspace(-2.7, 2.7,finsh_len)));
        mat_vec= ([zeros(1,z1) ones(1,(z2- z1) -finsh_len) finsh_vec zeros(1,lenz-z2)]);
    end
    
       
    %define eps_r
    epsr= ones(2,lenz)+(e_r-1)*mat_vec;
    
    %define condatory 
    sigma= sigma0*mat_vec;
    eaf= dt/2/eps*sigma./epsr;
    conset_for_les_E= ((1-eaf)./(1+eaf));
    conset_for_H= (0.5./epsr./(1+eaf));
    
     % magnetization difniton
    ms_i= zeros(3, z2-z1,Tsteps);
    ms_i(:,:,1)= ((ones(z2-z1,1).* ms0)/(z2-z1))'  ;
       
    % propagate
    
    for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
        E(:,4,n)=E(:,4,n)+Pulse(:,n); %inject the input signal
    
        H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
        H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));
    
        ms_i(:,:,n+1) = runge_kutte_LLG_1(ms_i(:,:,n),H(:,z1+1:z2,n) + H_ext' ,dt,alpha); % solving LLG step
        if close_system
            H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms_i(:,:,n)- ms_i(:,:,n+1); % integration line
        end
        E(2,1:lenz-1,n+1)=conset_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+conset_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
        E(1,1:lenz-1,n+1)=conset_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-conset_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
    end
    
    for n=3:1:Tsteps-1

       E(:,4,n)=E(:,4,n)+Pulse(:,n); %inject the input signal
       
       E(:,lenz,n)=E(:,lenz-1,n-2);%right boundery condition
       E(:,1,n)=E(:,2,n-2);%left boundery condition
    

       H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
       H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));
    
       ms_i(:,:,n+1) = runge_kutte_LLG_1(ms_i(:,:,n),H(:,z1+1:z2,n)+ H_ext' ,dt,alpha); % solving LLG step
       if close_system
           H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms_i(:,:,n)- ms_i(:,:,n+1); % integration line
       end
       E(2,1:lenz-1,n+1)=conset_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+conset_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
       E(1,1:lenz-1,n+1)=conset_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-conset_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));

    end
    Tsteps= Tsteps-1;

end