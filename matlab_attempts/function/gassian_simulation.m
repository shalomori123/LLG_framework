function  [E,H,ms_i,z1,z2,lenz,dz,W, Tsteps] = gassian_simulation(e_r,sigma0,ms0,alpha, max_field ...
    , pulse_width, fea,H_ext, close_system ,NameValueArgs)
    % Running one simulation of maxwell equation integrated with LLG
    % equations, using a gaussian pulse according to the paramters:
    % returns:
    % E,H,ms - as matrix (value [A/m] * z_loction * time_step) (E is normalized)
    % z1,z2- start and end indxes of matiral
    % lenz- amunt of disctizition len unit
    % dz - len unit space

    arguments
            e_r (2,1) double % 2 dimensional epsilon r of the material 
            sigma0 (2,1) double % 2 dimensional conductivity of the material [S/m]
            ms0 (1,3) double % initial state of magnetization [A/m]
            alpha (1,1) double % parameter for the LLG equation, between 0-1
            max_field (1,1) double % the maximum field of the pulse, [A/m]
            pulse_width (1,1) double % pulse width in cycle number the timecycal
            fea (1,1) double % phase differention of the pulse [rad]
            H_ext (1,3) double % 
            close_system  (1, 1) logical % True for including integration, false otherwize
            NameValueArgs.time_cut=40 % resulution parameter of how much to cut T/dt
            NameValueArgs.lamd0= 800e-9 % wavelength [m]
            NameValueArgs.len_simulation= 2*800e-9 % length of space [m]
            NameValueArgs.start_matiral= 800e-9 % start of material [m]
            NameValueArgs.len_matiral= 80e-9 % length of material [m]
            NameValueArgs.teta= 0 % initial angle of the polarization [rad]
    end
    
    
    % phisical constants;
    eps=8.854e-12; %[F/m]
    miu=4*pi*1e-7; %[H/m]
    c=1/sqrt(eps*miu);
    
    % sistem constents 
    
    lada0= NameValueArgs.lamd0;
    dt=(1/NameValueArgs.time_cut)*(lada0/c); % the recomended sampeling freq is ten point per wave length is (1/40)*(lada0/c);
    
    dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c
    W=2*pi*c/lada0;
    

    
    % grid and setup
    % the matrix should be as (value * z_loction * time_step)
    
    % spece defintion
    L= NameValueArgs.len_simulation;
    lenz=ceil(L/dz+1);

    % pulse defintion
    pulse_width_us= pulse_width*2*pi/W;
    Tsteps= ceil(4*lenz+(5*pulse_width_us/dt)) ;
    peak_location =2.5*pulse_width_us/dt/Tsteps;
    
    Pulse =create_gassian_envelope(Tsteps,dt,pulse_width_us,max_field,peak_location);
    Pulse= amplified(W,dt,Pulse,fea, NameValueArgs.teta) ;
    %plot(Pulse)
    Pulse= [Pulse(:,1)';Pulse(:,2)'; zeros(1,Tsteps+1)];

    %field start point
    H=zeros(3,lenz,Tsteps);
    E=zeros(3,lenz,Tsteps);
        
    % the material defniton
    z1=ceil(NameValueArgs.start_matiral/dz + 1); %location of material start
    z2= min([ceil(NameValueArgs.len_matiral/dz +z1), lenz]); %material end
    mat_vec= logical([zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)]);
        
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






