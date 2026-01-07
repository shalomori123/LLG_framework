function  [E,H,ms_i,z1,z2,lenz,dz,w, Tsteps] = faraday_simulation_complex(e_r,sigma0,ms0,alpha, field_amplitude ...
    , fie,H_external, close_system ,NameValueArgs)
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
            ms0 (1,1) double % initial state of magnetization [A/m]
            alpha (1,1) double % parameter for the LLG equation, between 0-1
            field_amplitude (1,1) double % the amplitod of the field  [A/m]
            fie (1,1) double % phase differention of the pulse [rad]
            H_external (1,1) double %the size of the exterbal H fild [A/m]
            close_system  (1, 1) logical % True for including integration, false otherwize
            NameValueArgs.time_cut=40 % resulution parameter of how much to cut T/dt
            NameValueArgs.lamd0= 800e-9 % wavelength [m]
            NameValueArgs.len_in_lemda= 6 % length of space in wavelengths
            NameValueArgs.start_matiral= 0.2 % relative start of material
            NameValueArgs.len_matiral= 3*800e-9 % length of material [m]
            NameValueArgs.simlation_time= 1e-12 %time for the simalion to run [sec]
    end
    
    
    % phisical constants;
    eps=8.854e-12; %[F/m]
    miu=4*pi*1e-7; %[H/m]
    c=1/sqrt(eps*miu);
    
    % sistem constents 
    
    lada0= NameValueArgs.lamd0;
    dt=(1/NameValueArgs.time_cut)*(lada0/c); % the recomended sampeling freq is ten point per wave length is (1/40)*(lada0/c);
    
    dz=2*c*dt;  % the 2 factor give as 2 stap in time of singl time of dz/c
    w=2*pi*c/lada0;
    

    
    % grid and setup
    % the matrix should be as (value * z_loction * time_step)
    
    % spece defintion
    L= NameValueArgs.len_in_lemda*lada0;
    lenz=ceil(L/dz+1);
    Tsteps= ceil(NameValueArgs.simlation_time/dt);
    amplitude_vec= field_amplitude* (0.5*erf(-2:0.1:(Tsteps)/10)+0.5);

    %field start point
    H=zeros(3,lenz,Tsteps);
    E=zeros(3,lenz,Tsteps);
        
    % the material defniton
    z1=ceil(NameValueArgs.start_matiral*lenz); %location of material start
    z2=ceil(NameValueArgs.len_matiral/dz +z1); %material end
    mat_vec= logical([zeros(1,z1) ones(1,z2-z1) zeros(1,lenz-z2)]);
        
    %define eps_r
    epsr= ones(2,lenz)+(e_r-1)*mat_vec;
    
    %define condatory 
    sigma= sigma0*mat_vec;
    eaf= dt/2/eps*sigma./epsr;
    const_for_les_E= ((1-eaf)./(1+eaf));
    const_for_H= (0.5./epsr./(1+eaf));
    
     % magnetization difniton
    ms_i= zeros(3, z2-z1,Tsteps);
    ms_i(3,:,1)= sign(H_external)*((ones(z2-z1,1)* ms0)/(z2-z1))  ;
    H_exte= [0; 0; H_external];
       
    % propagate
    
    for n=1:1:2    %Duplication of for loop was done to record two steps for boundry condition
        E(:,4,n)=E(:,4,n)+ amplitude_vec(n)*[exp(1i*w*n*dt); exp(1i*(w*n*dt+fie)); 0]; %inject the input signal
    
        H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
        H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));
    
        ms_i(:,:,n+1) = runge_kutte_LLG_1(ms_i(:,:,n),real(H(:,z1+1:z2,n))+H_exte,dt,alpha); % solving LLG step
        if close_system
            H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms_i(:,:,n)- ms_i(:,:,n+1); % integration line
        end
        E(2,1:lenz-1,n+1)=const_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+const_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
        E(1,1:lenz-1,n+1)=const_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-const_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));
    end
    
    for n=3:1:Tsteps-1

       E(:,4,n)=E(:,4,n)+amplitude_vec(n)*[exp(1i*w*n*dt); exp(1i*(w*n*dt+fie)); 0]; %inject the input signal
       
       E(:,lenz,n)=E(:,lenz-1,n-2);%right boundery condition
       E(:,1,n)=E(:,2,n-2);%left boundery condition
    

       H(1,2:lenz,n+1)=H(1,2:end,n)+0.5*(E(2,2:end,n)-E(2,1:lenz-1,n));
       H(2,2:lenz,n+1)=H(2,2:lenz,n)-0.5*(E(1,2:lenz,n)-E(1,1:lenz-1,n));
    
       ms_i(:,:,n+1) = runge_kutte_LLG_1(ms_i(:,:,n),real(H(:,z1+1:z2,n))+H_exte,dt,alpha); % solving LLG step
       if close_system
           H(:,z1+1:z2,n+1) = H(:,z1+1:z2,n+1) + ms_i(:,:,n)- ms_i(:,:,n+1); % integration line
       end
       E(2,1:lenz-1,n+1)=const_for_les_E(2,1:lenz-1).*E(2,1:lenz-1,n)+const_for_H(2,1:lenz-1).*(H(1,2:lenz,n+1)-H(1,1:lenz-1,n+1));
       E(1,1:lenz-1,n+1)=const_for_les_E(1,1:lenz-1).*E(1,1:lenz-1,n)-const_for_H(1,1:lenz-1).*(H(2,2:lenz,n+1)-H(2,1:lenz-1,n+1));

    end
end






