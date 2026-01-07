%%
clc
clear all
close all

% phisical constants;
eps=8.854e-12; %[F/m]
miu=4*pi*1e-7; %[H/m]
C=1/sqrt(eps*miu);
heta0 = sqrt(miu/eps);
%% gassian gen verfiction
%- check power in
% check 1- p =sum (heta0 * E^2)

Pin_vec = 2000:5000:100000;
P_gen = zeros(size(Pin_vec));
for i= 1:length(Pin_vec)
    [pulse ,dt, lamda0] = gassian_gen(Pin_vec(i), 0);
    P_gen(i) = heta0 * sum(sum(abs(pulse).^2))/2;
end
check1 = (Pin_vec - P_gen)./Pin_vec;
check1_mean_error = mean(check1)

%check 2 - does polarization change the power?
%          checks for what polriztion we make
polarization = linspace(-pi/2, pi/2,9);
P_gen = zeros(size(polarization));
Powerin= 100000;
hold on
for i= 1:length(polarization)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, polarization(i));
    P_gen(i) = heta0 * sum(sum(abs(pulse).^2))/2;
    fv= [pulse(1,:);pulse(2,:)];
    [tau,epsilon] =polellip(fv);
    t= dt * (1:length(pulse));
    plot(t,epsilon, "DisplayName" , num2str(i))

end
check2 = (Powerin - P_gen)/Powerin; % no :-)
check2_mean_error = mean(check2)


% check 3 - does the pulse_width change the power?
pulse_width = (0.2* 2.6685e-14): 0.5*2.6685e-14:2.6685e-13;
P_gen = zeros(size(pulse_width));
Powerin= 100000;
for i= 1:length(pulse_width)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, 0,"pulse_width",pulse_width(i));
    P_gen(i) = heta0 * sum(sum(abs(pulse).^2))/2;
end
check3 = (Powerin - P_gen)/Powerin; % not a lot, looks random :-)
check3_mean_error = mean(check3)
plot(pulse_width, check3)
%looks good, everything works


% check 4 - does the time_cut change the power?
time_cut =40:10:120;
P_gen = zeros(size(time_cut));
Powerin= 100000;
for i= 1:length(time_cut)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, 0,"time_cut",time_cut(i));
    P_gen(i) = heta0 * sum(sum(abs(pulse).^2))/2;
end
check4 = (Powerin - P_gen)/Powerin; % not a lot looks random :-)
check4_mean_error = mean(check4)
% not sagnificant.
%% void siumllation verifction 
% vector pointing anlize

eps_r = 1*  [1, 1] ;
sigma =0* [1 1];


ms0 = [0 0 1];
alpha =0.0001;
H_ext = [0 0 0];
close_system = 0 ;
edge_reflection = 1;

% - connection beetwin pointing vec-power.
figure
% check1 - power in effact on S
Pin_vec = 2000:2000:10000;
P_comp = zeros(size(Pin_vec));
hold on
for i= 1:length(Pin_vec)
    [pulse ,dt, lamda0] = gassian_gen(Pin_vec(i), pi/2);
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);
    E_real= heta0*E;
    S = cross(real(E_real),real(H));
    S_power = squeeze(sum(S(3,:,:),3));
    Z= (1:size(E,2));
    plot(Z,S_power, "DisplayName" , num2str(i))
    P_comp(i) = Pin_vec(i) - mean(S_power(5:end));
end
cheke1 = P_comp./Pin_vec 
% conclusion - const ratio
%%


% check2 - pulse width efact on S
pulse_width = linspace(0.05* 2.6685e-14, 2.6685e-14, 20);

P_comp = zeros(size(pulse_width));
Powerin  = 10^6;
for i= 1:length(pulse_width)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, 0,"pulse_width",pulse_width(i));
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);
    E_real= heta0*E;
    S = cross(real(E_real),real(H));
    S_power = squeeze(sum(S(3,:,:),3));
    Z= (1:size(E,2));
    P_comp(i) = Powerin - mean(S_power(5:end));
end
cheke2 = P_comp/Powerin
figure()
plot(10*pulse_width/(2.6685e-14),cheke2)
xlabel("pulse width [time cycle]")
ylabel("relative mistake") 
% conclusion - stablizes wider then 4 time cycles


% check3 - time cut effect on S
time_cut =80:10:160;

P_comp = zeros(size(time_cut));
Powerin  = 10^6;
for i= 1:length(time_cut)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, 0,"time_cut",time_cut(i));
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);
    E_real= heta0*E;
    S = cross(real(E_real),real(H));
    S_power = squeeze(sum(S(3,:,:),3));
    Z= (1:size(E,2));
    P_comp(i) = Powerin - mean(S_power(5:end));
end

cheke3 = P_comp/Powerin;

figure("Name","time cut efact on pointing_vec")
plot(log(time_cut),log(abs(cheke3)))
t_log = log(time_cut);
cheke3_log = log(abs(cheke3));
xlabel("time cut")
ylabel("relative mistake") 
% conclusion: log (relative mistake0) = -1.967* log( time_cut) +3.501


%% void siumllation verifction 
% measure wave len verfiction 

% using the calculations:
% - mean in time wave len---- look for when relevant
% - mean in space wave len- inside the measure_wave_length finction

[pulse ,dt, lamda0] = gassian_gen(Powerin, 0);
[E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);

E_real= heta0*E;
S = cross(real(E_real),real(H));
power_envlop = sum(heta0 * (abs(pulse).^2)/2,1);


% - when can we measure wave len-
%% check 1 - standart pulse using graph
Pin = 3000;
wave_len= zeros(1,size(E,3));
for n= 1: size(E,3)
    wave_len(n) = measure_wave_length(real(E(:,:,n)), [6 size(E,2)-1]);
end
t= 1: (size(E,3));
plot(t, wave_len)

% find stability indexes by graph

stabl_lim = [257 2144]; % ~ [exp(-5.5) exp(-8)] of max
very_stabl = [1390 1981]; % ~ [after_max  exp(-5)] of max
figure("Name","comper")
hold on
for ind = very_stabl
    plot(real(squeeze(E(1,:,ind))),DisplayName= num2str(ind))
end
hold off

% finding the power of S that should move throgh relevant indexes 
max_S = max(max(S(3,end,:)));
stabl_S_pass = [max(S(3,end, 258-10:258+10)) max(S(3,end,2144-10:2144+10))]/ max_S ; % s val 
vary_stabl_S_pass = [max(S(3,end, 1390-10:1390+10)) max(S(3,end,1981-10:1981+10))]/ max_S ; % noat the first is aftar max(S) 

ind_first = find(S(3,end,:) >= vary_stabl_S_pass(1)* max_S, 1, "last");
ind_last = find(S(3,5,:) >= vary_stabl_S_pass(2) * max_S, 1, "last");
posbel_lim = [ ind_first ind_last]

% conclusion: - measure most to times good between last 0.7*max_S to
%              last 0.023*max_S
%            - first good measure at first 4e-3*max_S and last good measure at last 4e-3*max_S
            

%% check 2 - width effact on first and last good measures

Pin = 10e5;
pulse_width = linspace(1,2, 12);
pass_first = zeros(1,12);
pass_last = zeros(1,12);
betwin_times = zeros(1,12);
for i = 1:12
    [pulse ,dt, lamda0] = gassian_gen(Pin, 0, "pulse_width", pulse_width(i)*2.6685e-15);
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection);
    wave_len= zeros(1,size(E,3));
    for n= 1: size(E,3)
        wave_len(n) = measure_wave_length(real(E(:,:,n)), [6 size(E,2)-1]);
    end
    ind = [ find(wave_len == 20,1,"first") find(wave_len == 20,1,"last")];
    if length(ind) == 2 
        
        E_real= heta0*E;
        S = cross(real(E_real),real(H));
        S_max= max(S(3,90,:));
        pass_first(i) = max(S(3,end, ind(1):ind(1)+20))/S_max;
        pass_last(i) = max(S(3,6, ind(2)-20:ind(2)))/S_max;
        betwin_times(i) = ind(2)- ind(1);
    end
end

hold on
plot(pulse_width, pass_last, DisplayName= "last")
plot(pulse_width, pass_first, DisplayName= "first")

% plot(pulse_width, betwin_times, DisplayName= "chak frim")

% conclusion: succeed to measure correctly from  width of 1.46 time cycal.
%             needed 0.12*max(S_z) to first good measure, can be good down
%             to 1e-3max(S_z)

%% dialctricl simllation varifction
% - conction btiwin S to power in and pulse width  .
% - when can we maserr wave len-  inside and outside matial
% - vectoer pointing retrun and tarnsport: - mache to theory
%                                          - end matiral effcts


eps_r = 4*  [1, 1] ;
sigma =0* [1 1];


ms0 = [0 0 1];
alpha =0.0001;
H_ext = [0 0 0];
close_system = 0 ;
edge_reflection = 1;


% - connection bitween S to power in and pulse width  .

% check1 - power in effact on S
Pin_vec = 2000:2000:10000;
S_pro_matiral = zeros(size(Pin_vec));
S_in_matiral = zeros(size(Pin_vec));
S_after_matiral = zeros(size(Pin_vec));

hold on
for i= 1:length(Pin_vec)
    [pulse ,dt, lamda0] = gassian_gen(Pin_vec(i), pi/2, "relation_time", 30*2.6685e-14 );
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection, ...
        "len_matiral",3*800e-9, "start_matiral",3*800e-9, "len_simulation",7*800e-9);
    E_real= heta0*E;
    S = cross(real(E_real),real(H));
    S_power = squeeze(sum(S(3,:,:),3));
    t= (1:size(E,2));
    plot(t,S_power, "DisplayName" , num2str(i+5))
    S_pro_matiral(i) = mean(S_power(z1-50:z1));
    S_in_matiral(i) = mean(S_power(z2-50: z2));
    S_after_matiral(i) = mean(S_power(z2+2:end));
    
    
end

cheke1 = [Pin_vec; S_pro_matiral; S_in_matiral; S_after_matiral] 

plot([z1 z1], [-1.2e4, 1.2e4],"k")
plot([z2 z2], [-1.2e4, 1.2e4],"k")
% conclusion: - before the material we can see a standing wave that returns
%               from the material. lower then inside material due to
%               reflections.
%               after material - const value as expected.
%               inside material - also stading wave according to internal
%               reflection from the edges.

%%
% chack 2 - pulse width effect on S
pulse_width = linspace(0.05* 2.6685e-14, 5*2.6685e-14, 10);

S_pro_matiral = zeros(size(pulse_width));
S_in_matiral = zeros(size(pulse_width));
S_after_matiral = zeros(size(pulse_width));

Powerin  = 10^6;
for i= 1:length(pulse_width)
    [pulse ,dt, lamda0] = gassian_gen(Powerin, 0,"pulse_width",pulse_width(i));
    [E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection, ...
        "len_matiral",3*800e-9, "start_matiral",3*800e-9, "len_simulation",7*800e-9);
    E_real= heta0*E;
    S = cross(real(E_real),real(H));
    S_power = squeeze(sum(S(3,:,:),3));
    Z= (1:size(E,2));
    S_pro_matiral(i) = mean(S_power(z1-50:z1));
    S_in_matiral(i) = mean(S_power(z2-50: z2));
    S_after_matiral(i) = mean(S_power(z2+2:end));
end


cheke2 = [ S_pro_matiral; S_in_matiral; S_after_matiral] 
plot(S_pro_matiral)
hold on
plot(S_after_matiral)

% plot([z1 z1], [-1.2e4, 1.2e4],"k")
% plot([z2 z2], [-1.2e4, 1.2e4],"k")
%conclusion:
% Warning!!!! it looks like the amount of energy change - S_pro and S_after behave the same way
%              - short pulses have more frequencies and different
%              reflections.
%              - again, S inside the material is bigger then pro material due
%              to reflections.


%%
% check 3 - long z distent

eps_r = 4*  [1, 1] ;
sigma =0* [1 1];
ms0 = [0 0 1];
alpha =0.0001;
H_ext = [0 0 0];
close_system = 0 ;
edge_reflection = 1;
T= 2.6685e-15;
pulse_width= 3*T;
lamda0= 8e-7;
start_matiral= 7* (pulse_width/T) * lamda0;
len_matiral = 20* lamda0;
len_simulation = len_matiral + start_matiral + lamda0;
Powerin  = 10^6;

[pulse ,dt, lamda0] = gassian_gen(Powerin, pi/2, "pulse_width",pulse_width);
[E, H, ms_i, z1, z2, landz, dz, Tsteps] = simulation(eps_r,sigma, ms0, alpha, H_ext, pulse, dt, close_system, edge_reflection, ...
    "len_matiral",len_matiral, "start_matiral",start_matiral, "len_simulation",len_simulation);
E_real= heta0*E;
S = cross(real(E_real),real(H));

pro = squeeze(S(3,6,:));
in_matrial = squeeze(S(3,z1+40,:));


ind_all_pass= find(in_matrial > 0, 1,"first");
arrival_time_to_matiral = 2 * ceil(start_matiral/dz); %% conset for void
passing_material_time = 4* ceil(len_matiral/dz - 40); %% factor 4 dipend on the matrial..
time_to_pass_before = ceil(5*pulse_width/dt) + 2 * ceil(start_matiral/dz);
S_in = sum(pro(1:arrival_time_to_matiral));
S_return = sum(pro(arrival_time_to_matiral:arrival_time_to_matiral+ passing_material_time));
S_pass = sum(in_matrial(1:arrival_time_to_matiral+ passing_material_time));
Transmition = S_pass/ S_in   
Return = abs(S_return/ S_in)
check1 = (S_in - S_pass + S_return)/Powerin
check2 = Transmition + Return

hold on
plot(squeeze(S(3,5,:)),DisplayName = "befor")
plot(arrival_time_to_matiral,0,"ok")

plot(arrival_time_to_matiral + passing_material_time,0,"or")
plot(arrival_time_to_matiral + 2* passing_material_time,0,"ob")
plot(squeeze(S(3,z1,:)),DisplayName = "matiral")
plot(squeeze(S(3,z2+1,:)),DisplayName = "after")
% conclusion - we able to measure power Transmition and Return when we have
% long enough void-space to seperate the generated pulse from the returnd
% pulse.
% In this script we define the start location of material to fit for
% these condition, and find the indexes to execute this seperation.
%%
hold on
plot(squeeze(sum(S(3,:,:),3)))
plot([z1 z1], [-1.2e6, 1.2e6],"k")
plot([z2 z2], [-1.2e6, 1.2e6],"k")

%%


lim=max(max(max(real(E))));
%Y_axis=5;

[D ,I, N]=size(E);
z=0:dz:dz*(I-1);
n= 600;
[x_s, y_s] = meshgrid(-lim:lim/3:lim);
z_s = dz*meshgrid(z1:(z2-z1)/6:z2);
x1_s= -ones(size(y_s))*lim;
z2_s= z2* ones(size(x_s))*dz;
plot3(z,E(1,:,n),E(2,:,n));

hold on
surf(z_s, x1_s,y_s, 'FaceAlpha',0.3 , 'EdgeColor', 'none', 'FaceColor', 'r');
% surf(z2_s,x_s,y_s, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');
axis([0 z(end) -1.1*lim 1.1*lim -1.1*lim 1.1*lim]);
               
xlabel('z')
ylabel('E_x')
zlabel('E_y')
            
            
           
%% condctivty simllation vatifction
% - all as dialctricl
% - בעלית אנרגיה לפי מרחק- איך למדוד ואם מתיאים לתאוריה של עומק חדירה
% - 

%% dialctricl condctivty simulation


%% open systme verfiction (to benj)
% - Dune
