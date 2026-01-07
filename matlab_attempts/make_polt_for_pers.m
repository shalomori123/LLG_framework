path= 'C:\maxwell-LLG\ws\simulation_result\Full_Width_Pulse\';
file = 'summrize.npy';


data = load_np([path file]);
data = data{1};
%%
%%alpha = {"Co_Au" : 0.02, "Co_Pt" : 0.025 , "Fe_Au": 0.02, "Fe_Pt" : 0.025, "Ni_Au": 0.05, "Ni_Pt" : 0.06}

magntic_name = ["Co", "Ni", "Fe"];
conductive_name = ["Au", "Pt"];
ms0 = 3e5;
open_mz = sum(data.mz(1:2:12,:),2)/ms0;
open_my = sum(data.my(1:2:12,:),2)/ms0;
a = [0.02 0.025 0.02 0.025 0.05 0.06];

subplot(2,1,1)
plot(a, open_mz,"*")
grid on
title("Total moment effect M_z")
xlabel("Gilbert Damping")
ylabel("M_z/M")
xlim([0.018 0.061])
ylim([-3.6 -2]*1e-7)


subplot(2,1,2)
plot(a, open_my,"*")
grid on
title("Total moment effect M_y")
xlabel("Gilbert Damping")
ylabel("M_y/M")
xlim([0.018 0.061])
ylim([-1.4 -0.8]*1e-6)


