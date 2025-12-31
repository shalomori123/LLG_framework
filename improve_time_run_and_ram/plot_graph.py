import numpy as np
import matplotlib.pyplot as plt
import simulation_cpu  as sim
import scipy.special as scp
import analyze_tool as an
from tqdm import tqdm

FWHM = 0.22e-12 #[sec]
# B_opt =2.4e5 #[A/m ] = 0.3[T]

# w= 2*np.pi/(mean_wave_lenght/sim.c)

# path = "C:\maxwell-LLG\ws\final_result\Return_Transmision_Included\\"
path = "C:\maxwell-LLG\ws\simulation_result\Return_Tranmision_Incloud\\"
dz = 2e-9/2

magntic_steps = 10e-9//dz
conductive_steps = 2e-9//dz


# pulse properties
sefe_start = 3
sigma_pulse = ((0.5/np.log(2))** 0.5) * FWHM # עם זה הרצנו את הכל
# sigma_pulse = 0.5* (np.log(2)**-0.5) * FWHM # זה הגרסא הנכונה לפי הנוסחא של FWHM
peak_enter_time= sefe_start * sigma_pulse
fie = np.pi/2

#material properties
magntic_name = ["Co", "Fe", "Ni"]
conductive_name = ["Au", "Pt"]
S_transmission = np.zeros((6,2))
S_return = np.zeros((6,2))
S_pass = np.zeros((6,2))
S_after_material = np.zeros((6,2))
S_bitween_materials = np.zeros((6,2))
i =0

for mag in range(3):
    for cond in range(2):
        for close_sistem in [False, True]:
            file_name = path + magntic_name[mag] +"_" + conductive_name[cond]
            if close_sistem:
                file_name += "_close.npy"
                c =0
            else:
                file_name+= "_open.npy"
                c = 1
            
            data = np.load(file_name, allow_pickle=True)
            data = data.item()
            data_dt = data["dt"]* data["time intervals"]
            matiral_l = data["matitral location"]
            save_indxes = list(data["locations"])


            z_enter = 0
            z_matirel = save_indxes.index(matiral_l[0])
            z_out = save_indxes.index(matiral_l[-1]+1)
            z_end_magnet = save_indxes.index(matiral_l[1])
            # print(f"exist indexes :{save_indxes[:]}")
            # print(f"using numbers : {[z_enter, z_matirel, z_out, z_end_magnet]}")
            
            E_file= data["E in locations"]

            
            void_steps = 2* data["matitral location"][0]/data["time intervals"]
            t1 = int(peak_enter_time/data_dt + void_steps)

            t_hm = int((sefe_start* 0.5* (np.log(2)**-0.5) * FWHM)/data_dt)
 
            E = np.zeros((np.size(E_file, 0), np.size(E_file,1), 3)) # (time, loction, dimantion)
            H = np.zeros((np.size(E_file, 0), np.size(E_file,1), 3)) # (time, loction, dimantion)

            E[:,:,0:2] = sim.heta* np. real(E_file)
            H[:,:,0:2] = np.real(data["H in locations"])

            S = np.cross(np.real(E),np.real(H))

            S_in_enter_point = S[:,z_enter,2]
            S_in_z1 = S[:,z_matirel,2]
             
            time = np.arange(0, np.size(E_file,0)) *data_dt




            S_transmission[i,c] = sum(S_in_enter_point[:t1])*data_dt
            S_return[i,c] = sum(S_in_enter_point[t1:])*data_dt
            S_pass[i,c] = sum(S_in_z1)*data_dt
            S_after_material[i,c] = sum(S[:,z_out,2])*data_dt
            S_bitween_materials[i,c] = sum(S[:,z_end_magnet,2])*data_dt
                        
        i += 1


transmition = S_pass/ S_transmission
Return = S_return/ S_transmission
print(-Return.T)
print(transmition.T)




####################### 
# er = [1.5, 1.6]
# sigma = [3000, 4000]
# B_opt =5*2.4e5 #[A/m ] = 0.3[T]



# mean_wave_lenght= 784e-9 # [m]
# FWHM = 0.4e-12 # sec
# w= 2*np.pi/(mean_wave_lenght/sim.c)


# # besic stuff
# time_cycal= mean_wave_lenght/sim.c
# dz = mean_wave_lenght/40
# dt = dz/(2*sim.c)

# magntic_steps = 1.5*mean_wave_lenght//dz
# conductive_steps = 0.5*mean_wave_lenght//dz

# sefe_start = 3
# sigma_pulse = 0.5* (np.log(2)**-0.5) * FWHM
# peak_enter_time= sefe_start * sigma_pulse
# fie = np.pi/2


# ms0 = 3e6*np.array([1, 0, 0 , 0 ,0, -0.0001]).reshape((2,3))

# z1 = 24
# z_end = int(z1+ magntic_steps + conductive_steps)
# z_indexes = [z1, int(z1+ magntic_steps), z_end ]

# void_steps = 2*z1

# exit_steps = int(time_cycal/ dt)


# d ={"z_ind" : z_indexes}

# # we run in sorter time but this what we need:
# simulation_time = int(3*(peak_enter_time//dt + void_steps +4* (magntic_steps + conductive_steps)) + exit_steps)

# save_loc = []
# save_time = np.linspace(simulation_time*0.5,simulation_time*0.9, 20).astype("int")

# data_open = sim.simulation(z_indexes, er , sigma, 0.5 , ms0,
#         B_opt, fie , FWHM,
#         save_loc, save_time, False,
#         dt, simulation_time, lamda=mean_wave_lenght, safety_start= sefe_start, time_interval= 7 )


# d["H_open"] =np.real(data_open["H times picture"])

# data_close = sim.simulation(z_indexes, er , sigma, 0.1 , ms0,
#         B_opt, fie , FWHM,
#         save_loc, save_time, True,
#         dt, simulation_time, lamda=mean_wave_lenght, safety_start= sefe_start, time_interval= 7 )



# d["H_close"] = np.real(data_close["H times picture"])


# np.save(path + "for_matlav.npy", d)