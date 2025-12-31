import numpy as np
import matplotlib.pyplot as plt
import improve_time_run_and_ram.simulation_cpu_30_10  as sim
import scipy.special as scp
import analyze_tool as an

# magnets matirals proprety

# cobalt
Co_er = 1
Co_sigma = 1.7e7 #https://periodictable.com/Elements/027/data.html
Co_er_complex= -16 + 23.3j 

# Iron
Fe_er = 1
Fe_sigma = 1e7

# nical
Ni_er = 1
Ni_sigma = 1.4e7 #https://periodictable.com/Elements/027/data.html 

alpae_renge = (0.025 ,0.0025) #genral range
gilbert_damping_time = 600e-12 #[sec] artical - concting to alpha in {I dont remmber how}
precession_frequency = 8.5e9 #[Hz] artical

# aditinal matirals proprety
Pt_er = 1
Pt_sigma = 7e6 # [1/oham*m]

Au_er = 1
Au_sigma = 1e6 # (spintronics lab)


alpha = {"Co_Au" : 0.02, "Co_Pt" : 0.025 , "Fe_Au": 0.02, "Fe_Pt" : 0.025, "Ni_Au": 0.05, "Ni_Pt" : 0.06}

# pulse proprty (in the article)
mean_wave_lenght= 784e-9 # [m]
# B_opt =  [0.02, 0.2 ] # [T]
# FWHM = 1.1e-12 #[sec]


FWHM = 3 *(mean_wave_lenght/sim.c) #[sec]
B_opt =2.4e5 #[A/m ] = 0.3[T]

w= 2*np.pi/(mean_wave_lenght/sim.c)


path = "C:\maxwell-LLG\ws\simulation_result\Full_Width_Pulse\\"




# besic stuff
lamda = 8e-7 
time_cycal= mean_wave_lenght/sim.c
dz = 2e-9/8
dt = dz/(2*sim.c)

# matiral length according to arcticle
magntic_steps = 10e-9//dz
conductive_steps = 2e-9//dz

# pulse properties
sefe_start = 3
sigma_pulse = 0.5* (np.log(2)**-0.5) * FWHM
peak_enter_time= sefe_start * sigma_pulse
fie = np.pi/2

#material properties
magntic_name = ["Co", "Fe", "Ni"]
conductive_name = ["Au", "Pt"]


ms0 = 3e5*np.array([1, 0, 0 , 0 ,0, 0.001]).reshape((2,3))

#simulation building blocks
z1 = 8
# z1 = 20
z_end = int(z1+ magntic_steps + conductive_steps)
z_indexes = [z1, int(z1+ conductive_steps), z_end ]


save_loc = np.arange(z1-1, z_end+1)
save_loc[0] = 4
save_time= []


void_steps = 2*z1

exit_steps = int(time_cycal/ dt)



simulation_time = int(2.5*(peak_enter_time//dt + void_steps +4* (magntic_steps + conductive_steps)) + exit_steps)


magntic_eps = [Co_er, Fe_er, Ni_er ]
conductive_eps = [  Au_er, Pt_er]

magntic_sigma = [Co_sigma, Fe_sigma, Ni_sigma ]
conductive_sigma = [ Au_sigma, Pt_sigma]



epsr_list = [magntic_eps[0], conductive_eps[0]]
sigma_list = [magntic_sigma[0], conductive_sigma[0]]
            

data = sim.simulation(z_indexes, epsr_list , sigma_list, alpha[magntic_name[0] +"_" + conductive_name[0]] , ms0,
        B_opt, fie , FWHM,
        save_loc, save_time, False,
        dt, simulation_time, lamda=mean_wave_lenght, safety_start= sefe_start, time_interval= 5 )
      





H = data["H in locations"]
plt.plot(np.real(H[:,1,0]))
plt.show()



key_name = [
            "H times picture",
            "E times picture",
            "ms time picture",
            "time picture",
            "H in locations",
            "E in locations",
            "ms in locations",
            "locations",
            "time intervals",
            "matitral location", "e_r" ,"conductivity", "gilbert damping factor","initial magnetization",
            "max magnetic field" ,"polarization phase","pulse width (FWHM)","mean wave lenght",
            "systen status", "dt", "safety_start"                
            ]


