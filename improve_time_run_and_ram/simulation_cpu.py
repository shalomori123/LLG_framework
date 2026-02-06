import numpy as np
from tqdm import tqdm
# import kernel_function

eps=8.854e-12 #[F/m]
miu=4*np.pi*1e-7 #[H/m]
c=1/(eps*miu)**0.5
heta = (miu/eps)**0.5

q = 1.60217646e-19    # Elementary charge [Coulombs]
miu = 4 * np.pi * 1e-7    # Magnetic permeability [H/m]
g = 2    # Landau factor
me = 9.1093821545e-31    # Electron mass [kg]
gma_factor = 1
gma = gma_factor * g * q / (2 * me)    # Gyromagnetic ratio [rad/sec*T]


def pulse_gen(max_feald: float, fea : float, FWHM: float, dt : float ,lamda = 800e-9 , teta = 0, start_sefty = 4 ):
    """
    input-
    :param max_feald: max value of the fileds [A/m]
    :param fea: phase defernt betwin E_x and E_y
    :param FWHM: time from helf power to helf of optical autoroltion [sec].
    :patm dt: time step size.
    :param lamda: wave len in void to the AC amplufind 
    :param teta: angle of the central axis of the polarization [rad]
    
    :return: function out_pulse(n) 
    """
    # move to time w
    w = 2* np.pi * c/lamda
    teta -= np.pi/4
    turn_m= np.array([np.cos(teta) , -np.sin(teta), np.sin(teta) ,np.cos(teta) ])
    turn_m = turn_m.reshape((2,2))
    pulse = max_feald * turn_m @ np.array([1 ,np.exp(-fea*1j)]).reshape(2)
    sigma = (0.5/(np.log(2)** 0.5)) * FWHM
    peak_location = start_sefty*sigma

    def out_pulse(n: int) -> np.array:
        """
        input-
        :patam n: index of step in time
        :return: 2D comlmx np.array of the inject pulse (E_x E_y) in this time step 
        """
        r = pulse  * np.exp(- (n*dt - peak_location)**2/( 2*sigma**2))* np.exp(1j*w*n*dt)
        return r
    return out_pulse 


def LLG_step(M: np.array, H: np.array, dt: float, alpha: float) -> np.array:
    """"""


    M0 = np.linalg.norm(M, axis=1, keepdims=True)

    
    gma_LL=gma/((1+alpha**2))
    LL_lamda=gma*alpha/(1+alpha**2)
    
   # Compute LLG terms
    an = -gma_LL * miu * np.cross(M, H, axis=1) - (LL_lamda * miu / M0) * np.cross(M, np.cross(M, H, axis=1), axis=1)
    
    bn = -gma_LL * miu * np.cross(M + (dt / 2) * an, H, axis=1) - (LL_lamda * miu / M0) * np.cross(M + (dt / 2) * an, np.cross(M + (dt / 2) * an, H, axis=1), axis=1)
    
    cn = -gma_LL * miu * np.cross(M + (dt / 2) * bn, H, axis=1) - (LL_lamda * miu / M0) * np.cross(M + (dt / 2) * bn, np.cross(M + (dt / 2) * bn, H, axis=1), axis=1)
    
    dn = -gma_LL * miu * np.cross(M + dt * cn, H, axis=1) - (LL_lamda * miu / M0) * np.cross(M + dt * cn, np.cross(M + dt * cn, H, axis=1), axis=1)
    
    new_M = M + (dt/6)*(an+2*bn+2*cn+dn)
    return new_M


def simulation(z_ind: list[int], e_r : list[float], sigma: list[float], alpha: float, m0: np.array,
               max_E: float, fea : float, FWHM: float,
               save_location, save_time, close_system: bool,
               dt: float, simulation_time_steps: int, lamda = 800e-9 , teta = 0, safety_start = 4, time_interval = 1):
    """
    inputs-
    :param z1: start index of material
    :param z2: end index of material
    :param e_r:  relative permittivity
    :param sigma: condactivity of matiral [Oham/m]
    :param alpha: deeiha param of moments
    :param m0: initial BIG MAGNETIZATION (sum of moments) configuration in [x,y,z] -> vector size 3
    :param max_E: max value of the fileds [A/m]
    :param fea: phase defernt betwin E_x and E_y
    :param FWHM: the whith of the gausian , time when the pulse is greater than half max power [sec].
    :patm save_location: interbel of index in space to save all value in the for all time.
    :patm save_time: interbel of index in time to save all value at the intaer space.
    :param close_system: wether will use the reflection of M on H on E or not.
    :param dt: time step size [s]
    :param simulation_time_steps: numer for time steps to do [int]
    :param lamda: amplifid wave lenght [m]
    :param teta: angle of the central axis of the polarization [rad]
    
    :return: two list in len 3
        first list - save parmaters in all space in the save_time stepes
        secnd list - saved paarmates in all time in the seve_location ind on space
        parm are: [H, E, msi ] ( msi given only in the relvent space if axis) 
                    all of them in the shep of [time_step]* [space_step]* [3 demntion]
    """


    # system defntion 
    dz = 2*c *dt
    lenz= z_ind[-1] + min(int(lamda/dz), int(0.5*(z_ind[-1]-z_ind[0])))
 

    # matrial proprty
    const_for_les_E = np.ones(lenz-1)
    const_for_H = 0.5*np.ones(lenz-1)
    for i in range(len(e_r)):
        eps_r = np.ones(z_ind[i+1] - z_ind[i]) 
        eps_r += (e_r[i]-1) 
        eaf = dt/2/eps* (sigma[i]) /eps_r
        const_for_les_E[z_ind[i]:z_ind[i+1]] = ((1-eaf)/(1+eaf))
        const_for_H[z_ind[i]:z_ind[i+1]] = (0.5/eps_r/(1+eaf))
    

    
    # plase for run the symolation
    pulse = pulse_gen(max_E, fea, FWHM, dt, lamda, teta, safety_start)

    # defing the worke spaec for the simltion
    E_next= np.zeros((lenz, 3)).astype(complex)
    E_current= np.zeros((lenz, 3)).astype(complex)
    right_boundry = [np.zeros((3)), np.zeros((3))]
    left_boundry= [np.zeros((3)), np.zeros((3))]
    H_next= np.zeros((lenz, 3)).astype(complex)
    H_current= np.zeros((lenz, 3)).astype(complex)
    
    ms_i_current= np.zeros((z_ind[-1] -z_ind[0], 3))
    # intial the beging magntiztion
    for i in range(len(e_r)):
        ms_us = m0[i,:]
        st_ind = z_ind[i]- z_ind[0]
        end_ind = z_ind[i+1] - z_ind[0]
        for k in range(3):
            ms_i_current[st_ind:end_ind,k]= ms_us[k]/(end_ind-st_ind)


    # defing rutenting array:
    # for space indxes
    time_save = int(simulation_time_steps/time_interval + 1)
    E_space_return = np.zeros((time_save, len(save_location), 2)).astype(complex)
    H_space_return = np.zeros((time_save, len(save_location), 2)).astype(complex)
    msi_space_return = np.zeros((time_save, z_ind[-1]-z_ind[0],3)).astype(complex)

    # for time indxes
    E_time_return = np.zeros((len(save_time), lenz, 2)).astype(complex)
    H_time_return = np.zeros((len(save_time), lenz, 2)).astype(complex)
    msi_time_return = np.zeros((len(save_time), z_ind[-1]-z_ind[0] , 3)).astype(complex)
    index_to_save_in =0
    
    # propegate

    for n in tqdm(range(0, simulation_time_steps), desc= "simulation pross: "):
        E_current[ 3,0:2] += pulse(n) # inject the input signal

        # set bondry condition for not replction
        r = right_boundry.pop(0)
        l = left_boundry.pop(0)

        E_current[ -1,:] = r
        E_current[0,:] = l
    
        # H eqution
        H_next[1:, 0]= H_current[1:, 0]+ 0.5*(E_current[1:,1] - E_current[:-1, 1] )        
        H_next[1:, 1]= H_current[1:, 1]- 0.5*(E_current[1:, 0] - E_current[:-1, 0] )
        # LLG and cang in the ms_i        
        ms_i_next = LLG_step(ms_i_current, np.real(H_current[z_ind[0]:z_ind[-1],:]), dt, alpha)
        if close_system: # intagrtion line for close sistem B = H + M
            H_next[ z_ind[0]:z_ind[-1], :] += ms_i_current - ms_i_next
        # E eqution
        E_next[ :-1, 1] = const_for_les_E*E_current[:-1,1] + const_for_H* (H_next[1:,0] - H_next[ :-1, 0])
        E_next[ :-1, 0] = const_for_les_E*E_current[:-1, 0] - const_for_H* (H_next[1:, 1] - H_next[ :-1, 1])

     
        # get redy to next time step
        right_boundry.append(np.copy(E_current[ -2,:]))
        left_boundry.append(np.copy(E_current[ 1,:]))
        E_current = np.copy(E_next)
        H_current = np.copy(H_next)
        ms_i_current = np.copy(ms_i_next)

        # save to return
        # save picture of parameters in time
        if n in save_time:
            H_time_return[index_to_save_in, :,:] = H_current[:,:2]
            E_time_return[index_to_save_in, :,:] = E_current[:,:2]
            msi_time_return[index_to_save_in, :,:] = ms_i_current

            index_to_save_in += 1

        # save the parameters in speceific location in all times
        if not n % time_interval :
            p =int(n/time_interval)
            H_space_return[p, :, :] = H_current[save_location, :2]
            E_space_return[p, :, :] = E_current[save_location, :2]
            msi_space_return[p, :, :] = ms_i_current
    # time_r = H_time_return, E_time_return, msi_time_return
    # spece_r = H_space_return, E_space_return, msi_space_return
    all_data = {
                "H times picture": H_time_return,
                "E times picture": E_time_return,
                "ms time picture": msi_time_return,
                "time picture": save_time,
                "H in locations": H_space_return,
                "E in locations": E_space_return,
                "ms in locations": msi_space_return,
                "locations": save_location,
                "time intervals": time_interval,
                "matitral location": z_ind, "e_r" : e_r ,"conductivity": sigma, "gilbert damping factor": alpha,"initial magnetization": m0,
                "max magnetic field" :max_E ,"polarization phase": fea,"pulse width (FWHM)": FWHM, "mean wave lenght" : lamda,
                "systen status": close_system, "dt": dt, "safety_start": safety_start
                }       
    return  all_data

 

if __name__ == '__main__':
    print("i do this")
    lamda = 8e-7 
    time_cycal= lamda/c
    bwhp = 4*time_cycal 
    max_E = 3e5
    fea = np.pi/2


    
    # e_r = 2
    # sigma = 0
    # alpha = 0
    # ms0= np.array([0, 0 , 1])*1e4
    # len_matiral = 8*lamda
    

    # time_cut = 80
    # dt = time_cycal/time_cut
    # dz = 2*c *dt


    # # spaece and time defntion
    # z1 = int(5*bwhp/dt) +3
    # z2 = z1 + int(len_matiral/dz) 
    # lenz= z2 + time_cut//2
    # simulation_time_steps = int(4*z1 + 5*(e_r**0.5)*(z2- z1)+time_cut)
    
    # save_location = [4, z1, z2-1, z2+2]
    # save_time  = [2*bwhp//dt, bwhp//dt +2*z1 ]



    # time_res, specr_res = simulation(z1, z2, e_r , sigma, alpha, ms0,
    #            max_E, fea , bwhp,
    #            save_location, save_time, False,
    #            dt, simulation_time_steps, lamda )


    # x= 1
