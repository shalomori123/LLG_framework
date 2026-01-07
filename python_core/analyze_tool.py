import numpy as np
# from numba import jit
import matplotlib.pyplot as plt
# import simulation_cpu  as sim
from matplotlib.animation import FuncAnimation




def time_plot(E):
    E_dimntion = len(E)
    max_E = 1.1*np.max(np.real(E))
    def update(frame):
        ax.clear()
        ax.plot(np.real(E[frame,:,0]))
        # y_lim = max(max(max(E_t)))
        # ax.set_xlim(-10, 10)  # Adjust x-axis limits as needed
        ax.set_ylim(-max_E, max_E)  # Adjust y-axis limits as needed
        ax.set_title(f"Frame {frame}")


    fig, ax = plt.subplots()
    anim = FuncAnimation(fig, update, frames=E_dimntion, interval=100)
    plt.show()

