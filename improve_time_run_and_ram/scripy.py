import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Example data (sine wave)
time_points = np.linspace(0, 10, 100)
data = np.sin(time_points)

# Initialize the plot
fig, ax = plt.subplots()

# Plotting function for updating the plot frame by frame
def update(frame):
    ax.clear()
    ax.plot(time_points[:frame], data[:frame])
    ax.set_xlim(0, 10)  # Adjust x-axis limits as needed
    ax.set_ylim(-1, 1)  # Adjust y-axis limits as needed
    ax.set_title(f"Frame {frame}")

# Create the animation
anim = FuncAnimation(fig, update, frames=len(time_points), interval=100)

# Display the plot
plt.show()