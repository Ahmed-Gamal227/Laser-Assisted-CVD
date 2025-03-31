import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

def load_data(filename):
    """Load CSV data with blank line separators"""
    with open(filename) as f:
        blocks = []
        current_block = []
        for line in f:
            if line.strip() == '':
                if current_block:
                    blocks.append(np.loadtxt(current_block, delimiter=','))
                    current_block = []
            else:
                current_block.append(line)
        if current_block:  # Add last block if file doesn't end with blank line
            blocks.append(np.loadtxt(current_block, delimiter=','))
    return blocks

# Load simulation data
temperature_data = load_data('temperature.csv')
thickness_data = load_data('thickness.csv')

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Laser Assisted CVD Simulation Results', fontsize=16)

# Plot initial frames
im1 = ax1.imshow(temperature_data[0], cmap='hot', origin='lower', 
                extent=[0, 100, 0, 100], vmin=300, vmax=np.max(temperature_data[0]))
ax1.set_title('Temperature Field (K)')
ax1.set_xlabel('X Position (μm)')
ax1.set_ylabel('Y Position (μm)')

divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im1, cax=cax1)

im2 = ax2.imshow(thickness_data[0]*1e9, cmap='viridis', origin='lower', 
                extent=[0, 100, 0, 100], vmin=0, vmax=np.max(thickness_data[-1])*1e9)
ax2.set_title('Deposited Thickness (nm)')
ax2.set_xlabel('X Position (μm)')
ax2.set_ylabel('Y Position (μm)')

divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im2, cax=cax2, label='Thickness (nm)')

# Animation update function
def update(frame):
    im1.set_array(temperature_data[frame])
    im1.set_clim(vmin=300, vmax=np.max(temperature_data[frame]))
    
    im2.set_array(thickness_data[frame]*1e9)
    im2.set_clim(vmin=0, vmax=np.max(thickness_data[frame])*1e9)
    
    fig.suptitle(f'Laser Assisted CVD Simulation - Time Step {frame}/{len(temperature_data)-1}', fontsize=16)
    return im1, im2

# Create animation
ani = FuncAnimation(fig, update, frames=len(temperature_data), 
                   interval=200, blit=False)

plt.tight_layout()
plt.show()

# Optionally save the animation
# ani.save('lcvd_simulation.mp4', writer='ffmpeg', fps=5, dpi=200)import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

def load_data(filename):
    """Load CSV data with blank line separators"""
    with open(filename) as f:
        blocks = []
        current_block = []
        for line in f:
            if line.strip() == '':
                if current_block:
                    blocks.append(np.loadtxt(current_block, delimiter=','))
                    current_block = []
            else:
                current_block.append(line)
        if current_block:  # Add last block if file doesn't end with blank line
            blocks.append(np.loadtxt(current_block, delimiter=','))
    return blocks

# Load simulation data
temperature_data = load_data('temperature.csv')
thickness_data = load_data('thickness.csv')

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Laser Assisted CVD Simulation Results', fontsize=16)

# Plot initial frames
im1 = ax1.imshow(temperature_data[0], cmap='hot', origin='lower', 
                extent=[0, 100, 0, 100], vmin=300, vmax=np.max(temperature_data[0]))
ax1.set_title('Temperature Field (K)')
ax1.set_xlabel('X Position (μm)')
ax1.set_ylabel('Y Position (μm)')

divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im1, cax=cax1)

im2 = ax2.imshow(thickness_data[0]*1e9, cmap='viridis', origin='lower', 
                extent=[0, 100, 0, 100], vmin=0, vmax=np.max(thickness_data[-1])*1e9)
ax2.set_title('Deposited Thickness (nm)')
ax2.set_xlabel('X Position (μm)')
ax2.set_ylabel('Y Position (μm)')

divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im2, cax=cax2, label='Thickness (nm)')

# Animation update function
def update(frame):
    im1.set_array(temperature_data[frame])
    im1.set_clim(vmin=300, vmax=np.max(temperature_data[frame]))
    
    im2.set_array(thickness_data[frame]*1e9)
    im2.set_clim(vmin=0, vmax=np.max(thickness_data[frame])*1e9)
    
    fig.suptitle(f'Laser Assisted CVD Simulation - Time Step {frame}/{len(temperature_data)-1}', fontsize=16)
    return im1, im2

# Create animation
ani = FuncAnimation(fig, update, frames=len(temperature_data), 
                   interval=200, blit=False)

plt.tight_layout()
plt.show()

# Optionally save the animation
# ani.save('lcvd_simulation.mp4', writer='ffmpeg', fps=5, dpi=200)