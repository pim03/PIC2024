import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def amplitude_phase_plot(Mode_Efield, ax):
    """
    Plots the amplitude and phase of a mode field.
    Amplitude determines the transparency and phase determines the color.
    
    Parameters:
    - Mode_Efield: The complex electric field mode to be visualized.
    - ax: The matplotlib axis on which to plot.
    """
    # Calculate amplitude and normalize
    Amp = np.abs(Mode_Efield)
    max_Amp = np.max(Amp)
    if max_Amp > 0:
        Amp /= max_Amp  # Normalize amplitude to 1 for transparency scaling
    else:
        Amp[:] = 0  # Set all to 0 if the max amplitude is 0 to avoid NaN issues

    # Calculate phase
    Ang = np.angle(Mode_Efield)
    
    # Set the color map for the phase and plot with transparency
    cmap = plt.get_cmap('hsv')
    norm = Normalize(vmin=-np.pi, vmax=np.pi)  # Normalize phase from -pi to pi
    img = ax.imshow(Ang, cmap=cmap, norm=norm)
    img.set_alpha(Amp)  # Set transparency based on amplitude
    
    # Set the plot background color to black
    ax.set_facecolor('black')
    
    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Show colorbar for phase
    plt.colorbar(img, ax=ax, orientation='vertical', label='Phase')

    return img
