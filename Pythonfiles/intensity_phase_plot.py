import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

def intensity_phase_plot(Mode_Efield):
    """
    Displays the transverse electric field's phase with a colormap and
    the intensity overlayed as alpha (transparency) values, along with a
    black background.
    """
    Int = np.abs(Mode_Efield)**2
    Int = Int / np.max(Int)
    Ang = np.angle(Mode_Efield)
    
    # Create HSV image
    hsv_image = np.zeros((*Ang.shape, 3))  # Initialize HSV image with the same size as Ang
    hsv_image[..., 0] = (Ang + np.pi) / (2 * np.pi)  # Hue from the phase
    hsv_image[..., 1] = 1  # Full saturation
    hsv_image[..., 2] = Int  # Value from the intensity
    
    # Convert HSV to RGB
    rgb_image = hsv_to_rgb(hsv_image)
    
    # Display the image
    plt.imshow(rgb_image)
    plt.gca().set_facecolor((0, 0, 0))  # Set background to black
    plt.axis('off')  # Get rid of axes
    plt.show()

