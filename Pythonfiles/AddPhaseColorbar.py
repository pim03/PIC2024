import matplotlib.pyplot as plt
import numpy as np

def add_colorbar(mappable, fontsize=26, fontname='Times New Roman'):
    """
    Adds a colorbar to a plot based on the given mappable object with customized ticks, labels, and font properties.

    Parameters:
    - mappable: The mappable object associated with a plotted dataset, typically returned from a plot command.
    - fontsize: Font size of the colorbar ticks.
    - fontname: Font name of the colorbar ticks.
    """
    cbar = plt.colorbar(mappable)
    cbar.set_ticks([-np.pi, 0, np.pi])
    cbar.set_ticklabels(['$-\pi$', '$0$', '$\pi$'])

    # Set the font size for the colorbar ticks. For font family, you may need to adjust depending on matplotlib's support.
    cbar.ax.tick_params(labelsize=fontsize)
    for label in cbar.ax.get_yticklabels():
        label.set_family(fontname)
