import numpy as np
from scipy.special import genlaguerre

def GenModesLG(ModeTypes, waist, Rad, Angle):
    """
    This function calculates Laguerre-Gauss modes.
    The created mode indices should be listed in the matrix "ModeTypes",
    where the first column should include the OAM values and the second
    column should include the corresponding p-values.
    - waist should be the wanted beam waist used in calculation
    - Rad should be the wanted 2D radial coordinate grid
    - Angle contains the corresponding azimuthal coordinates

    Output is a 3D matrix containing all of the input mode structures stacked
    together in the 3rd dimension of the matrix
    """
    ModeNum = ModeTypes.shape[0]
    GridDim = Rad.shape
    LGout = np.zeros((GridDim[0], GridDim[1], ModeNum), dtype=complex)
    for ModeInd in range(ModeNum):  # Loop through all of the listed modes
        # Extract OAM and p-value
        l = ModeTypes[ModeInd, 0]
        p = ModeTypes[ModeInd, 1]
        # Calculate the Mode using the Laguerre polynomial function:
        LGTemp = ((np.sqrt(2) * Rad / waist) ** (abs(l))) * np.exp(-(Rad ** 2) / waist ** 2) \
                 * np.polyval(genlaguerre(p, abs(l)), 2 * (Rad ** 2) / waist ** 2) \
                 * np.exp(1j * l * Angle)
        # Normalize the field:
        NormA = LGTemp * np.conj(LGTemp)
        SumOfA = np.sum(NormA)
        A = np.sqrt(np.real(SumOfA))
        # Store the value to the output:
        LGout[:, :, ModeInd] = LGTemp / A
    return LGout

