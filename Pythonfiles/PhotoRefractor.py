import numpy as np

def PhotoRefractor(Beam, Nonlinearity=1):
    """
    Apply a photorefractive nonlinearity to the beam.
    
    Parameters:
    - Beam: The complex beam field.
    - Nonlinearity: The flag to indicate if the nonlinearity should be applied.
    
    Returns:
    - BeamProp: The beam after the nonlinearity has been applied.
    """
    
    if Nonlinearity == 1:
        # Photorefractive crystal nonlinearity
        NL = 2 * np.pi * np.abs(Beam)**2 / (1 + np.abs(Beam)**2)
        BeamProp = Beam * np.exp(-1j * NL)

        # Uncomment below for saturable absorber nonlinearity
        # Isat = 1e-4  # Saturation intensity
        # Tmin = 0.8  # Transmission for zero incident intensity
        # Tmax = 0.99  # Transmission at saturation intensity
        # NL = Tmin + np.abs(BeamProp)**2 / Isat * (Tmax - Tmin)  # Nonlinear transmission matrix
        # NL[NL > 1] = 1  # Impose that maximum transmission is unity
        # BeamProp = BeamProp * NL
    
    else:
        # If nonlinearity is not to be applied, return the input beam unaltered
        BeamProp = Beam
    
    return BeamProp
