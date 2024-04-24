import numpy as np

def SplitStepProp(Beam, KZ, distance):
    """
    Propagate the beam over a certain distance using the split-step Fourier method.

    Parameters:
    - Beam: The complex beam field at the input plane.
    - KZ: The propagation constant matrix.
    - distance: The distance over which to propagate.

    Returns:
    - BeamProp: The beam field at the output plane after propagation.
    """
    print('before split size prop: ',Beam.shape)
    print('before split step prop: ',Beam[0,0])

    # Take the 2D Fourier transform of the input beam field
    BeamFFT = np.fft.fft2(Beam)
    
    # Apply the propagation factor in the frequency domain
    BeamK = BeamFFT * np.exp(-1j * (KZ * distance))
    
    # Inverse Fourier transform to get the propagated beam field
    BeamProp = np.fft.ifft2(BeamK)
    
    return BeamProp
