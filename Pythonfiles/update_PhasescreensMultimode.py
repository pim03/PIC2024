# import numpy as np

# def update_PhasescreensMultimode(Beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue):
#     """
#     Update phase screens using the current Beam and BeamBack values.
#     Arguments:
#     - Beam: Current forward propagating beam array
#     - BeamBack: Current backward propagating beam array (conjugated initially)
#     - Hologram: Current state of the hologram
#     - PhScrInd: Phase screen index being updated
#     - nx, ny: Dimensions of the phase screen
#     - Nmodes: Number of modes
#     - CoarsePixelation: Boolean to determine if the SLM grid should be made coarser
#     - NPixelsCoarse: Number of pixels in the coarser grid
#     - MaxPhaseValue: Maximum phase value for discretization
#     """
#     DeltaHologram = np.zeros((nx, ny, Nmodes), dtype=complex)

#     for jmodes in range(Nmodes):
#         # Overlap = BeamBack[:, :, PhScrInd, jmodes] * Beam[:, :, PhScrInd, jmodes] * np.exp(1j * Hologram[:, :, PhScrInd])
#         Overlap = BeamBack[:, :, PhScrInd, jmodes] * Beam[:, :, PhScrInd, jmodes] * np.exp(1j * Hologram[:, :, PhScrInd, jmodes])
#         AvePhase = np.angle(np.mean(Overlap))
#         DeltaHologram[:, :, jmodes] += Overlap * np.exp(-1j * AvePhase)

#     DeltaHoloSum = -np.angle(np.sum(DeltaHologram, axis=2))
#     DeltaHoloSum = DeltaHoloSum[:, :, np.newaxis]

#     Hologram[:, :, PhScrInd] += DeltaHoloSum
#     Hologram[:, :, PhScrInd] = np.mod(Hologram[:, :, PhScrInd], 2 * np.pi)
#     # Hologram[:, :, PhScrInd] = np.mod(np.real(Hologram[:, :, PhScrInd]), 2 * np.pi)


#     if CoarsePixelation:
#         NCellsTiles = nx // NPixelsCoarse
#         for x in range(0, nx, NCellsTiles):
#             for y in range(0, ny, NCellsTiles):
#                 tile = Hologram[x:x + NCellsTiles, y:y + NCellsTiles, PhScrInd]
#                 Hologram[x:x + NCellsTiles, y:y + NCellsTiles, PhScrInd] = np.mean(tile)

#     max_phase = np.max(Hologram[:, :, PhScrInd])
#     if max_phase > 0:
#         Hologram[:, :, PhScrInd] = np.floor(Hologram[:, :, PhScrInd] / max_phase * MaxPhaseValue)
#         Hologram[:, :, PhScrInd] = (Hologram[:, :, PhScrInd] / MaxPhaseValue) * 2 * np.pi

#     return Hologram


import numpy as np

def update_PhasescreensMultimode(Beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue):
    """
    Update phase screens using the current Beam and BeamBack values.
    Arguments:
    - Beam: Current forward propagating beam array
    - BeamBack: Current backward propagating beam array
    - Hologram: Current state of the hologram
    - PhScrInd: Phase screen index being updated
    - nx, ny: Dimensions of the phase screen
    - Nmodes: Number of modes
    - CoarsePixelation: Boolean to determine if the SLM grid should be made coarser
    - NPixelsCoarse: Number of pixels in the coarser grid
    - MaxPhaseValue: Maximum phase value for discretization
    """
    oldhologram = Hologram
    DeltaHologram = np.zeros((nx, ny, Nmodes), dtype=complex)
    # print(f'Beam: ',Beam[0,0,PhScrInd,:])
    # print(f'BeamBack: ',BeamBack[0,0,PhScrInd,:])
    # print(f'Hologram: ',Hologram)
    for jmodes in range(Nmodes):
        Overlap = BeamBack[:, :, PhScrInd, jmodes] * Beam[:, :, PhScrInd, jmodes] * np.exp(1j * Hologram[:, :, PhScrInd, jmodes])
        # print(f'overlap: ',Overlap)
        AvePhase = np.angle(np.mean(Overlap))
        # print(f'avephase: ',AvePhase)
        DeltaHologram[:, :, jmodes] += Overlap * np.exp(-1j * AvePhase)
        # print(f'deltaHologram: ',DeltaHologram)

    DeltaHoloSum = -np.angle(np.sum(DeltaHologram, axis=2))

    DeltaHoloSum = DeltaHoloSum[:, :, np.newaxis]
    # print('DeltaHoloSum: ', DeltaHoloSum)
    Hologram[:, :, PhScrInd] += DeltaHoloSum
    # Hologram[:, :, PhScrInd] = np.mod(Hologram[:, :, PhScrInd], 2 * np.pi)

    # Ensure the Hologram is in the correct data type
    # Hologram[:, :, PhScrInd] = np.mod(Hologram[:, :, PhScrInd].astype(np.float64), 2 * np.pi)
    Hologram[:, :, PhScrInd] = np.mod(np.real(Hologram[:, :, PhScrInd]), 2 * np.pi)

    # print(f'Hologram: ',Hologram)

    if CoarsePixelation:
        NCellsTiles = nx // NPixelsCoarse
        for x in range(0, nx, NCellsTiles):
            for y in range(0, ny, NCellsTiles):
                tile = Hologram[x:x + NCellsTiles, y:y + NCellsTiles, PhScrInd]
                Hologram[x:x + NCellsTiles, y:y + NCellsTiles, PhScrInd] = np.mean(tile)

    max_phase = np.max(Hologram[:, :, PhScrInd])
    
    if max_phase > 0:
        # Hologram[:, :, PhScrInd] = np.floor(Hologram[:, :, PhScrInd] / max_phase * MaxPhaseValue)
        # Assuming Hologram should be real because it represents phase information
        Hologram[:, :, PhScrInd] = np.floor(np.real(Hologram[:, :, PhScrInd] / max_phase * MaxPhaseValue))
        Hologram[:, :, PhScrInd] = Hologram[:, :, PhScrInd] / MaxPhaseValue * 2 * np.pi

    # print(f'Hologram: ',Hologram)
    print('old and new hologram are the same: ', np.allclose(oldhologram, Hologram))
    return Hologram
