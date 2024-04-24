import numpy as np
import re
from SplitStepProp import SplitStepProp
from PhotoRefractor import PhotoRefractor
import matplotlib.pyplot as plt
from amplitude_phase_plot import amplitude_phase_plot

def elements(Beam, BeamBack, Arch_elements, Hologram, KZ, L1, ElementShow, PhScrInd, start, finish):
    if start < finish:
        counter_beam = PhScrInd -1
        add_counter = 1
        active_beam = Beam
    else:
        counter_beam = PhScrInd -1
        active_beam = BeamBack
        add_counter = -1
    # print('active beam: ',active_beam[0,0,0])
    for i in range(start, finish + add_counter, add_counter):
        element = Arch_elements[i]
        print('\n antes do if active_beam[:,:,counter_beam]: ',active_beam[0,0,counter_beam])
        print('iteração: ',i)
        print('elemento: ',element)
        # Simple Propagator
        if element.startswith('P'):
            print('entrei no simple prop na iteração: ',i)
            match = float(element[1:])
            if match: #ta menos 9.11j em vez de mais
                PropDist = match # group(1) captures the first set of parentheses
                print('counter_beam: ',counter_beam)
                print('dentro do if active_beam[:,0,counter_beam]: ',active_beam[0,0,counter_beam])
                active_beam[:, :, counter_beam] = SplitStepProp(active_beam[:, :, counter_beam], KZ, PropDist * 1e-3)
            else:
                print(f"No numeric distance found in element '{element}'")
            print('after simple propagator: ',active_beam[0,0,0])
        # Lens
        elif element.startswith('L'):
            print('entrei na lens na iteração: ',i)
            print('antes da lente active_beam[:,0,counter_beam]: ',active_beam[0,0,counter_beam])
            # print('L1: ',L1)
            active_beam[:, :, counter_beam] = active_beam[:, :, counter_beam] * L1
            print('counter beam: ',counter_beam)
            print('depois da lente active_beam[:,0,counter_beam]: ',active_beam[0,0,counter_beam])


        # Phasemask
        elif element.startswith('H') and (counter_beam != 1 or start < finish): # this makes sure that for bacwards propagation we dont generate the last part
            phase_mask_index = int(re.search(r'\d+', element).group()) -1 
            active_beam[:, :, counter_beam + add_counter] = active_beam[:, :, counter_beam] * np.exp(-1j * Hologram[:, :, phase_mask_index])
            counter_beam += add_counter

        elif element.startswith('N'):
            active_beam[:, :, counter_beam] = PhotoRefractor(active_beam[:, :, counter_beam])

        if ElementShow:
            plt.figure()
            ax = plt.gca()
            amplitude_phase_plot(active_beam[:, :, counter_beam], ax)
            plt.title(f'Iteration: {i}, Output after: {add_counter}, Element: {element}, Phase Screen Num: {counter_beam}')
            plt.show()

    return active_beam
