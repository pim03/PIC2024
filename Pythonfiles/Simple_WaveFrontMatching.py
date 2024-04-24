import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import Normalize

from GenModesLG import GenModesLG
from amplitude_phase_plot import amplitude_phase_plot
from elements import elements
from update_PhasescreensMultimode import update_PhasescreensMultimode


########################## Constants and Parameters ############################

# Nonlinear layer between Hologram planes
Nonlinearity = 0 #1 = add nonlinear layers in propagation (SplitStepProp.m)

# Regarding optimization:
ConvStop = 0 #1 = stop iterations if convergence criterion is met
# Convergence accuracy (number of decimals):
convLim = 4
# Number of previous iterations to consider when checking for convergence:
convValue = 3
# Maximum number of iterations:
IterMax = 200
    
# SLM specifications:
PixelSize = 8e-6 # width = height(m) --> it sets also the x,y grid step for the spatial modes
#PhaseScreenNum = 1 # Number of Phase screens used in optimization
MaxPhaseValue = 255 # Number of discrete phase shift values in the 2pi range
CoarsePixelation = 0 #1 = make the SLM grid coarser
NPixelsCoarse = 256

# Setup parameters:
lambda_ = 1550e-9 # Wavelength (m)
w_in = 0.2e-3 # Input beam radius (m) 3-state superpos.
w_out = 0.2e-3 # Output beam radius (m) 3-state superpos.
f1 = 31 #Lens focal distance
PropagationHologram = 1e-2 #Distance between holograms

# Other parameters
# Simulation window size (bigger helps with unwanted "computational box reflections"
# in split step propagation):
WidthX = 4e-3 # SLM window size or arbitrary size (m)
WidthY = 4e-3 # SLM window size or arbitrary size (m)

#System Architecture   # Px propagation over x mm  #Hx beam goes through mask x  # L make beam go through lens with previously defined focal lenght # separated by space character
#Arch = [strcat('P', string(f1*1000),' L P', string(f1*1000-PropagationHologram/2), ' H1 P', string(PropagationHologram), ' H2 P', string(f1*1000-PropagationHologram/2), ' L P400')] # always in [mm]
Arch = f'P{f1} L P{f1} H1 P{f1} L P{f1}' # Definition for 1 Hologram
#Arch = 'P300 L P300 H1 P300 L P300'

# Splitting the string into a list of elements
Arch_elements = Arch.split()
print('Current Architecture: ', Arch_elements, '\n')

# Finding the positions of elements that start with 'H' and 'L'
hologram_position = [i for i, element in enumerate(Arch_elements) if element.startswith('H')]
lens_position = [i for i, element in enumerate(Arch_elements) if element.startswith('L')]

# Counting the number of occurrences of 'H' in the string
PhaseScreenNum = Arch.count('H')


######################## Flags controlling optimization:#######################

# When are the phase screens updated: (Both can be used at the same time)
# (This can slightly affect which phase screen removes phase singularities
# [first or last usually])
# Update when propagating backwards:
BackUpdate = True
# Update when propagating forwards:
ForwardUpdate = True

# Plot beam forward and backward at each element
ElementShow = False

########################## Input and output modes:##############################



# Input: (Only LG modes here)
# First column is OAM (azimuthal) index, second is radial index
# Example: 
# InputModes = np.array([[-1, 0]])  # Here the mode is an OAM=-1 mode without radial structure

InputModes = np.array([[0, 0]])  # if size(InputModes, 2) > 2 --> superposition of states

InputSuperposCoeffs = np.array([1, 1])  # Coefficients of the superposition --> [[a*superpos1 modeA + b*superpos1 modeB];[c*superpos2 modeA + d*superpos2 modeB]]

# Normalization
NsuperposInput = InputModes.shape[1] // 2
InputSuperposCoeffs = InputSuperposCoeffs * (1 / np.sqrt(NsuperposInput)) # normalization

# Output mode
OutputModes = np.array([[1, 0]])

Nmodes = InputModes.shape[0]

################### Calculate simulation parameters from inputs ################

# Grid:
nx = round(WidthX / PixelSize)  # Amount of pixels (x-dir)
ny = round(WidthY / PixelSize)  # Amount of pixels (y-dir)
x = np.linspace(-WidthX / 2, WidthX / 2, nx, endpoint=False)
y = np.linspace(-WidthY / 2, WidthY / 2, ny, endpoint=False)
X, Y = np.meshgrid(x, y)  # Cartesian grid

# Grid in cylindrical coordinates:
Rad = np.sqrt(X**2 + Y**2)
Angle = np.angle(X + 1j * Y) + np.pi  # Matrix with all angles starting left-center

# K-space grid used in split step propagation:
kx = (np.mod(np.arange(nx) + nx // 2, nx) - nx // 2) * 2 * np.pi / WidthX
ky = (np.mod(np.arange(ny) + ny // 2, ny) - ny // 2) * 2 * np.pi / WidthY
KX, KY = np.meshgrid(kx, ky)
# Matrix for k-space propagation direction components:
KZ = np.sqrt((2 * np.pi / lambda_)**2 - (KX**2 + KY**2))

# Fourier Lens definition
L1 = np.exp(1j * 2 * np.pi / lambda_ * (f1 - np.sqrt(X**2 + Y**2 + f1**2)))
# Visualization
plt.figure()
plt.imshow(np.angle(L1), extent=(x.min(), x.max(), y.min(), y.max()))
plt.colorbar()  # Show color bar for reference
plt.title('Phase of Fourier Lens')
plt.xlabel('X position')
plt.ylabel('Y position')
plt.show()

########################### Calculate Used modes ###############################

# Calculate input and output modes
ModesIn = np.zeros((nx, ny, Nmodes), dtype=complex)
ModesOut = np.zeros((nx, ny, Nmodes), dtype=complex)

print("Calculating modes \n")

# Plotting the modes
fig, axs = plt.subplots(Nmodes, 2, figsize=(10, 5 * Nmodes))  # Adjust the figure size as needed

for jmodes in range(Nmodes):
    for jsuperposInput in range(Nmodes):
        mode_indices = InputModes[jmodes, jsuperposInput * 2:(jsuperposInput + 1) * 2]
        ModesIn[:, :, jmodes] += InputSuperposCoeffs[jmodes] * GenModesLG(np.array([mode_indices]), w_in, Rad, Angle)[:, :, 0]
    # Define output modes at last phase mask plane, after going through non linear layer    
    ModesOut[:, :, jmodes] = GenModesLG(np.array([OutputModes[jmodes]]), w_out, Rad, Angle)[:, :, 0]

    # Input modes
    ax = axs[jmodes, 0] if Nmodes > 1 else axs[0]
    img = amplitude_phase_plot(ModesIn[:, :, jmodes], ax)
    ax.set_title('Input Mode')
    # fig.colorbar(img, ax=ax, orientation='vertical', label='Phase')

    # Output modes
    ax = axs[jmodes, 1] if Nmodes > 1 else axs[1]
    img = amplitude_phase_plot(ModesOut[:, :, jmodes], ax)
    ax.set_title('Output Mode')
    # fig.colorbar(img, ax=ax, orientation='vertical', label='Phase')

fig.tight_layout()  # Adjust layout to prevent overlap
plt.show()

############# Propagate the input and output before starting optimization #####

print("Initial propagation \n \n")
# Assuming the following variables are defined earlier in your Python script:
# nx, ny, Nmodes, ModesIn, ModesOut, PhaseScreenNum, Arch_elements, KZ, L1, ElementShow
# Also assuming that hologram_position is a list of indices where hologram elements are located

# Initialize Beam and BeamBack for each mode
Beam = np.zeros((nx, nx, PhaseScreenNum + 1, Nmodes), dtype=complex)
BeamBack = np.zeros((nx, nx, PhaseScreenNum, Nmodes), dtype=complex)
Hologram = np.zeros((nx, nx, PhaseScreenNum, Nmodes), dtype=complex)

print('\n \n Before Initial Propagation')
print(f'Beam: ',Beam[0,0,0,:])
print(f'BeamBack: ',BeamBack[0,0,0,:])

# print('Modes In: ',ModesIn[0,0,:])
# Propagating the input forward and backward:
for jmodes in range(Nmodes):
    # Copy the input mode to the first "slice" of the Beam array
    Beam[:, :, 0, jmodes] = ModesIn[:, :, jmodes]

    # Forward propagation
    print(f'beam shape: ', Beam.shape)
    Beam[:, :, :, jmodes] = elements(
        Beam=Beam[:, :, :, jmodes],
        BeamBack=None,  # Not used in forward propagation
        Arch_elements=Arch_elements,
        Hologram=Hologram[:, :, :, jmodes],
        KZ=KZ,
        L1=L1,
        ElementShow=ElementShow,
        PhScrInd=1,  # Starting at the first phase screen index for forward propagation
        start=0,
        finish=hologram_position[-1] + 1
    )
    print(f'Beam after elements: ', Beam[0,0,0,:])
    # Prepare the output mode for backward propagation with the negative phase
    BeamBack[:, :, PhaseScreenNum-1, jmodes] = np.abs(ModesOut[:, :, jmodes]) * np.exp(-1j * np.angle(ModesOut[:, :, jmodes]))
    # Backward propagation
    BeamBack[:, :, :, jmodes] = elements(
        Beam=None,  # Not used in backward propagation
        BeamBack=BeamBack[:, :, :, jmodes],
        Arch_elements=Arch_elements,
        Hologram=Hologram[:, :, :, jmodes],
        KZ=KZ,
        L1=L1,
        ElementShow=ElementShow,
        PhScrInd=PhaseScreenNum,  # Starting at the last phase screen index for backward propagation
        start=len(Arch_elements) - 1,
        finish=hologram_position[0] + 1
    )

print('\n \n After Initial Propagation')
print(f'Beam: ',Beam[0,0,0,:])
print(f'BeamBack: ',BeamBack[0,0,0,:])

# If you need to plot the initial amplitude and phase of the beam, you can do so with the amplitude_phase_plot function.

################### Creating a Button to stop optimization ####################
# def stop_optimization():
#     print("Optimization stopped")
#     window.destroy()  # This will close the window

# # Create a new window
# window = tk.Tk()
# window.title("Stop Optimization")
# window.geometry("100x50+300+300")  # Size and position: Width x Height + X_offset + Y_offset

# # Create a stop button
# stop_button = tk.Button(window, text="Stop WFM", command=stop_optimization)
# stop_button.pack(expand=True, fill=tk.BOTH)  # Expand the button to fill the window

# # Start the application
# window.mainloop()


# Define the main application window
root = tk.Tk()
root.title("WaveFront Matching Control")

# Add a button to control the optimization process
stop_pressed = False

def on_stop():
    global stop_pressed
    stop_pressed = True
    messagebox.showinfo("Process Info", "Optimization stopped by user.")

stop_button = tk.Button(root, text="Stop WFM", command=on_stop)
stop_button.pack()

# Set up plotting within the Tkinter GUI
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack()


########################### WaveFront Matching #################################

print('Wave Front Matching \n \n')

# Initialize arrays for storing results
resultoverlap = np.nan * np.ones((Nmodes, IterMax))
lines = []

# Initial propagation (Forward and backward setup)
for jmodes in range(Nmodes):
    Beam[:, :, 0, jmodes] = ModesIn[:, :, jmodes]
    BeamBack[:, :, -1, jmodes] = np.conjugate(ModesOut[:, :, jmodes])

IterationCount = 0

# while IterationCount < IterMax and not stop_pressed:
#     IterationCount += 1
#     print(f"Current iteration: {IterationCount}")

#     # Propagate and update in the forward direction
#     for PhScrInd in range(PhaseScreenNum):
#         if ForwardUpdate:
#             Hologram = update_PhasescreensMultimode(Beam, BeamBack, Hologram, 0, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

#         if PhScrInd == PhaseScreenNum:
#             Beam[:,:,:,jmodes] = elements(Beam[:,:,:,jmodes], BeamBack[:,:,:,jmodes], Arch_elements, Hologram[:,:,:,jmodes], KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[-1])# , Nonlinearity)
#         else:
#             Beam[:,:,:, jmodes] = elements(Beam[:,:,:,jmodes], BeamBack[:,:,:,jmodes], Arch_elements, Hologram[:,:,:,jmodes], KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[PhScrInd+1]) #, Nonlinearity)

#     # Propagate and update in the backward direction
#     for PhScrInd in range(PhaseScreenNum-1, -1, -1):
#         if BackUpdate:
#             # print(f"Updating phase screen {PhScrInd} in the backward direction")
#             Hologram = update_PhasescreensMultimode(Beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

#         for jmodes in range(Nmodes):
#             if PhScrInd > 1:
#                 BeamBack[:, :, :, jmodes] = elements(Beam, BeamBack, Arch_elements, Hologram, KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[PhScrInd-1], Nonlinearity)

    # Beam[:, :, :, jmodes] = elements(
    #     Beam=Beam[:, :, :, jmodes],
    #     BeamBack=None,  # Not used in forward propagation
    #     Arch_elements=Arch_elements,
    #     Hologram=Hologram[:, :, :, jmodes],
    #     KZ=KZ,
    #     L1=L1,
    #     ElementShow=ElementShow,
    #     PhScrInd=1,  # Starting at the first phase screen index for forward propagation
    #     start=1,
    #     finish=hologram_position[-1] + 1
    # )

while IterationCount < IterMax and not stop_pressed:


    IterationCount += 1
    print(f"Current iteration: {IterationCount}")

    # Forward and backward propagation
    for direction in ['forward', 'backward']:
        if direction == 'forward':
            current_beam = Beam
            start_index, end_index, step = 0, PhaseScreenNum, 1
        else:
            current_beam = BeamBack
            start_index, end_index, step = PhaseScreenNum-1, -1, -1

        for PhScrInd in range(start_index, end_index, step):
            for jmodes in range(Nmodes):
                current_beam[:, :, :, jmodes] = elements(
                    current_beam[:, :, :, jmodes], 
                    current_beam[:, :, :, jmodes], 
                    Arch_elements, 
                    Hologram[:, :, :, jmodes], 
                    KZ, 
                    L1, 
                    ElementShow, 
                    PhScrInd, 
                    hologram_position[PhScrInd], 
                    lens_position[-1 if PhScrInd == PhaseScreenNum-1 else PhScrInd+1])
            Hologram = update_PhasescreensMultimode(current_beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

    # IterationCount += 1
    # print(f"Current iteration: {IterationCount}")

    # # Forward and backward propagation
    # for direction in ['forward', 'backward']:
    #     if direction == 'forward':
    #         current_beam = Beam
    #         start_index, end_index, step = 0, PhaseScreenNum, 1
    #     else:
    #         current_beam = BeamBack
    #         start_index, end_index, step = PhaseScreenNum-1, -1, -1
    #     for PhScrInd in range (start_index, end_index, step):
    #         for jmodes in range(Nmodes):
    #             if direction == 'forward':
    #                 current_beam[:, :, :, jmodes] = elements(
    #                     current_beam[:, :, :, jmodes],
    #                     None, 
    #                     Arch_elements, 
    #                     Hologram[:, :, :, jmodes], 
    #                     KZ, 
    #                     L1, 
    #                     ElementShow, 
    #                     PhScrInd, 
    #                     hologram_position[PhScrInd], 
    #                     lens_position[-1 if PhScrInd == PhaseScreenNum-1 else PhScrInd+1]
    #                     )
    #             else:
    #                 print('aqui')
    #                 # print('is current_beam: ', current_beam[:, :, :, jmodes])
    #                 current_beam[:, :, :, jmodes], = elements(
    #                     None, 
    #                     current_beam[:, :, :, jmodes], 
    #                     Arch_elements, 
    #                     Hologram[:, :, :, jmodes], 
    #                     KZ, L1, 
    #                     ElementShow, 
    #                     PhScrInd, 
    #                     hologram_position[PhScrInd], 
    #                     hologram_position[-1]-1)
    #                     # lens_position[PhScrInd-1])
    #         # current_beam = elements(current_beam, None, Arch_elements, Hologram, KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[-1 if PhScrInd == PhaseScreenNum-1 else PhScrInd+1])
    #         Hologram = update_PhasescreensMultimode(current_beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)
             
    # Beam[:, :, :, jmodes] = elements(
    #     Beam=Beam[:, :, :, jmodes],
    #     BeamBack=None,  # Not used in forward propagation
    #     Arch_elements=Arch_elements,
    #     Hologram=Hologram[:, :, :, jmodes],
    #     KZ=KZ,
    #     L1=L1,
    #     ElementShow=ElementShow,
    #     PhScrInd=1,  # Starting at the first phase screen index for forward propagation
    #     start=1,
    #     finish=hologram_position[-1] + 1
    # )

    # Update overlap calculations and plot
    for jmodes in range(Nmodes):
        Overlap_to_output = np.real(np.abs(np.sum(np.conjugate(Beam[:, :, -1, jmodes]) * ModesOut[:, :, jmodes]))**2)
        resultoverlap[jmodes, IterationCount-1] = Overlap_to_output

        if lines:
            lines[jmodes].set_ydata(resultoverlap[jmodes, :])
        else:
            line, = ax.plot(resultoverlap[jmodes, :], label=f'Overlap Mode {jmodes+1}')
            lines.append(line)

    ax.legend(loc='lower right')
    ax.set_ylabel('Overlap')
    ax.set_xlabel('Iteration')
    ax.grid(True)
    canvas.draw()
    root.update_idletasks()

    # Check for convergence
    if IterationCount > convValue and np.all(np.round(resultoverlap[:, IterationCount-1], convLim) <= np.round(np.mean(resultoverlap[:, IterationCount-convValue:IterationCount], axis=1), convLim)):
        print("Convergence criteria met.")
        break

# Ensure GUI cleanup
root.destroy()



# while IterationCount < IterMax and not stop_pressed:
#     IterationCount += 1
#     print(f"Current iteration: {IterationCount}")

#     # Propagate and update in the forward direction
#     for PhScrInd in range(PhaseScreenNum):
#         if ForwardUpdate:
#             print(f"Updating phase screen {PhScrInd} in the forward direction")

#             update_PhasescreensMultimode(Beam, BeamBack, Hologram, 0, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

#         if PhScrInd == PhaseScreenNum:
#             Beam[:,:,:,jmodes] = elements(Beam[:,:,:,jmodes], BeamBack[:,:,:,jmodes], Arch_elements, Hologram[:,:,:,jmodes], KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[-1])# , Nonlinearity)
#         else:
#             Beam[:,:,:, jmodes] = elements(Beam[:,:,:,jmodes], BeamBack[:,:,:,jmodes], Arch_elements, Hologram[:,:,:,jmodes], KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[PhScrInd+1]) #, Nonlinearity)

#     # Propagate and update in the backward direction
#     for PhScrInd in range(PhaseScreenNum-1, -1, -1):
#         if BackUpdate:
#             print(f"Updating phase screen {PhScrInd} in the backward direction")
#             update_PhasescreensMultimode(Beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

#         for jmodes in range(Nmodes):
#             if PhScrInd > 1:
#                 BeamBack[:, :, :, jmodes] = elements(Beam, BeamBack, Arch_elements, Hologram, KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[PhScrInd-1], Nonlinearity)

#     # Beam[:, :, :, jmodes] = elements(
#     #     Beam=Beam[:, :, :, jmodes],
#     #     BeamBack=None,  # Not used in forward propagation
#     #     Arch_elements=Arch_elements,
#     #     Hologram=Hologram[:, :, :, jmodes],
#     #     KZ=KZ,
#     #     L1=L1,
#     #     ElementShow=ElementShow,
#     #     PhScrInd=1,  # Starting at the first phase screen index for forward propagation
#     #     start=1,
#     #     finish=hologram_position[-1] + 1
#     # )

# # while IterationCount < IterMax and not stop_pressed:
# #     IterationCount += 1
# #     print(f"Current iteration: {IterationCount}")

# #     # Forward and backward propagation
# #     for direction in ['forward', 'backward']:
# #         if direction == 'forward':
# #             current_beam = Beam
# #             start_index, end_index, step = 0, PhaseScreenNum, 1
# #         else:
# #             current_beam = BeamBack
# #             start_index, end_index, step = PhaseScreenNum-1, -1, -1

# #         for PhScrInd in range(start_index, end_index, step):
# #             elements(current_beam, None, Arch_elements, Hologram, KZ, L1, ElementShow, PhScrInd, hologram_position[PhScrInd], lens_position[-1 if PhScrInd == PhaseScreenNum-1 else PhScrInd+1])
# #             update_PhasescreensMultimode(current_beam, BeamBack, Hologram, PhScrInd, nx, ny, Nmodes, CoarsePixelation, NPixelsCoarse, MaxPhaseValue)

#     # Update overlap calculations and plot
#     for jmodes in range(Nmodes):
#         Overlap_to_output = np.real(np.abs(np.sum(np.conjugate(Beam[:, :, -1, jmodes]) * ModesOut[:, :, jmodes]))**2)
#         resultoverlap[jmodes, IterationCount-1] = Overlap_to_output

#         if lines:
#             lines[jmodes].set_ydata(resultoverlap[jmodes, :])
#         else:
#             line, = ax.plot(resultoverlap[jmodes, :], label=f'Overlap Mode {jmodes+1}')
#             lines.append(line)

#     ax.legend(loc='lower right')
#     ax.set_ylabel('Overlap')
#     ax.set_xlabel('Iteration')
#     ax.grid(True)
#     canvas.draw()
#     root.update_idletasks()

#     # Check for convergence
#     if IterationCount > convValue and np.all(np.round(resultoverlap[:, IterationCount-1], convLim) <= np.round(np.mean(resultoverlap[:, IterationCount-convValue:IterationCount], axis=1), convLim)):
#         print("Convergence criteria met.")
#         break


# print('hologram: ', Hologram)

# Display the phase screens
plt.figure()
for PhScrInd in range(PhaseScreenNum):
    ax = plt.subplot(1, PhaseScreenNum, PhScrInd + 1)
    im = plt.imshow(np.angle(Hologram[:, :, PhScrInd]), cmap='hsv')  # Display phase of the hologram
    plt.axis('square')
    plt.title(f'Hologram {PhScrInd + 1}')
    plt.xticks([])  # Remove x ticks
    plt.yticks([])  # Remove y ticks

# Adjusting the colorbar to show phase values
plt.colorbar(im, ax=ax, orientation='vertical', ticks=[0, np.pi, 2*np.pi], format='$%0.2f$')
plt.gca().set_position(ax.get_position())  # Reset the position of the last axis to original

# Maximizing the figure window - not directly possible like MATLAB, so using `plt.get_current_fig_manager()`
plt.get_current_fig_manager().window.state('zoomed')

# Plot the simulated conversion of input and output
plt.figure(figsize=(10, 5 * Nmodes))  # Adjust size as needed
for jmodes in range(Nmodes):
    # BeamBefore = ModesOut[:, :, jmodes]
    ax = plt.subplot(Nmodes, 2, 2 * jmodes + 1)
    # amplitude_phase_plot(BeamBefore, ax)
    # plt.title('Initial mode')
    # plt.axis('square')

    # Input modes - this is working
    BeamBefore = ModesOut[:, :, jmodes]
    #ax = axs[jmodes, 0] if Nmodes > 1 else axs[0]
    img = amplitude_phase_plot(BeamBefore, ax)
    ax.set_title('Initial Mode')
    # fig.colorbar(img, ax=ax, orientation='vertical', label='Phase')


    ax = plt.subplot(Nmodes, 2, 2 * (jmodes + 1))
    BeamAfter = Beam[:, :, PhaseScreenNum, jmodes]
    # amplitude_phase_plot(BeamAfter, ax)
    # plt.title('Simulated output')
    # plt.axis('square')
    img = amplitude_phase_plot(BeamAfter, ax)
    ax.set_title('Simulated Output')
    # fig.colorbar(img, ax=ax, orientation='vertical', label='Phase')


plt.tight_layout()
plt.show()

# Display final average overlap
#Overlap_to_output = np.real(np.sum(np.conjugate(Beam[:, :, -1, :]) * ModesOut, axis=(0, 1)))  # Summing over spatial dimensions
# print(f"\nOverlap to output: {Overlap_to_output}")
# final_average_overlap = np.trace(Overlap_to_output) / Nmodes  # Assuming you want the average per mode
# nota: sera que o overlap nao esta a ser definido como uma matriz no codigo
# matlab? parece que nao mas depois eles fazem o trace para obter o overlap medio
# no entanto esse valor da igual ao overlap normal
print(f"\nFinal average overlap is {Overlap_to_output}\n")
