import numpy as np
import matplotlib.pyplot as plt
import AddPhaseColorbar
import LaguerreGen
import amplitude_phase_plot
import intensity_phase_plot
import LG
import SplitStepProp
from scipy import special
import math
from matplotlib.colors import LinearSegmentedColormap

# if 'C:\\Users\\mpiccardo\\' in os.listdir('C:\\Users\\'):
#     fproot = 'C:\\Users\\mpiccardo\\'
# else:
#     fproot = 'C:\\Users\\marco\\'

# RGBVincent = np.load(fproot + '\\Dropbox (Harvard University)\\Postdoc\\Vincent\\Near-field landscape\\RGB_Vincent.npy')

# Normalize and create colormap
RGBVincent = np.array([
    [[255, 0, 0], [0, 255, 0], [0, 0, 255]],  # Red, Green, Blue
    [[255, 255, 0], [255, 0, 255], [0, 255, 255]]  # Yellow, Magenta, Cyan
]) / 255.0  # Normalize to [0, 1]

# Create a custom colormap from the RGB array
CT = LinearSegmentedColormap.from_list("custom_cmap", RGBVincent.reshape(-1, 3))


# CT = RGBVincent
# CT = np.flipud(RGBVincent)
Path = ['P7']
N_apx = 1
N_apy = 1
L_ap = 3
d_ap = 0.3
xw = 5
yw = xw
xw = int(xw/d_ap)*d_ap
yw = int(yw/d_ap)*d_ap
lambda_ = 1064e-6
Npoints = 30
dx = d_ap/Npoints
x = np.arange(-xw/2, xw/2+dx, dx)
Nx = len(x)
dy = d_ap/Npoints
y = np.arange(-yw/2, yw/2+dy, dy)
Ny = len(y)
X, Y = np.meshgrid(x, y)
n = 1
k = n*2*np.pi/lambda_
L_J = 1
phasetype = 15
arraytype = 0
defect = 0
L_def = 1
nomask = 0
tictactoe = 0
gaussian = 0
FFTtype = 2
Tukeywindow = 0
Tukeywinfrac = 0
Fouriermask = 1
intphase_plot = 0
FFmasktype = 1
R_M_F = 2
pixelating = 0
xpixels = 10*lambda_
Npixels = int(np.ceil(xpixels/dx))
azi_sec = 0
N_azi_sec = 1
theta_azi_sec = 90/180*np.pi
zT = 2*d_ap**2/lambda_
print('Defining mask...')
if nomask:
    M = np.ones((Nx, Ny))
else:
    M = np.zeros((Nx, Ny))
for i in range(N_apx):
    for j in range(N_apy):
        if arraytype == 0:
            X_AP = (i-1)*d_ap - (N_apx/2)*d_ap + d_ap/2
            Y_AP = (j-1)*d_ap - (N_apy/2)*d_ap + d_ap/2
        elif arraytype == 1:
            X_AP = (i-1 + 0.5*np.mod(j,2))*d_ap - (N_apx/2)*d_ap + d_ap/2
            Y_AP = (j-1)*d_ap*np.sqrt(3)/2 - (N_apy/2)*d_ap + d_ap/2
        elif arraytype == 2:
            X_AP = (i-1)*d_ap - (N_apx/2)*d_ap
            Y_AP = (j-1)*d_ap - (N_apy/2)*d_ap
        RHO = ((X-X_AP)**2 + (Y-Y_AP)**2)**0.5
        PHI = np.arctan2((Y-Y_AP),(X-X_AP))
        if phasetype == 0:
            M[ RHO <= L_ap/2 ] = 1
        elif phasetype == 1:
            M[ RHO <= L_ap/2 ] = 1*np.exp(1j*np.mod(i+j,2)*np.pi)
        elif phasetype == 2:
            M[ RHO <= L_ap/2 ] = 1*np.exp(1j*np.random.rand(1)*2*np.pi)
        elif phasetype == 3:
            if tictactoe:
                if i == 3 or i == N_apx-2 or j == 3 or j == N_apy-2:
                    M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])
            else:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])
        elif phasetype == 4:
            M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*np.mod(i+j,2)*np.pi)
        elif phasetype == 5:
            M[ RHO <= L_ap/2 ] = (np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ]) + np.exp(1j*3*L_J*PHI[ RHO <= L_ap/2 ]))
        elif phasetype == 6:
            M[ RHO <= L_ap/2 ] = np.exp(1j*np.pi/2)
        elif phasetype == 7:
            w0 = d_ap/8
            M[ RHO <= L_ap/2 ] = (RHO[ RHO <= L_ap/2 ] *np.sqrt(2)/w0)**abs(L_J)*(2*RHO[ RHO <= L_ap/2 ]**2/w0**2)*np.exp(-(RHO[ RHO <= L_ap/2 ] /w0)**2)*np.exp(1j*np.mod(i+j,2)*np.pi)
        elif phasetype == 8:
            if arraytype == 0:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*((-np.mod(i,2)+np.mod(j,2))*np.pi/2 + np.mod(i,2)*np.mod(j,2)*np.pi + np.pi/3))
            elif arraytype == 1:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*(np.mod(j,3)*2*np.pi/3 + np.pi/3))
        elif phasetype == 9:
            if arraytype == 0:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*((np.mod(i,2)-np.mod(j,2))*np.pi/2 + np.mod(i,2)*np.mod(j,2)*np.pi + np.pi/3))
            elif arraytype == 1:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*(-np.mod(j,3)*2*np.pi/3 + np.pi/3))
        elif phasetype == 10:
            PHIGLOBAL = np.arctan2((Y_AP),(X_AP))
            LGLOBAL = 1
            M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*LGLOBAL*PHIGLOBAL)
        elif phasetype == 11:
            PHIGLOBAL = np.arctan2((Y_AP),(X_AP))
            RHOGLOBAL = ((X_AP)**2 + (Y_AP)**2)**0.5
            LGLOBALext = 2
            LGLOBALint = 1
            if RHOGLOBAL <= N_apx*d_ap/5:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*LGLOBALint*PHIGLOBAL)
            else:
                M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*LGLOBALext*PHIGLOBAL)
        elif phasetype == 12:
            w0 = L_ap/5
            HG_L_index = 1
            HG_M_index = 0
            C_HG = np.sqrt(2/(np.pi*math.factorial(HG_L_index)*math.factorial(HG_M_index)))*2**(-(HG_L_index+HG_M_index)/2)
            M[ RHO <= L_ap/2 ] = C_HG * special.hermite(HG_L_index,np.sqrt(2)*(X[ RHO <= L_ap/2 ]-X_AP)/w0) * special.hermite(HG_M_index,np.sqrt(2)*(Y[ RHO <= L_ap/2 ]-Y_AP)/w0)*np.exp(-RHO[ RHO <= L_ap/2 ]**2/w0**2)
        elif phasetype == 13:
            M[ RHO <= L_ap/2 ] = (np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ]) + np.exp(-1j*L_J*PHI[ RHO <= L_ap/2 ]))*np.exp(1j*np.mod(i+j,2)*np.pi)
        elif phasetype == 14:
            M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*np.exp(1j*np.random.rand(1)*2*np.pi)
        elif phasetype == 15:
            w0 = 0.1
            M[ RHO <= L_ap/2 ] = np.exp(1j*L_J*PHI[ RHO <= L_ap/2 ])*(RHO[ RHO <= L_ap/2 ]/w0*np.sqrt(2))**abs(L_J)*np.exp(-RHO[ RHO <= L_ap/2 ]**2/(w0)**2)
        if gaussian:
            w0 = 0.65
            M[ RHO <= L_ap/2 ] = np.exp(-RHO[ RHO <= L_ap/2 ]**2/(w0)**2)*M[ RHO <= L_ap/2 ]
        if defect and i == round((N_apx+1)/2) and j == round((N_apy+1)/2):
            M[ RHO <= L_ap/2 ] = np.exp(1j*L_def*PHI[ RHO <= L_ap/2 ])*(RHO[ RHO <= L_ap/2 ]/w0*np.sqrt(2))**abs(L_def)*np.exp(-RHO[ RHO <= L_ap/2 ]**2/(w0)**2)
        if azi_sec:
            for jazi in range(N_azi_sec):
                if phasetype == 8 or phasetype == 9:
                    vortex_rot_angle = np.mod((-np.mod(i,2)+np.mod(j,2))*np.pi/2 + np.mod(i,2)*np.mod(j,2)*np.pi,2*np.pi)
                    M[ np.mod(PHI,2*np.pi) > (jazi-1)*2*np.pi/jazi + vortex_rot_angle and np.mod(PHI,2*np.pi) < (jazi-1)*2*np.pi/jazi + theta_azi_sec + vortex_rot_angle and RHO <= L_ap/2]  = 0
                else:
                    M[ np.mod(PHI,2*np.pi) > (jazi-1)*2*np.pi/jazi and np.mod(PHI,2*np.pi) < (jazi-1)*2*np.pi/jazi + theta_azi_sec and RHO <= L_ap/2]  = 0
print('Mask defined')
if pixelating:
    Tiles = np.array_split(M, Npixels, axis=0)
    M = np.block([[np.mean(tile) for tile in Tiles]]*Npixels)
if Tukeywindow:
    Tukeywx = np.hanning(len(x))
    Tukeywy = np.hanning(len(y))
    TUKEYWX, TUKEYWY = np.meshgrid(Tukeywx, Tukeywy)
    TUKEYW = TUKEYWX*TUKEYWY
else:
    TUKEYW = np.ones((X.shape[0], X.shape[1]))
xwlim = xw*(1-Tukeywinfrac)
ywlim = yw*(1-Tukeywinfrac)
fig = plt.figure(figsize=(15, 10))
ax1 = fig.add_subplot(231, projection='3d')
ax1.plot_surface(X, Y, np.abs(M)**2, cmap=CT)
ax1.view_init(0, 90)
ax1.set_xlabel('x (mm)')
ax1.set_ylabel('y (mm)')
ax1.set_zlabel('Intensity')
ax1.set_xlim([-xwlim/2, xwlim/2])
ax1.set_ylim([-ywlim/2, ywlim/2])
ax1.set_title('Mask intensity')
ax2 = fig.add_subplot(234, projection='3d')
ax2.plot_surface(X, Y, np.angle(M), cmap=CT)
ax2.view_init(0, 90)
ax2.set_xlabel('x (mm)')
ax2.set_ylabel('y (mm)')
ax2.set_zlabel('Phase')
ax2.set_xlim([-xwlim/2, xwlim/2])
ax2.set_ylim([-ywlim/2, ywlim/2])
ax2.set_title('Mask phase')
f1 = 150
L1 = np.exp(1j*2*np.pi/lambda_*(f1 - (X**2 + Y**2 + f1**2)**0.5))
L1[ X**2 + Y**2 > (xw/2)**2 ] = 0
NA1 = n*xw/(2*f1)
diff_lim1 = lambda_/(2*NA1)
f2 = 10
L2 = np.exp(1j*2*np.pi/lambda_*(f2 - (X**2 + Y**2 + f2**2)**0.5))
L2[ X**2 + Y**2 > (xw/2)**2 ] = 0
NA2 = n*xw/(2*f2)
diff_lim2 = lambda_/(2*NA2)
lambda_DG = 1064e-6
theta_DG = 0.5/180*np.pi
d_DG = lambda_DG/np.sin(theta_DG)
DG = (np.ones((X.shape[0], X.shape[1]))-np.random.rand(X.shape[0], X.shape[1])/1.5)*np.exp(1j*np.random.rand(X.shape[0], X.shape[1])*2*np.pi/10)
Tiles = np.array_split(DG, Npixels, axis=0)
DG = np.block([[np.mean(tile) for tile in Tiles]]*Npixels)
RHO = (X**2 + Y**2)**0.5
rho0 = 0.8
fAL = 5
AL = np.exp(1j*2*np.pi/lambda_*(fAL - ((RHO - rho0)**2 + fAL**2)**0.5))
AL[ X**2 + Y**2 > (xw/2)**2 ] = 0
if Fouriermask:
    L_MF = 1
    PHI_MF = np.arctan2(Y,X)
    M_F = np.exp(1j*L_MF*PHI_MF)
else:
    M_F = np.ones((M.shape[0], M.shape[1]))
L_g = 7
n_g = 1.82
d_g = 146
M_g = np.zeros((Nx, Ny))
M_g[ X**2 + Y**2 <= (L_g/2)**2 ] = 1
dx_defect = 0.1
Np_defect = round(dx_defect/dx)
DM = np.ones((X.shape[0], X.shape[1]))
DM[round(X.shape[0]/2)-round(Np_defect/2):round(X.shape[0]/2)+round(Np_defect/2),:] = 0
dkx = 2*np.pi*1/dx/Nx
kx = np.arange(-2*np.pi*1/dx/2, 2*np.pi*1/dx/2-dkx, dkx)
dky = 2*np.pi*1/dy/Ny
ky = np.arange(-2*np.pi*1/dy/2, 2*np.pi*1/dy/2-dky, dky)
KX, KY = np.meshgrid(kx, ky)
z1 = f1
z2 = f2
z3 = f2 - z2
z4 = f2
z5 = zT/2
z6 = fAL
z7 = 300
H1 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z1)
H1 = np.fft.ifftshift(H1)
H2 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z2)
H2 = np.fft.ifftshift(H2)
H3 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z3)
H3 = np.fft.ifftshift(H3)
H4 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z4)
H4 = np.fft.ifftshift(H4)
H5 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z5)
H5 = np.fft.ifftshift(H5)
H6 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z6)
H6 = np.fft.ifftshift(H6)
H7 = np.exp(1j*np.sqrt(k**2 - (KX**2 + KY**2))*z7)
H7 = np.fft.ifftshift(H7)
E = M.copy()

for jelem in Path:
    elements(jelem)

ax3 = fig.add_subplot(232, projection='3d')
ax3.plot_surface(X, Y, np.abs(E)**2, cmap=CT)
ax3.view_init(0, 90)
ax3.set_xlabel('x (mm)')
ax3.set_ylabel('y (mm)')
ax3.set_zlabel('Intensity')
ax3.set_xlim([-xwlim/2, xwlim/2])
ax3.set_ylim([-ywlim/2, ywlim/2])
ax3.set_title('Propagated field intensity')
ax4 = fig.add_subplot(235, projection='3d')
ax4.plot_surface(X, Y, np.angle(E), cmap=CT)
ax4.view_init(0, 90)
ax4.set_xlabel('x (mm)')
ax4.set_ylabel('y (mm)')
ax4.set_zlabel('Phase')
ax4.set_xlim([-xwlim/2, xwlim/2])
ax4.set_ylim([-ywlim/2, ywlim/2])
ax4.set_title('Propagated field phase')
if FFTtype == 0:
    E_FFT = np.fft.fftshift(np.fft.fft2(E))
    ax5 = fig.add_subplot(233, projection='3d')
    ax5.plot_surface(X, Y, np.abs(E_FFT)**2, cmap=CT)
    ax5.set_title('Propagated field FFT int.')
    ax5.view_init(0, 90)
    ax5.set_xlabel('x (mm)')
    ax5.set_ylabel('y (mm)')
    ax5.set_zlabel('Intensity')
    ax5.set_xlim([-xwlim/2, xwlim/2])
    ax5.set_ylim([-ywlim/2, ywlim/2])
    ax6 = fig.add_subplot(236, projection='3d')
    ax6.plot_surface(X, Y, np.angle(E_FFT), cmap=CT)
    ax6.set_title('Propagated field FFT ph.')
    ax6.view_init(0, 90)
    ax6.set_xlabel('x (mm)')
    ax6.set_ylabel('y (mm)')
    ax6.set_zlabel('Phase')
    ax6.set_xlim([-xwlim/2, xwlim/2])
    ax6.set_ylim([-ywlim/2, ywlim/2])
else:
    pass
if intphase_plot:
    plt.figure()
    pm = np.mod(np.angle(E),2*np.pi)/(2*np.pi)
    sat = np.ones(pm.shape)
    mm = np.abs(E)**2/np.max(np.max(np.abs(E)**2))
    hsv_im = np.dstack((pm, sat, mm))
    rgb_im = plt.cm.hsv(hsv_im)
    plt.imshow(rgb_im, extent=[-xwlim/2, xwlim/2, -ywlim/2, ywlim/2], aspect='auto')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Intensity-phase')

def elements(jelement):
    global E
    if jelement == 'M_g':
        E = M_g*E
    elif jelement == 'G':
        G = G0/(1+np.abs(E)**2/Isat)
        E = G*E
    elif jelement == 'L1':
        E = L1*E
    elif jelement == 'L2':
        E = L2*E
    elif jelement == 'AL':
        E = AL*E
    elif jelement == 'DG':
        E = DG*E
    elif jelement == 'S':
        E = SL*SL*E
    elif jelement == 'M':
        E = E*M
    elif jelement == 'DM':
        E = E*DM
    elif jelement == 'M_F':
        E = M_F*E
    elif jelement == 'Ir':
        E = Ir*E
    elif jelement == 'OC':
        E = R*E
    elif jelement == 'Mi':
        E = np.fliplr(E)
    elif jelement == 'FFT':
        E = np.fft.fftshift(np.fft.fft2(E))
    elif jelement == 'IFFT':
        E = np.fft.ifft2(np.fft.ifftshift(E))
    elif jelement == 'P1':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H1*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P2':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H2*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P3':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H3*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P4':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H4*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P5':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H5*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P6':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H6*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P7':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H7*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    elif jelement == 'P8':
        E_FFT = np.fft.fft2(E)
        E2_FFT = H8*E_FFT
        E2 = np.fft.ifft2(E2_FFT)
        E = E2*TUKEYW
    else:
        print(jelement, 'does not exist. Please correct')

def Hermite_pol(n, x=None):
    h = np.zeros(n+1)
    n_fact = np.math.factorial(n)
    for m in range(0, int(np.floor(n/2))+1):
        h[2*m+1] = n_fact * (-1)**m / (np.math.factorial(m) * np.math.factorial(n-2*m)) * 2**(n-2*m)
    if x is not None:
        h = np.polyval(h, x)
    return h


