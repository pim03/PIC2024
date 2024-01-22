clear all
% close all

if exist('C:\Users\mpiccardo\', 'dir')
   fproot = 'C:\Users\mpiccardo\';
else
   fproot = 'C:\Users\marco\';
end

load([fproot,'\Dropbox (Harvard University)\Postdoc\Vincent\Near-field landscape\RGB_Vincent.mat'])
CT = RGBVincent;
CT = flipud(RGBVincent);

% Path = {'DG','P1','L1','P1','P1','L1','P1'};
% Path = {'DG','P1','L1','P1'};

% Path = {'DG','P1','L1','P1','M_F','P1','L1','P1'};
% Path = {'DG','P1','L1','P1','M_F'};

% Path = {'P1','L1','P1','P1','L1','P1'};
Path = {'P7'};

% Path = {'P1','L1','P1','M_F','P1','L1','P1'};
% Path = {'P1','L1','P1','M_F','P1','L1','P1','Mi','P1','L1','P1','P1','L1','P1'};
% Path = {'P1','L1','P1','P1','L1','P1','Mi','P1','L1','P1','P1','L1','P1'};

%% Mask parameters

N_apx = 1; %number of apertures in x direction
N_apy = 1; %number of apertures in y direction
L_ap = 3; % [mm] aperture diameter
d_ap = 0.3; %[mm] distance between the mask apertures

%% Spatial parameters

xw = 5; %[mm] x width of window
yw = xw; %[mm] y width of window

xw = fix(xw/d_ap)*d_ap; % Adjusts window width to get best FFT with periodic mask
yw = fix(yw/d_ap)*d_ap; % Adjusts window width to get best FFT with periodic mask

%% Field parameters

lambda = 1064e-6; %[mm] wavelength
Npoints = 30;
dx = d_ap/Npoints;
x = [-xw/2:dx:xw/2];
Nx = length(x);
dy = d_ap/Npoints;
y = [-yw/2:dy:yw/2];
Ny = length(y);
[X,Y] = meshgrid(x,y);

n = 1; %refractive index of propagation medium
k = n*2*pi/lambda; %wavenumber

%% Define mask

L_J = 1;
phasetype = 15; %0 = in phase, 1 = out of phase, 2 = random, 3 = OAM in phase, 4 = OAM out of phase, 5 = test, 6 = phase mask,
                %7 = donuts with no charge out of phase, 8 = vortex CW, 9 = vortex CCW, 10 = in phase OAM with global OAM phase,
                %11 = in phase OAM with global concentric OAM phases, 12 = HG, 13 = superposition of opposite OAM modes out of phase
                %14 = OAM random phase, 15 = p=0 OAM mode
arraytype = 0; %0 = square with aperture on (0,0) coord; 1 = triangular; 2 = square with no aperture on (0,0) coord;
defect = 0; %if 1 it puts (some) defect(s) in the array
L_def = 1; %topological charge of the defect J-plate
nomask = 0; %if 1 it removes the mask: useful to simulate array of metasurfaces without physical mask
tictactoe = 0; %creates a tic tac toe illumination pattern for phasetype = 3 (test)
gaussian = 0; %gaussian amplitude distribution in each aperture
FFTtype = 2; %0, FFT of propagated field; 1, FFT of mask; 2, do not calculate FFT
Tukeywindow = 0; %0 = no absorbing layer at boundary; 1 = with absorbing frame represented by a Tukey window
Tukeywinfrac = 0; %fraction of window that is used for the absorbing layer, e.g. 0 = no absorbing layer, 0.5 = half window is absorbing
Fouriermask = 1; %1 = Fourier mask
intphase_plot = 0; %1 = plot figure combining intensity and phase
FFmasktype = 1; %0 = negative FF mask (inner circle blocks light); 1 = positive FF mask (inner circle lets light pass); 2 = load a FF mask
R_M_F = 2; %[mm] radius of the Fourier mask
pixelating = 0; % 1 = discretize spatially
xpixels = 10*lambda; %[mm] size of 1 pixel
Npixels = ceil(xpixels/dx);
azi_sec = 0; % 1 = switch on azimuthal sector mask
N_azi_sec = 1; %number of azi sectors
theta_azi_sec = 90/180*pi; %[rad] width of the azi sector

zT = 2*d_ap^2/lambda; %Talbot length

disp('Defining mask...')

%Mask from array
if nomask
    M = ones(Nx,Ny);
else
    M = zeros(Nx,Ny);
end

for i=1:N_apx
    for j=1:N_apy
        switch arraytype
            case 0
                X_AP(i,j) = (i-1)*d_ap - (N_apx/2)*d_ap + d_ap/2;
                Y_AP(i,j) = (j-1)*d_ap - (N_apy/2)*d_ap + d_ap/2;
            case 1
                X_AP(i,j) = (i-1 + 0.5*mod(j,2))*d_ap - (N_apx/2)*d_ap + d_ap/2;
                Y_AP(i,j) = (j-1)*d_ap*sqrt(3)/2 - (N_apy/2)*d_ap + d_ap/2;
            case 2
                X_AP(i,j) = (i-1)*d_ap - (N_apx/2)*d_ap;
                Y_AP(i,j) = (j-1)*d_ap - (N_apy/2)*d_ap;
        end
        
        RHO = ((X-X_AP(i,j)).^2 + (Y-Y_AP(i,j)).^2).^0.5;
        PHI = atan2((Y-Y_AP(i,j)),(X-X_AP(i,j)));
        
        switch phasetype
            case 0
                M( RHO <= L_ap/2 ) = 1;
            case 1
                M( RHO <= L_ap/2 ) = 1*exp(1i*mod(i+j,2)*pi);    
            case 2
                M( RHO <= L_ap/2 ) = 1*exp(1i*rand(1)*2*pi);
            case 3
                if tictactoe
                    if i == 3 || i == N_apx-2 || j == 3 || j == N_apy-2
                        M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ));
                    end
                else
                    M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ));
                end
            case 4
                M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*mod(i+j,2)*pi);
            case 5
                M( RHO <= L_ap/2 ) = (exp(1i*L_J*PHI( RHO <= L_ap/2 )) + exp(1i*3*L_J*PHI( RHO <= L_ap/2 ))); %exp(-RHO( RHO <= L_ap/2 ).^2/0.05^2).*
            case 6
                M( RHO <= L_ap/2 ) = exp(1i*pi/2);
%                 M( RHO <= L_ap/2 ) = exp(1i*76/180*pi);
            case 7
                w0 = d_ap/8; %beam waist of the embedded Gaussian
                M( RHO <= L_ap/2 ) = (RHO( RHO <= L_ap/2 ) *sqrt(2)/w0).^abs(L_J).*(2*RHO( RHO <= L_ap/2 ).^2/w0^2).*exp(-(RHO( RHO <= L_ap/2 ) /w0).^2)*exp(1i*mod(i+j,2)*pi);
            case 8
                switch arraytype
                    case 0
                        M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*((-mod(i,2)+mod(j,2))*pi/2 + mod(i,2)*mod(j,2)*pi + pi/3)); % "+ pi/3" is an arbitrary phase shift to distinguish phase of aperture from background
                    case 1
                        M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*(mod(j,3)*2*pi/3 + pi/3)); % "+ pi/3" is an arbitrary phase shift to distinguish phase of aperture from background
                end
            case 9
                switch arraytype
                    case 0
                        M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*((mod(i,2)-mod(j,2))*pi/2 + mod(i,2)*mod(j,2)*pi + pi/3)); % "+ pi/3" is an arbitrary phase shift to distinguish phase of aperture from background
                    case 1
                        M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*(-mod(j,3)*2*pi/3 + pi/3)); % "+ pi/3" is an arbitrary phase shift to distinguish phase of aperture from background
                end
            case 10
                PHIGLOBAL = atan2((Y_AP(i,j)),(X_AP(i,j)));
                LGLOBAL = 1;
                M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*LGLOBAL*PHIGLOBAL);
            case 11 %concentric global charges
                PHIGLOBAL = atan2((Y_AP(i,j)),(X_AP(i,j)));
                RHOGLOBAL = ((X_AP(i,j)).^2 + (Y_AP(i,j)).^2).^0.5;
                LGLOBALext = 2;
                LGLOBALint = 1;
                if RHOGLOBAL <= N_apx*d_ap/5
                    M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*LGLOBALint*PHIGLOBAL);
                else
                    M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*LGLOBALext*PHIGLOBAL);
                end
            case 12
                %Based on Eqs. 1 and 6 in M. W. Beijersbergen et al., Optics Communications 96 (1993) 123-132
                %https://doi.org/10.1016/0030-4018(93)90535-D
                w0 = L_ap/5; %beam waist
                HG_L_index = 1;
                HG_M_index = 0;
                C_HG = sqrt(2/(pi*factorial(HG_L_index)*factorial(HG_M_index)))*2^(-(HG_L_index+HG_M_index)/2);
                M( RHO <= L_ap/2 ) = C_HG*Hermite_pol(HG_L_index,sqrt(2)*(X( RHO <= L_ap/2 )-X_AP(i,j))/w0).*Hermite_pol(HG_M_index,sqrt(2)*(Y( RHO <= L_ap/2 )-Y_AP(i,j))/w0).*exp(-RHO( RHO <= L_ap/2 ).^2/w0^2);
            case 13
                M( RHO <= L_ap/2 ) = (exp(1i*L_J*PHI( RHO <= L_ap/2 )) + exp(-1i*L_J*PHI( RHO <= L_ap/2 )))*exp(1i*mod(i+j,2)*pi);
            case 14
                M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 ))*exp(1i*rand(1)*2*pi);
            case 15
                w0 = 0.1; %beam waist
                M( RHO <= L_ap/2 ) = exp(1i*L_J*PHI( RHO <= L_ap/2 )).*(RHO( RHO <= L_ap/2 )/w0*sqrt(2)).^abs(L_J).*exp(-RHO( RHO <= L_ap/2 ).^2/(w0)^2);
        end
        
        if gaussian
            w0 = 0.65; %[mm] beam waist
            M( RHO <= L_ap/2 ) = exp(-RHO( RHO <= L_ap/2 ).^2/(w0)^2).*M( RHO <= L_ap/2 );
        end
        
        if defect && i == round((N_apx+1)/2) && j == round((N_apy+1)/2)  %single defect at center
            M( RHO <= L_ap/2 ) = exp(1i*L_def*PHI( RHO <= L_ap/2 )).*(RHO( RHO <= L_ap/2 )/w0*sqrt(2)).^abs(L_def).*exp(-RHO( RHO <= L_ap/2 ).^2/(w0)^2);
%             M( RHO <= L_ap/2 ) = exp(1i*L_def*PHI( RHO <= L_ap/2 )); %defect with different charge L_def
        end
        
        if azi_sec
            for jazi = 1:N_azi_sec
                if phasetype == 8 || phasetype == 9
                    vortex_rot_angle = mod((-mod(i,2)+mod(j,2))*pi/2 + mod(i,2)*mod(j,2)*pi,2*pi);
                    M( mod(PHI,2*pi) > (jazi-1)*2*pi/jazi + vortex_rot_angle & mod(PHI,2*pi) < (jazi-1)*2*pi/jazi + theta_azi_sec + vortex_rot_angle & RHO <= L_ap/2)  = 0;
                else
                   M( mod(PHI,2*pi) > (jazi-1)*2*pi/jazi & mod(PHI,2*pi) < (jazi-1)*2*pi/jazi + theta_azi_sec & RHO <= L_ap/2)  = 0;
                end
            end
        end
        
    end
end

disp('Mask defined')

%% Negative inverse of mask
% M(M==0) = 2;
% M(M==1) = 0;
% M(M==2) = 1;

%% Spatial sampling

if pixelating
%     Tiles = mat2tiles(angle(M), [Npixels,Npixels]);
%     PixelPhase = cell2mat( cellfun(@(c) c*0+mean(c(:)) , Tiles, 'uni',0  ));
%     M = exp(1i*PixelPhase);
    Tiles = mat2tiles(M, [Npixels,Npixels]);
    M = cell2mat( cellfun(@(c) c*0+mean(c(:)) , Tiles, 'uni',0  ));
end

%% Tukey window

if Tukeywindow
    Tukeywx = window(@tukeywin,length(x),Tukeywinfrac);
    Tukeywy = window(@tukeywin,length(y),Tukeywinfrac);
    [TUKEYWX,TUKEYWY] = meshgrid(Tukeywx,Tukeywy);
    TUKEYW = TUKEYWX.*TUKEYWY;
else
    TUKEYW = ones(size(X,1),size(X,2));
end

xwlim = xw*(1-Tukeywinfrac); %size of window without absorbing frame for plot limits
ywlim = yw*(1-Tukeywinfrac); %size of window without absorbing frame for plot limits

%% Plot mask

figure
set(gcf, 'Position', get(0, 'Screensize'));
sax1=subplot(231);
surf(X,Y,abs(M).^2)
view([0,90])
shading flat
axis square
xlim([-xwlim/2 xwlim/2])
ylim([-ywlim/2 ywlim/2])
xlabel('x (mm)')
ylabel('y (mm)')
colorbar
colormap(CT)
title('Mask intensity')

sax3=subplot(234);
surf(X,Y,angle(M))
view([0,90])
shading flat
axis square
xlim([-xwlim/2 xwlim/2])
ylim([-ywlim/2 ywlim/2])
xlabel('x (mm)')
ylabel('y (mm)')
colorbar
colormap(CT)
title('Mask phase')

%% Lenses

f1 = 150; %[mm] lens focal length
L1 = exp(1i*2*pi/lambda*(f1 - (X.^2 + Y.^2 + f1^2).^0.5)); %spherical lens
L1( X.^2 + Y.^2 > (xw/2)^2 ) = 0; %To make the lens circular
NA1 = n*xw/(2*f1); %numerical aperture
diff_lim1 = lambda/(2*NA1); %diffraction limit

f2 = 10; %[mm] lens focal length
L2 = exp(1i*2*pi/lambda*(f2 - (X.^2 + Y.^2 + f2^2).^0.5));
L2( X.^2 + Y.^2 > (xw/2)^2 ) = 0; %To make the lens circular
NA2 = n*xw/(2*f2); %numerical aperture
diff_lim2 = lambda/(2*NA2); %diffraction limit

%% Diffraction grating

lambda_DG = 1064e-6; %[mm] central design wavelength
theta_DG = 0.5/180*pi; %diffraction angle of first order for central design wavelength
d_DG = lambda_DG/sin(theta_DG); %period of grating

% DG = besselj(0,2*pi/d_DG*(X.^2+Y.^2).^0.5); %circular Bessel diffraction grating
% DG = exp(1i*2*pi/d_DG*(X.^2+Y.^2).^0.5); %axicon grating
% T_PVB = lambda*f1/(pi*w0); %thickness of perfect vortex beam
% R_PVB = f1*lambda/d_DG; %radius of perfect vortex beam

% DG = sin(0,2*pi/d_DG*(X.^2+Y.^2).^0.5); %circular sinusoidal diffraction grating
% DG = exp(1i*sin(2*pi/d_DG*(X.^2+Y.^2).^0.5)); %circular diffraction grating
% DG = abs(sin(2*pi/d_DG*X)); %1D diffraction grating
% DG = sin(2*pi/d_DG*X); %1D sinusoidal diffraction grating
% DG = exp(1i*2*pi/d_DG*Y); %1D diffraction grating
DG = (ones(size(X,1),size(X,2))-rand(size(X,1),size(X,2))/1.5).*exp(1i*rand(size(X,1),size(X,2))*2*pi/10);

Tiles = mat2tiles(DG, [Npixels,Npixels]);
DG = cell2mat( cellfun(@(c) c*0+mean(c(:)) , Tiles, 'uni',0  ));


%% Annular lens

RHO = (X.^2 + Y.^2).^0.5;
rho0 = 0.8; %target radius of the transformed PVB
fAL = 5; %[mm] lens focal length
AL = exp(1i*2*pi/lambda*(fAL - ((RHO - rho0).^2 + fAL^2).^0.5));
AL( X.^2 + Y.^2 > (xw/2)^2 ) = 0; %To make the lens circular

%% Fourier mask
if Fouriermask
    L_MF = 1; %charge of Fourier mask
    PHI_MF = atan2(Y,X);
    M_F = exp(1i*L_MF*PHI_MF);
else
    M_F = ones(size(M,1),size(M,2));
end

% if sum(ismember(Path,'M_F'))
%     switch FFmasktype
%         case 0 %negative mask
%             M_F = ones(Nx,Ny);
%             M_F( (X.^2 + Y.^2).^0.5 < R_M_F ) = 0;
%         case 1 %positive mask
%             M_F = zeros(Nx,Ny);
%             M_F( (X.^2 + Y.^2).^0.5 < R_M_F ) = 1;
%         case 2 %load mask
%             FF_filepath = 'C:\Users\mpiccardo\Dropbox (Harvard University)\IIT\Lab\Talbot\FoxLi\masks\';
%             FF_filename = 'HG10_dots_quarter';
%             load([FF_filepath,FF_filename])
%             M_F = FF_mask;
%             if not(size(FF_mask,1) == size(M,1))
%                 disp('FF mask size is not suitable for this simulation.')
%                 disp('Regenerate mask making sure xw, Npoints, Tukeywinfrac are the same as in this simulation.')
%                 return
%             end
%     end
% else
%     M_F = zeros(Nx,Ny);
% end

%% Gain medium mask

L_g = 7; %[mm] diameter of the gain medium
n_g = 1.82; %index of the gain medium (Nd:YAG --> 1.82)
d_g = 146; %[mm] length of gain rod

M_g = zeros(Nx,Ny);
M_g( X.^2 + Y.^2 <= (L_g/2)^2 ) = 1; %aperture mask of the gain rod

%% Defect mask mimicking obscuration line of wedge mirror

dx_defect = 0.1; %[mm] line defect width
Np_defect = round(dx_defect/dx);

DM = ones(size(X,1),size(X,2));
DM(round(size(X,1)/2)-round(Np_defect/2):round(size(X,1)/2)+round(Np_defect/2),:) = 0;

% figure
% set(gcf,'Position',[1,41,2194.28571428571,1122.85714285714])
% surf(X,Y,DM)
% view([0,90])
% shading flat
% axis square
% xlim([-xw/2 xw/2])
% ylim([-yw/2 yw/2])
% xlabel('x (mm)')
% ylabel('y (mm)')
% colorbar
% title('Defect mask amplitude')

%% Propagate field

dkx = 2*pi*1/dx/Nx;
kx = -2*pi*1/dx/2:dkx:2*pi*1/dx/2-dkx;
dky = 2*pi*1/dy/Ny;
ky = -2*pi*1/dy/2:dky:2*pi*1/dy/2-dky;
[KX,KY] = meshgrid(kx,ky);

z1 = f1; % [mm] propagation distance 1
z2 = f2; % [mm]
z3 = f2 - z2; % [mm]
z4 = f2; % [mm]
z5 = zT/2; % [mm]
z6 = fAL; % [mm]
z7 = 300; % [mm]

%Propagator 1
H1 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z1); %optical transfer function to position z
H1 = ifftshift(H1);

%Propagator 2
H2 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z2); %optical transfer function to position z
H2 = ifftshift(H2);

%Propagator 3
H3 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z3); %optical transfer function to position z
H3 = ifftshift(H3);

%Propagator 4
H4 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z4); %optical transfer function to position z
H4 = ifftshift(H4);

%Propagator 5
H5 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z5); %optical transfer function to position z
H5 = ifftshift(H5);

%Propagator 6
H6 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z6); %optical transfer function to position z
H6 = ifftshift(H6);

%Propagator 7
H7 = exp(1i*sqrt(k^2 - (KX.^2 + KY.^2))*z7); %optical transfer function to position z
H7 = ifftshift(H7);

% Path = {'P1','L1','P1','P4','L2','P3','M_g','P2','M_g'};
% Path = {'P5'}; %Just Talbot propagation

% Path = {'DG', 'P1', 'L1', 'P1','M_F'}; %Diffraction grating, then propagation
% Path = {'DG', 'P1', 'L1', 'P1', 'M_F', 'L1', 'P1'}; %Diffraction grating, then metalens
% Path = {'DG','P1','L1', 'P1', 'M_F', 'L1', 'P1'}; %Diffraction grating, then metalens
% Path = {'P1','L1','P1'}; %FFT with lens
% Path = {'FFT'}; %Fourier transform
% Path = {'P5','DM','P6'}; %Talbot propagation with defect
% Path = {'P1','L1','P1','P4','L2','P3','P2'};

% if any(strcmp(Path,'DG'))
%     if d_DG/dx < 8
%        disp('x sampling is not fine enough for grating')
%        return
%     end
% end

E = M;

for jelem = 1:length(Path) 
    elements(Path{jelem});
end

sax2 = subplot(232);
surf(X,Y,abs(E).^2)
view([0,90])
shading flat
axis square
xlim([-xwlim/2 xwlim/2])
ylim([-ywlim/2 ywlim/2])
xlabel('x (mm)')
ylabel('y (mm)')
colorbar
colormap(CT)
title('Propagated field intensity')

sax4 = subplot(235);
surf(X,Y,angle(E))
view([0,90])
shading flat
axis square
xlim([-xwlim/2 xwlim/2])
ylim([-ywlim/2 ywlim/2])
xlabel('x (mm)')
ylabel('y (mm)')
colorbar
colormap(CT)
title('Propagated field phase')
% linkaxes([sax1,sax2,sax3,sax4],'xy')

switch FFTtype
    case 0 %FFT of propagated field
        E_FFT = fftshift(fft2(E));
        sax5 = subplot(233);
        surf(X,Y,abs(E_FFT).^2)
        title('Propagated field FFT int.')
        set(gcf,'Position',[1,41,2194.28571428571,1122.85714285714])
        view([0,90])
        shading flat
        axis square
        xlim([-xwlim/2 xwlim/2])
        ylim([-ywlim/2 ywlim/2])
        colorbar
        colormap(CT)
        
        sax6 = subplot(236);
        surf(X,Y,angle(E_FFT))
        title('Propagated field FFT ph.')
        view([0,90])
        shading flat
        axis square
        xlim([-xwlim/2 xwlim/2])
        ylim([-ywlim/2 ywlim/2])
        colorbar
        colormap(CT)
        linkaxes([sax1,sax2,sax3,sax4,sax5,sax6],'xy')
    case 1 %FFT of mask
        M_FFT = fftshift(fft2(M));
        sax5 = subplot(233);
        surf(X,Y,abs(M_FFT).^2)
        title('Mask FFT int.')
        view([0,90])
        shading flat
        axis square
        xlim([-xwlim/2 xwlim/2])
        ylim([-ywlim/2 ywlim/2])
        colorbar
        colormap(CT)

        sax6 = subplot(236);
        surf(X,Y,angle(M_FFT))
        title('Mask FFT ph.')
        view([0,90])
        shading flat
        axis square
        xlim([-xwlim/2 xwlim/2])
        ylim([-ywlim/2 ywlim/2])
        colorbar
        colormap(CT)
        linkaxes([sax1,sax2,sax3,sax4,sax5,sax6],'xy')
    case 2
        %not calculating FFT
end

if intphase_plot
    subplot(233)
    pm = mod(angle(E),2*pi)/(2*pi); %phase matrix normalized between 0 and 1
    sat = ones(size(pm)); %saturation matrix
    mm = abs(E).^2/max(max(abs(E).^2)); %intensity matrix normalized between 0 and 1
    hsv_im = cat(3,pm,sat,mm);
    rgb_im = hsv2rgb(hsv_im);
    image(x,y,rgb_im)
    xlim([-xwlim/2 xwlim/2])
    ylim([-ywlim/2 ywlim/2])
    axis square
    view([0 -90])
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Intensity-phase')
end

%% Cavity elements

function elements(jelement)

switch jelement
    case 'M_g' %Gain mask
        evalin('base',['E = M_g.*E;'])
    case 'G' %Gain
        evalin('base',['G = G0./(1+abs(E).^2/Isat);'])
        evalin('base',['E = G.*E;']) %To check: should the field or intensity be multiplied with gain?
    case 'L1' %Lens 1
        evalin('base',['E = L1.*E;'])
    case 'L2' %Lens 2
        evalin('base',['E = L2.*E;'])
    case 'AL' %Annular lens
        evalin('base',['E = AL.*E;'])
    case 'DG' %Lens 2
        evalin('base',['E = DG.*E;'])
    case 'S' %Spherical mirror
        evalin('base',['E = SL.*SL.*E;'])
    case 'M' %Mask
        evalin('base',['E = E.*M;'])
    case 'DM' %Defect mask
        evalin('base',['E = E.*DM;'])
    case 'M_F' %Fourier mask
        evalin('base',['E = M_F.*E;'])
    case 'Ir' %Iris mask
        evalin('base',['E = Ir.*E;'])
    case 'OC' %Output coupler
        evalin('base',['E = R*E;'])
    case 'Mi' %flat mirror
        evalin('base',['E = fliplr(E);'])
    case 'FFT' %Fourier transform
        evalin('base',['E = fftshift(fft2(E));'])
    case 'IFFT' %Inverse Fourier transform
        evalin('base',['E = ifft2(ifftshift(E));'])
    case 'P1' %Propagation 1
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H1.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P2' %Propagation 2
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H2.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P3' %Propagation 3
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H3.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P4' %Propagation 4
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H4.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P5' %Propagation 5
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H5.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P6' %Propagation 6
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H6.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P7' %Propagation 7
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H7.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    case 'P8' %Propagation 8
        evalin('base',['E_FFT = fft2(E);']) %Fourier transform of initial field
        evalin('base',['E2_FFT = H8.*E_FFT;'])
        evalin('base',['E2 = ifft2(E2_FFT);'])
        evalin('base',['E = E2.*TUKEYW;'])
    otherwise
        fprintf([jelement,' does not exist. Please correct']);
end

end

%% Hermite functions

%Based on "hermite_formula.png"
function h = Hermite_pol(n,x)
    h = zeros(1,n+1);
    n_fact = factorial(n);
    for m=0:floor(n/2) 
        h(2*m+1) = n_fact * (-1)^m / (factorial(m) * factorial(n-2*m)) * 2^(n-2*m); 
        % The argument of h is a vector of length n+1 whose elements are the coefficients
        % (in descending powers) of an nth-degree polynomial. Thus the polynomial
        % coefficient of x^(n-2*m) will be at position n+1-(n-2*m)=2*m+1 in vector h
    end
    if exist('x','var') 
        h = polyval (h, x); 
    end
end


%% Export Far Field mask example

% export_FF_mask(E,'C:\Users\mpiccardo\Dropbox (Harvard University)\IIT\Lab\Talbot\Results\FF masks\','HG10_dots')
