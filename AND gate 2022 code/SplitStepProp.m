function BeamProp = SplitStepProp(Beam,KZ,distance)
%%  Propagate beam via FFT
% - "Beam" is the transverse structure
% - "KZ" is a grid of wave vector components in the propagation direction
%   (this is a k-space grid and is dependant on the wavelength)
%   (Note: This grid needs to be in the fftshifted form where the k-space 
%   origin is in the corners of the grid)
% - "distance" is the propagation distance we want to simulate the propagation over
% The output is the "Beam"'s transverse structure after a propagation "distance"

if evalin( 'base', 'Nonlinearity == 1' )
    BeamFFT = fft2(Beam);
    BeamK = BeamFFT.*exp(-1i*(KZ*distance/2)); %distance/2 --> nonlinear layer between holograms  
    BeamProp = ifft2(BeamK);
    
    %Saturable absorber
%     Isat = 1e-4; %saturation intensity
%     Tmin = 0.8; %transmission for zero incident intensity
%     Tmax = 0.99; %transmission at saturation intensity
%     NL = Tmin + abs(BeamProp).^2/Isat*(Tmax - Tmin); %nonlinear transmission matrix
%     NL(NL > 1) = 1; %impose that maximum transmission is unity
%     BeamProp = BeamProp.*NL;

    %Photorefractive crystal
    NL = pi*abs(BeamProp).^2/(1 + abs(BeamProp).^2);
    BeamProp = BeamProp.*exp(-1i*NL);
    
    BeamFFT = fft2(BeamProp);
    BeamK = BeamFFT.*exp(-1i*(KZ*distance/2));    
    BeamProp = ifft2(BeamK);
else
    BeamFFT = fft2(Beam);
    BeamK = BeamFFT.*exp(-1i*(KZ*distance));    
    BeamProp = ifft2(BeamK);
end