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
    Isat = 0.5e-4; %saturation intensity - optical intensity (power per unit area) that it takes in a steady state to reduce the absorption to half of its unbleached value
    Tmin = 0.8; %transmission for zero incident intensity
    Tmax = 0.99; %transmission at saturation intensity
    alpha_loss = 0.01; %difference from total transmition
    %NL = Tmin + abs(BeamProp).^2/Isat*(Tmax - Tmin); %nonlinear transmission matrix
    %NL(NL > 1) = 1; %impose that maximum transmission is unity
    NL = 1 - (Tmax-Tmin)*exp(-abs(BeamProp).^2/Isat) - alpha_loss; %from paper "https://www.researchgate.net/publication/336731959_Study_of_a_Graphene_Saturable_Absorber_Film_Fabricated_by_the_Optical_Deposition_Method"
    BeamProp = BeamProp.*NL;

    %Photorefractive crystal
    %NL = pi*abs(BeamProp).^2/(1 + abs(BeamProp).^2);
    %BeamProp = BeamProp.*exp(-1i*NL);
    
    BeamFFT = fft2(BeamProp);
    BeamK = BeamFFT.*exp(-1i*(KZ*distance/2));    
    BeamProp = ifft2(BeamK);
else
    BeamFFT = fft2(Beam);
    BeamK = BeamFFT.*exp(-1i*(KZ*distance));    
    BeamProp = ifft2(BeamK);
end