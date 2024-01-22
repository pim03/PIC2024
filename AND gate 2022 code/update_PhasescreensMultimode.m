%% Subscript for updating phasescreens 
% using current Beam and BeamBack values.
% Not a function since Beam and BeamBack are huge, and we're not sure
% if matlab uses pointers in functions.

DeltaHologram = zeros(nx,ny,Nmodes);

for jmodes=1:Nmodes
    % Calculate overlap o_kii (BeamBack was conjugated in the
    % initialization)
    Overlap = BeamBack(:,:,PhScrInd,jmodes).*Beam(:,:,PhScrInd,jmodes).*exp(1i*Hologram(:,:,PhScrInd));
    % Average phase
    AvePhase = mean(angle(Overlap));
    % Values for updating the Hologram (For a single mode pair):
    DeltaHologram = DeltaHologram + Overlap.*exp(-1i*AvePhase);
end

% Calculate the total change required for this hologram:
DeltaHoloSum = -angle(sum(DeltaHologram,3));
% New Hologram
Hologram(:,:,PhScrInd) = Hologram(:,:,PhScrInd)+DeltaHoloSum;

% Include phase resolution
Hologram(:,:,PhScrInd) = mod(Hologram(:,:,PhScrInd),2*pi);

% Make the SLM grid coarser
if CoarsePixelation
    NCellsTiles = round(size(X,1)/NPixelsCoarse);
    Tiles = mat2tiles(Hologram(:,:,PhScrInd), [NCellsTiles,NCellsTiles]);
    PixelPhase = cell2mat( cellfun(@(c) c*0+mean(c(:)) , Tiles, 'uni',0  ));
    Hologram(:,:,PhScrInd) = PixelPhase;
end

% Normalize, scale to max phase value and discretize
Hologram2(:,:,PhScrInd) = floor(Hologram(:,:,PhScrInd)/max(max(Hologram(:,:,PhScrInd)))*MaxPhaseValue);
% Back to radians
Hologram(:,:,PhScrInd) = Hologram2(:,:,PhScrInd)/max(max(Hologram2(:,:,PhScrInd)))*2*pi;


