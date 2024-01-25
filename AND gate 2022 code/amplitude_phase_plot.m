function amplitude_phase_plot(Mode_Efield)
%% Amplitude_phase_plot
% Displays the transverse electric field's phase with a colormap and 
% the amplitude overlayed as alpha (transparency) values, along with a 
% black background.

    Amp = abs(Mode_Efield);
    Amp = Amp/max(max(Amp));
    Ang = angle(Mode_Efield);
    
    % Colormap for phase:
    colormap(hsv)
    imag = imagesc(Ang,[-pi pi]);
    
    %colormap(hsv) % Optional colormap 
    set(imag, 'AlphaData', Amp) %overlay with amplitude
    set(gca,'color',[0 0 0]) %set background black
    
    % Get rid of axes:
    set(gca,'xtick',[])
    set(gca,'ytick',[]) 

end