clear all
close all

cmap = cmocean('dense',255);

%% World parameters

Np = 128;
Niter = 10000000;
% Niter = 1;
Nstep_plot = 1; %after how many iterations it updates the world plot

x = linspace(-1,1,Np);
y = x;
[X,Y] = meshgrid(x,y);
RHO = (X.^2 + Y.^2).^0.5;
PHI = atan2(Y,X);

%% Input state

% A_noise = rand(size(X)); %input random matrix
% A_noise = conv2(A_noise,ones(5,5),'same');
% w_in = 0.3; x_in = 0; y_in = 0; %Gaussian parameters
% Ain = exp(-(((X-x_in).^2 + (Y-y_in).^2).^0.5).^2/(w_in)^2) + 0.05*A_noise; %input Gaussian

% Ain = rand(size(X)); %input random matrix
% Ain = conv2(Ain,ones(3,3),'same');
% Ain = Ain/max(max(Ain));

Ain = RandomPatches(zeros(Np,Np), 13, floor(Np/6));

w_sum = 0.1; %weight for the sum of initial and coupled state (Lenia: Deltat)

%% Kernel

% alpha = 4;
% K = exp(alpha*(1-1./(4*RHO.*(1-RHO)))); % (Lenia: Gaussian bump (raised to the power alpha))
% % K = (4*RHO.*(1-RHO)).^alpha; %
% K(RHO > 1) = 0;

w_K = 0.04; %beam waist
K = (RHO/w_K*sqrt(2)).*exp(-RHO.^2/(w_K)^2); %Kernel defined as LG mode
K = K/sum(sum(K)); %normalize Kernel using sum, then FFT becomes normalized to 1
fftK = fft2(K);

fh = figure;
colormap(cmap);
fh.WindowState = 'maximized';
pause(1)
subplot(231)
imagesc(x,y,K)
axis square; title('Kernel'); colorbar;
subplot(232)
imagesc(x,y,abs(fftshift(fftK)))
axis square; title('FFT(Kernel)'); colorbar

%% Nonlinear map    

w_map = 0.09; % (Lenia: sigma)
rho_map = 0.3; % (Lenia: mu)
% w_map = 0.017; % (Lenia: sigma)
% rho_map = 0.15; % (Lenia: mu)
amp_val = linspace(0,1,1000);
g_map = 2*exp(-(amp_val-rho_map).^2/(2*w_map^2)) - 1; % (Lenia: Delta d(n))
% g_map = exp(-(amp_val-rho_map).^2/(2*w_map^2)); % Special map [0,1]
subplot(233)
plot(amp_val,g_map); xlabel('Intensity'); ylabel('Field amplitude'); axis square

%% Iterations

subplot(234)
phandle1 = surf(x,y,Ain);
A = Ain;
phandle1.ZDataSource = 'abs(A)';
axis square; view([0 90]); shading flat; title('World'); colorbar
caxis([0 1])

subplot(235)
B = zeros(size(A));
phandle2 = surf(x,y,B);
phandle2.ZDataSource = 'B';
axis square; view([0 90]); shading flat; title('Potential'); colorbar
caxis([0 1])

subplot(236)
G = zeros(size(A));
phandle3 = surf(x,y,G);
phandle3.ZDataSource = 'G';
axis square; view([0 90]); shading flat; title('Field'); colorbar
caxis([0 1])

for jiter=1:Niter
    if mod(jiter,Nstep_plot) == 0
        disp(num2str(jiter))
        refreshdata(phandle1)
        refreshdata(phandle2)
        refreshdata(phandle3)
        pause(0.1)
    end

    B = ifft2(fftK.*fft2(A)); %potential

    % B = circshift(real(B), [floor(Np/2)+1, floor(Np/2)+1]); %it's in Lenia code but to me it makes no sense as the ifft2 after fft2 is not shifted!

    % B = B/max(max(abs(B))); %normalize convoluted output --> not in Lenia!
   
    G = 2*exp(-abs(B-rho_map).^2/(2*w_map^2)) - 1; %nonlinear map
    % G = B; %linear map

    A = A + w_sum*G; %sum weighted change to initial state

    A = min(max(abs(A),0),1); %clip A to range [0,1]

end

refreshdata(phandle1)


%% Functions

function world = RandomPatches(world, R, border)
    randSize = floor(R * 0.9);
    SIZE = size(world, 1);
    range = [border SIZE-border-randSize];
    world = zeros(SIZE, SIZE);
    for k = 1:30
        rands = rand(randSize, randSize) * (rand()*0.5+0.5);
        r = randi(range);
        c = randi(range);
        world((1:randSize)+r, (1:randSize)+c) = rands(1:randSize, 1:randSize);
    end
end

