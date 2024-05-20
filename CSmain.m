% Ralf Mouthaan
% University of Adelaide
% May 2024
% 
% Fork of Zhang's Compressive DHM code.
% Original paper: Zhang et al. "Twin-Image-Free Holography: A Compressive
% Sensing Approach", PRL 2018.
% 
% I'm aiming to refactor Zhang's code to understand how it works.

%%
 
close all; clear variables; clc;
addpath('./Functions');

%% User-Defined Parameters

nx = 500; ny = 500;
lambda = 0.532;  % wavelength (um)
z = 12000;  % distance from detector to first reconstructed plane (um)
dx = 4; % um

k0 = 1/lambda;
x = (-nx/2:nx/2-1)*dx;
kx = (-nx/2:nx/2-1)/(nx*dx);

%% Read in image

f = imread('cell.jpg');
f = rgb2gray(f);
f = imresize(f, [nx, ny]);
f = 1 - im2double(f);

figure;
imshow(f);
drawnow;

%% Propagation kernel

Phase = GenerateKernel_ASM(k0, kx, z);

figure;
imagesc(angle(Phase));
title('Phase of kernel');
axis image;
colormap(hot); 
colorbar;
drawnow;

%% Field measurement and backpropagation

 % Backpropagation of illumination
E=ones(nx,ny); 
E=fftshift(fft2(E));
E=E.*conj(Phase);
E=ifft2(ifftshift(E));

% Propagation of object field
S=f.*E;
S=fftshift(fft2(S));
S=S.*Phase;
S=ifft2(ifftshift(S));

% Diffracted field
g = abs(S);
figure;
imshow(abs(g));
title('Diffracted field')

g = C2V(g, nx, ny);

%% Function Definitions

tau = 0.005; 
tolA = 1e-6;
iterations = 200;

A = @(fimg) ForwardPropagation(fimg,E,nx,ny,Phase);  % forward propagation operator
AT = @(fimg) BackwardPropagation(fimg,E,nx,ny,Phase);  % backward propagation operator
Psi = @(f, lambda) Denoise(f, lambda, nx, ny);
Phi = @(f) TotalVariance(f, nx, ny);

%% TwIST algorithm

f_reconstruct = ...
    TwIST(g, A, tau, 'AT', AT, 'Psi', Psi, 'Phi', Phi, ...
    'Initialization',2, 'Monotone', 1, 'StopCriterion',1,...
    'MaxIterA', iterations, 'MinIterA', iterations, ...
    'ToleranceA', tolA, 'Verbose', 1);

f_reconstruct = V2C(f_reconstruct, nx, ny);

figure;
imshow(abs(f_reconstruct))

fprintf('Err = %0.5f\n', sum(sum(abs(abs(f) - abs(f_reconstruct))))/nx/ny);

%% Functions

function Kernel = GenerateKernel_ASM(k0, kx, z)
    
    term = k0^2 - kx.^2 - kx.'.^2;
    term(term < 0) = 0;
    Kernel = exp(1i*2*pi*z*sqrt(term));

end
function F = V2C(F, nx, ny)

    assert(length(F) == nx*ny*2);

    F = F(1:nx*ny) + 1i*F(nx*ny+1:nx*ny*2);
    F = reshape(F, nx, ny);

end
function F = C2V(F, nx, ny)

    assert(size(F, 1) == nx);
    assert(size(F, 2) == ny);

    F = F(:);
    F=[real(F); imag(F)];

end
function F = ForwardPropagation(F, E, Nx, Ny, Kernel)

    F = V2C(F, Nx, Ny);
    
    %S = S.*E;
    F = fftshift(fft2(F));
    F = F.*Kernel;
    F = ifft2(ifftshift(F));
    F = abs(F);
    
    F = C2V(F, Nx, Ny);

end
function F = BackwardPropagation(F, E, Nx, Ny, Kernel)

    F = V2C(F, Nx, Ny);
    
    F = ifftshift(ifft2(F));
    F = conj(Kernel).*F;
    F = fft2(ifftshift(F));
    F = conj(E).*F;
    
    F = C2V(F, Nx, Ny);

end
function RetVal = TotalVariance(F, nx, ny)

    F=reshape(F,nx,2*ny);
    RetVal = Grad(F);
    RetVal=sum(abs(RetVal(:)));

end
function F = Denoise(F, lambda, nx, ny)

    % This is a denoising algorithm. I do not know where it comes from or why
    % it works, but it does seem to work. Something to do with Rudin–Osher–Fatemi?
    
    F = reshape(F,nx,2*ny);
    lambda = lambda*0.5;
    tau = 0.05; % Note, this tau is different in value to the one defined in the main script.
                 % But, I don't know if they have the same significance.
    
    pn = zeros(nx,2*ny,2);
    div_pn = zeros(nx,2*ny);
    b = zeros(nx, 2*ny);
    
    for i=1:4
    
        a = -Grad(div_pn - F./lambda);
        b(:,:,1) = sqrt(a(:,:,1).^2 + a(:,:,2).^2);
        b(:,:,2) = b(:,:,1);
        pn = (pn+tau.*a)./(1.0+tau.*b);
        div_pn = -Div(pn);
        
    end
    
    F = F - lambda.*div_pn;
    
    F = reshape(F,2*nx*ny,1);

end
function RetVal = Grad(F)

    [nx, ny] = size(F);
    RetVal= zeros(nx, ny, 2);

    RetVal(:,:,1) = F - circshift(F,[-1 0]);
    RetVal(:,:,2) = F - circshift(F,[0 -1]);

    RetVal(nx,:,1) = 0.0;
    RetVal(:,ny,2) = 0.0;
    RetVal(:,ny/2,1) = 0.0;
    RetVal(:,ny/2,2) = 0.0;

end
function RetVal = Div(F)

    % We seem to require that we do x_{n+1} - x_n rather than x_{n} - x{n-1}
    % but I do not understand why
    
    Fx = F(:,:,1);
    Fx = circshift(Fx, [1 0]) - Fx;
    
    Fy = F(:,:,2);
    Fy = circshift(Fy, [0 1]) - Fy;

    RetVal = Fx + Fy;

end