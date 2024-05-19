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
Nx=nx; Ny=ny*2;
lambda=0.532;  % wavelength (um)
z=12000;  % distance from detector to first reconstructed plane (um)
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
g = im2double(g);
figure;
imshow(abs(g));
title('Diffracted field')

g=MyC2V(g(:));

%% Function Definitions

tau = 0.005; 
tolA = 1e-6;
iterations = 200;

A = @(fimg) MyForwardOperatorPropagation(fimg,E,nx,ny,Phase);  % forward propagation operator
AT = @(fimg) MyAdjointOperatorPropagation(fimg,E,nx,ny,Phase);  % backward propagation operator
Psi = @(f,lambda) Denoise(f,lambda, Nx, Ny);
Phi = @(f) TV(f, Nx, Ny);

%% TwIST algorithm

f_reconstruct = ...
    TwIST(g, A, tau, 'AT', AT, 'Psi', Psi, 'Phi', Phi, ...
    'Initialization',2, 'Monotone', 1, 'StopCriterion',1,...
    'MaxIterA', iterations, 'MinIterA', iterations, ...
    'ToleranceA', tolA, 'Verbose', 1);

f_reconstruct = reshape(MyV2C(f_reconstruct), nx, ny);

figure;
imshow(abs(f_reconstruct))

fprintf('Err = %0.5f\n', sum(sum(abs(f - f_reconstruct)))/nx/ny);

function Kernel = GenerateKernel_ASM(k0, kx, z)
    
    term = k0^2 - kx.^2 - kx.'.^2;
    term(term < 0) = 0;
    Kernel = exp(1i*2*pi*z*sqrt(term));

end
function x = V2C(x)

    nx = length(x);
    x = x(1:nx/2) + 1i*x(nx/2+1:nx);
    nx = sqrt(length(x));
    x = reshape(x, nx, nx);

end
function x = C2V(x)

    x = x(:);
    x=[real(x); imag(x)];

end
function y = TV(x, nx, ny)

    x=reshape(x,nx,ny);
    TV = Grad(x);
    y=sum(abs(TV(:)));

end
