% Ralf Mouthaan
% University of Adelaide
% May 2024
% 
% Fork of Zhang's Compressive DHM code.
% Original paper: Zhang et al. "Twin-Image-Free Holography: A Compressive
% Sensing Approach", PRL 2018.
% 
% I'm aiming to refactor Zhang's code to understand how it works.
 
close all;
clear variables;
clc;

addpath('./Functions');

f = imread('cell.jpg');
f = rgb2gray(f);
f = imresize(f,[500,500]);
f = 1-im2double(f);

figure;
imshow(f);

%% Parameters (1)
nx=size(f,1);  % data size
ny=size(f,2);
nz=1;
lambda=0.532;  % wavelength (um)
k = 2*pi/lambda;
detector_size=4;  % pixel pitch (um)
sensor_size=nx*detector_size;  % detector size (um)
z=12000;  % distance from detector to first reconstructed plane (um)
deltaX=detector_size;
deltaY=detector_size;
Nx=nx;
Ny=ny*nz*2;
Nz=1;

dx = deltaX;
k0 = 1/lambda; % I think this should be 2pi/lambda
x = (-nx/2:nx/2-1)*dx;
kx = (-nx/2:nx/2-1)/(nx*dx);

%% Propagation kernel (2)

Phase = GenerateKernel_ASM(k0, kx, z);

figure;
imagesc(angle(Phase));
title('Phase of kernel');
axis image;
colormap(hot); 
colorbar;
drawnow;

 % Propagation of illumination
E=ones(nx,ny); 
E=fftshift(fft2(E));
E=E.*conj(Phase);
E=ifft2(ifftshift(E));

%% Field measurement and backpropagation (3)

% Propagation of object field
S=f.*E;
S=fftshift(fft2(S));
S=S.*Phase;
S=ifft2(ifftshift(S));
S=(S+1).*conj(S+1);

% Propagation of reference field
S1 = ones(nx,ny);
S1=S1.*E;
S1=fftshift(fft2(S1));
S1=S1.*Phase;
S1=ifft2(ifftshift(S1));
S1=(S1+1).*conj(S1+1);

% Diffracted field
g = S./S1; % normalized 
g = im2double(g);
figure;
imshow(abs(g));
title('Diffracted field')

g=MyC2V(g(:));

%% Propagation operator (4)
A = @(fimg) MyForwardOperatorPropagation(fimg,E,nx,ny,nz,Phase);  % forward propagation operator
AT = @(fimg) MyAdjointOperatorPropagation(fimg,E,nx,ny,nz,Phase);  % backward propagation operator

%% TwIST algorithm (5)
tau = 0.005; 
piter = 4;
tolA = 1e-6;
iterations = 200;

Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);
Phi = @(f) MyTVphi(f,Nx,Ny,Nz);

[f_reconstruct,~,obj_twist,...
    times_twist,~,mse_twist]= ...
    TwIST(g,A,tau,...
    'AT', AT, ...
    'Psi', Psi, ...
    'Phi',Phi, ...
    'Initialization',2,...
    'Monotone',1,...
    'StopCriterion',1,...
    'MaxIterA',iterations,...
    'MinIterA',iterations,...
    'ToleranceA',tolA,...
    'Verbose', 1);


f_reconstruct=reshape(MyV2C(f_reconstruct),nx,ny,nz);
re=(abs(f_reconstruct));
re = im2double(re);
figure,imshow(re,[])

fprintf('Err = %0.5f\n', sum(sum(abs(f - re)))/nx/ny);

function Kernel = GenerateKernel_ASM(k0, kx, z)
    
    term = k0^2 - kx.^2 - kx.'.^2;
    term(term < 0) = 0;
    Kernel = exp(1i*2*pi*z*sqrt(term));

end
