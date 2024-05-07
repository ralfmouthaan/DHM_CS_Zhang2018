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

% inout signal f 
figure;
imshow(f);


%% Paprameters (1)
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

%% Propagation kernel (2)
Phase=MyMakingPhase(nx,ny,z,lambda,deltaX,deltaY);
figure;imagesc(plotdatacube(angle(Phase)));title('Phase of kernel');axis image;drawnow;
axis off; colormap(hot); colorbar;
E0=ones(nx,ny);  % illumination light
E=MyFieldsPropagation(E0,nx,ny,nz,Phase);  % propagation of illumination light

%% Field measurement and backpropagation (3)
cEs=zeros(nx,ny,nz);
Es=f.*E;
for i=1:nz
    cEs(:,:,i)=fftshift(fft2(Es(:,:,i)));
end
cEsp=sum(cEs.*Phase,3);
S=(ifft2(ifftshift(cEsp)));


f1 = ones(nx,ny);
Es1=f1.*E;
for i=1:nz
    cEs1(:,:,i)=fftshift(fft2(Es1(:,:,i)));
end
cEsp1=sum(cEs1.*Phase,3);
S1=(ifft2(ifftshift(cEsp1)));

% squared field
s=(S+1).*conj(S+1);
s1=(S1+1).*conj(S1+1);
%  diffracted field
g = s./s1; % normalized 
g = im2double(g);
figure;imshow(abs(g),[]);title('Diffracted field')
g=MyC2V(g(:));
transf=MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase);
transf=reshape(MyV2C(transf),nx,ny,nz);
figure;imshow(plotdatacube(abs(transf)),[],'border','tight')

%% Propagation operator (4)
A = @(f_twist) MyForwardOperatorPropagation(f_twist,E,nx,ny,nz,Phase);  % forward propagation operator
AT = @(g) MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase);  % backward propagation operator

%% TwIST algorithm (5)
tau = 0.005; 
piter = 4;
tolA = 1e-6;
iterations = 200;

Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);
Phi = @(f) MyTVphi(f,Nx,Ny,Nz);

[f_reconstruct,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
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
