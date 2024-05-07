function S=MyForwardOperatorPropagation(S,E,Nx,Ny,Nz,phase)

S=reshape(MyV2C(S),Nx,Ny);

S=S.*E;
S=fftshift(fft2(S));
S = S.*phase;
S=real((ifft2(ifftshift(S))));

S=MyC2V(S(:));
