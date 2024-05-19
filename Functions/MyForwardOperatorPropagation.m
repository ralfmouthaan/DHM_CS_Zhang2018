function S=MyForwardOperatorPropagation(S,E,Nx,Ny,phase)

S=reshape(MyV2C(S),Nx,Ny);

%S = S.*E;
S=fftshift(fft2(S));
S = S.*phase;
S=ifft2(ifftshift(S));
S = abs(S);

S=MyC2V(S(:));
