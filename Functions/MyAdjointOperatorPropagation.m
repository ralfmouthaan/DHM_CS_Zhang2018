function S=MyAdjointOperatorPropagation(S,E,Nx,Ny,phase)

S=reshape(MyV2C(S),Nx,Ny);

S=ifftshift(ifft2(S));
S=conj(phase).*S;
S=fft2(ifftshift(S));
S=conj(E).*S;

S=MyC2V(S(:));
