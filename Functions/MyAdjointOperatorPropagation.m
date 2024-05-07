function eta=MyAdjointOperatorPropagation(S,E,Nx,Ny,Nz,phase)

S=reshape(MyV2C(S),Nx,Ny);

cEsp=ifftshift(ifft2(real(S)));
cEs=conj(phase).*cEsp;
eta=fft2(ifftshift(cEs));
eta=conj(E).*eta;

eta=MyC2V(eta(:));
