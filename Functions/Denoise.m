function y=Denoise(F,lambda,tau,nx,ny,nz)

% This is a denoising algorithm. I do not know where it comes from or why
% it works, but it does seem to work. Something to do with Rudin–Osher–Fatemi?

F = reshape(F,nx,ny,nz);
lambda = lambda*0.5;

[nx,ny] = size(F);
pn = zeros(nx,ny,2);
div_pn = zeros(nx,ny);
b = zeros(nx, ny);

for i=1:4

    a = -Grad(div_pn-F./lambda);
    b(:,:,1) = sqrt(a(:,:,1).^2 + a(:,:,2).^2);
    b(:,:,2) = b(:,:,1);
    pn = (pn+tau.*a)./(1.0+tau.*b);
    div_pn = -Div(pn);
    
end

F = F - lambda.*div_pn;

y = reshape(F,nx*ny*nz,1);