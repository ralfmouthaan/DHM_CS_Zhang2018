function p=MyProjectionTV(g,tau,lam,iter)

% This is a denoising algorithm. I do not know where it comes from or why
% it works, but it does seem to work.
% Something to do with Rudin–Osher–Fatemi?

[nx,ny] = size(g);
pn = zeros(nx,ny,2);
div_pn = zeros(nx,ny);
b = zeros(nx, ny);

for i=1:iter

    a = -Grad(div_pn-g./lam);

    b(:,:,1) = sqrt(a(:,:,1).^2 + a(:,:,2).^2);
    b(:,:,2) = b(:,:,1);
    pn = (pn+tau.*a)./(1.0+tau.*b);

    div_pn = -Div(pn);
    
end

p = -lam.*Div(pn);

