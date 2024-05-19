function p=MyProjectionTV(g,tau,lam,iter)

[nx,ny]=size(g);
pn=zeros(nx,ny,2);
div_pn=zeros(nx,ny);
b=pn;

for i=1:iter

    a = Grad(div_pn-g./lam);

    b(:,:,1)=sqrt(a(:,:,1).^2+a(:,:,2).^2);
    b(:,:,2)=b(:,:,1);
    pn=(pn+tau.*a)./(1.0+tau.*b);

    div_pn=MyDiv3D(pn);
    
end

p=lam.*MyDiv3D(pn);

