function y=MyTVphi(x,Nvx,Nvy)

x=reshape(x,Nvx,Nvy);

[nx,ny]=size(x);
TV=zeros(nx,ny,3);

TV(:,:,1)=circshift(x,[-1 0 0])-x;
TV(nx,:,1)=0.0;

TV(:,:,2)=circshift(x,[0 -1 0])-x;
TV(:,ny,2)=0.0;

y=sum(abs(TV(:)));

end