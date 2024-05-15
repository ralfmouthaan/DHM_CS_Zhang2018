function TV = MyTV3D_conv(x)

[nx,ny]=size(x);
TV=zeros(nx,ny,3);
TV(:,:,1)=circshift(x,[-1 0])-x;
TV(nx,:,1)=0.0;

TV(:,:,2)=circshift(x,[0 -1])-x;
TV(:,ny,2)=0.0;
