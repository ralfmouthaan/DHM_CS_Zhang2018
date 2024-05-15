function y = Grad(x)

    [nx, ny] = size(x);
    y = zeros(nx, ny, 3);

    y(:,:,1)=circshift(x,[-1 0 0])-x;
    y(nx,:,1)=0.0;
    
    y(:,:,2)=circshift(x,[0 -1 0])-x;
    y(:,ny,2)=0.0;

end