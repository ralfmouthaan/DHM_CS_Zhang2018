function y = Grad(x)

    [nx, ny] = size(x);
    y= zeros(nx, ny, 2);

    y(:,:,1)=circshift(x,[-1 0])-x;
    y(:,:,2)=circshift(x,[0 -1])-x;

    y(nx,:,1) = 0.0;
    y(:,ny,2) = 0.0;
    y(:,ny/2,1) = 0.0;
    y(:,ny/2,2) = 0.0;

end