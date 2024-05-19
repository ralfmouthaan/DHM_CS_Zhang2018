function y = Grad(x)

    [nx, ny] = size(x);
    y= zeros(nx, ny, 2);

    y(:,:,1) = x - circshift(x,[-1 0]);
    y(:,:,2) = x - circshift(x,[0 -1]);

    y(nx,:,1) = 0.0;
    y(:,ny,2) = 0.0;
    y(:,ny/2,1) = 0.0;
    y(:,ny/2,2) = 0.0;

end