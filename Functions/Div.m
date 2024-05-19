function y=Div(F)

% We seem to require that we do x_{n+1} - x_n rather than x_{n} - x{n-1}
% but I do not understand why

x = F(:,:,1);
yx = circshift(x, [1 0]) - x;

y = F(:,:,2);
yy = circshift(y, [0 1]) - y;

y=yx+yy;