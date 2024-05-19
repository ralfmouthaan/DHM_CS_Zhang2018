function y=MyDiv3D(TV)

% g = Grad(TV(:,:,1));
% y = -g(:,:,1) - g(:,:,2);

[nx, ny, nz] = size(TV);

x = TV(:,:,1);
x_shift = circshift(x, [1 0]);
yx = x - x_shift;
%yx(1,:) = x(1,:);
yx(nx,:)= -x_shift(nx,:);

y = TV(:,:,2);
y_shift = circshift(y, [0 1]);
yy = y - y_shift;
%yy(:,1) = y(:,1);
yy(:,ny) = -y_shift(:,ny);

y=yx+yy;
