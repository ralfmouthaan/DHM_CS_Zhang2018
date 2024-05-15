function y=MyDiv3D(TV)

n=size(TV);

x = TV(:,:,1);
x_shift = circshift(x, [1 0]);
yx = x - x_shift;
yx(1,:) = x(1,:);
yx(n(1),:)= - x_shift(n(1),:);


y = TV(:,:,2);
y_shift = circshift(y, [0 1]);
yy = y - y_shift;
yy(:,1) = y(:,1);
yy(:,n(2)) = -y_shift(:,n(2));

y=yx+yy;
