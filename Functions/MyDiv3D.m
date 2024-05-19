function y=MyDiv3D(TV)

x = TV(:,:,1);
yx = circshift(x, [1 0]) - x;

y = TV(:,:,2);
yy = circshift(y, [0 1]) - y;

y=yx+yy;