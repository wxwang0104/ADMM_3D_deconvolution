function [paratrj] = LS_solver(u1, A, Adx, Ady, Adz, im, max_it)
%% Least squares fitting
% Adx, Ady, Adz are first order corrections of the 3D PSF in x, y and z.
% im is the original 2D image
% max_it specifies the max number of iterations
[center, elx, ely, elz] = local_3Dmax(u1);
[a,b,c] = size(center);
elx = elx(center>0);
ely = ely(center>0);
elz = elz(center>0);
I = center(center>0);
q = find(center>0);
[iy,ix,iz] = ind2sub(size(center), q);
n = numel(I);
paratrj = zeros(n, 4*max_it+4);
paratrj(:,1:4) = [I, elx+ix, ely+iy, elz+iz];
m = numel(im(:));
up = 3;
img = zeros(size(center(:,:,1)));
img(ceil(up/2):up:end, ceil(up/2):up:end) = im;
X = zeros(m*up*up,4*n);%iy2(i):3:end,ix2(i):3:end
for i = 1:numel(I)
    % 
    ta = A(:,:,end+1-iz(i));
    point = zeros(size(ta));
    point(iy(i),ix(i)) = 1;
    ta = fftshift(ifft2(fft2(ta).*fft2(point)));
    X(:,i) = ta(:);
    % x direction correction
    ta_x = -Adx(:,:,end+1-iz(i));
    ta_x = fftshift(ifft2(fft2(ta_x).*fft2(point)));
    X(:,n+i) = ta_x(:);
    % y direction correction
    ta_y = -Ady(:,:,end+1-iz(i));
    ta_y = fftshift(ifft2(fft2(ta_y).*fft2(point)));
    X(:,2*n+i) = ta_y(:);
    % z direction correction
    ta_z = -Adz(:,:,end+1-iz(i));
    ta_z = fftshift(ifft2(fft2(ta_z).*fft2(point)));
    X(:,3*n+i) = ta_z(:);
end
paras = (X'*X)\(X'*img(:));
paras = reshape(paras, n, 4);
paras(:,2) = paras(:,2)./paras(:,1);
paras(:,3) = paras(:,3)./paras(:,1);
paras(:,4) = paras(:,4)./paras(:,1);
paratrj(:,5:8) = [paras(:,1), paras(:,2)+ix, paras(:,3)+iy, paras(:,4)+iz];
paras = paratrj(:,5:8);
ix = round(paras(:,2));
iy = round(paras(:,3));
iz = round(paras(:,4));
%% repeat the least squares fitting for max_it times
all_n = n;
all_inx = 1:all_n;
I = paras(:,1);
ratio = 0.05;
thread = max(I)*ratio;
for s = 2:max_it
    inx1 = all_inx(iz>c|iz<1|ix>b|ix<1|iy>a|iy<1 | paras(:,1)<=thread);
    all_inx(inx1) = 0;
    paras(inx1, :) = [];
    
    ix = round(paras(:,2));
    iy = round(paras(:,3));
    iz = round(paras(:,4));
    
    n = numel(iy);
    X = zeros(m*up*up,4*n);
    for i = 1:n
        ta = A(:,:,end+1-iz(i));
        point = zeros(size(ta));
        point(iy(i),ix(i)) = 1;
        ta = fftshift(ifft2(fft2(ta).*fft2(point)));
        X(:,i) = ta(:);
        %
        ta_x = -Adx(:,:,end+1-iz(i));
        ta_x = fftshift(ifft2(fft2(ta_x).*fft2(point)));
        X(:,n+i) = ta_x(:);

        ta_y = -Ady(:,:,end+1-iz(i));
        ta_y = fftshift(ifft2(fft2(ta_y).*fft2(point)));
        X(:,2*n+i) = ta_y(:);

        ta_z = -Adz(:,:,end+1-iz(i));
        ta_z = fftshift(ifft2(fft2(ta_z).*fft2(point)));
        X(:,3*n+i) = ta_z(:);
    end
    paras = (X'*X)\(X'*img(:));
    paras = reshape(paras, n, 4);
    paras(:,2) = paras(:,2)./paras(:,1);
    paras(:,3) = paras(:,3)./paras(:,1);
    paras(:,4) = paras(:,4)./paras(:,1);
    
    paratrj(all_inx(all_inx>0), s*4+1:s*4+4)=[paras(:,1),paras(:,2)+ix,paras(:,3)+iy,paras(:,4)+iz];
    paras = paratrj(:,s*4+1:s*4+4);
    I = paras(:,1);
    
    thread = max(I)*ratio;
    ix = round(paras(:,2));
    iy = round(paras(:,3));
    iz = round(paras(:,4));
    all_inx = 1:all_n;
end