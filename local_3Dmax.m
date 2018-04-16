function [center, dx, dy, dz] = local_3Dmax(u1)
% clustered non-zero elements are considered as identified particles 
% and 3D averaged positions are calculated
% center, dx, dy, dz are all 3D matrix with the same size as u1
% center records center position of each particle
% dx, dy, dz records position corrections in x, y, and z dimensions.
[a,b,c] = size(u1);
[I, q] = max(u1(:));
[X, Y, Z] = meshgrid(1:b, 1:a, 1:c);
center = zeros(a,b,c);
dx = center;
dy = center;
dz = center;
while I>0
    [yc, xc, zc] = ind2sub([a,b,c], q);
    inx3 = zeros(a,b,c);
    inx3(q) = 1;
    inx3 = label_neibs(u1, inx3, yc, xc, zc, a, b, c);
    inx3 = inx3>0;
    photons = sum(u1(inx3));
    elx = sum(X(inx3).*u1(inx3))/photons;
    ely = sum(Y(inx3).*u1(inx3))/photons;
    elz = sum(Z(inx3).*u1(inx3))/photons;
    center(round(ely),round(elx),round(elz)) = photons;
    dx(round(ely),round(elx),round(elz)) = elx-round(elx);
    dy(round(ely),round(elx),round(elz)) = ely-round(ely);
    dz(round(ely),round(elx),round(elz)) = elz-round(elz);
    u1(inx3) = 0;
    [I, q] = max(u1(:));
end
thread = max(center(:))*0.05;
dx(center<=thread) = 0;
dy(center<=thread) = 0;
dz(center<=thread) = 0;
center(center<=thread) = 0;
end

function inx3 = label_neibs(h2d, inx3, yc, xc, zc, a, b, c)
[dx, dy, dz] = meshgrid(xc-1:xc+1,yc-1:yc+1, zc-1:zc+1);
dx(14) = []; dy(14) = []; dz(14) = [];
inx = (dx>0&dx<=b)&(dy>0&dy<=a)&(dz>0&dz<=c);
dx = dx(inx); dy = dy(inx); dz = dz(inx);

for k = 1:numel(dy)
    if(inx3(dy(k),dx(k),dz(k))==0 && h2d(dy(k),dx(k),dz(k))>0 && h2d(dy(k),dx(k),dz(k))<=h2d(yc,xc,zc))
        inx3(dy(k),dx(k),dz(k)) = 1;
        inx3 = label_neibs(h2d, inx3, dy(k),dx(k),dz(k),a,b,c);
    end
end
end