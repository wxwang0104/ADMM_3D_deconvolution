function [paratrj] = LS_solver2(paras, A, Adx, Ady, Adz, im, max_it)
%% repeat the LS fitting use determined 'true' particles
[a,b,c] = size(A);
I = paras(:,1);
ix = round(paras(:,2));
iy = round(paras(:,3));
iz = round(paras(:,4));
n = numel(I);
paratrj = zeros(n, 4*max_it+4);
paratrj(:,1:4) = paras;
m = numel(im(:));
up = 3;
img = zeros(size(A(:,:,1)));
img(ceil(up/2):up:end, ceil(up/2):up:end) = im;
X = zeros(m*up*up,4*n);%iy2(i):3:end,ix2(i):3:end
for i = 1:numel(I)
    ta = zeros(a,b);
    ta_x = ta;
    ta_y = ta;
    ta_z = ta;
    ty1 = max(iy(i)-a/2,1);
    ty2 = min(iy(i)-a/2-1+a,a);
    tx1 = max(ix(i)-b/2,1);
    tx2 = min(ix(i)-b/2-1+b,b);
    py1 = max(a/2+1-iy(i)+1,1);
    py2 = min(a/2+1-iy(i)+a,a);
    px1 = max(b/2+1-ix(i)+1,1);
    px2 = min(b/2+1-ix(i)+b,b);
    ta(ty1:ty2,tx1:tx2) = A(py1:py2,px1:px2,end+1-iz(i));
    X(:,i) = ta(:);
    
    ta_x(ty1:ty2,tx1:tx2) = -Adx(py1:py2,px1:px2,end+1-iz(i));
    X(:,n+i) = ta_x(:);

    ta_y(ty1:ty2,tx1:tx2) = -Ady(py1:py2,px1:px2,end+1-iz(i));
    X(:,2*n+i) = ta_y(:);

    ta_z(ty1:ty2,tx1:tx2) = -Adz(py1:py2,px1:px2,end+1-iz(i));
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
%% 
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
        ta = zeros(a,b);
        ta_x = ta;
        ta_y = ta;
        ta_z = ta;
        ty1 = max(iy(i)-a/2,1);
        ty2 = min(iy(i)-a/2-1+a,a);
        tx1 = max(ix(i)-b/2,1);
        tx2 = min(ix(i)-b/2-1+b,b);
        py1 = max(a/2+1-iy(i)+1,1);
        py2 = min(a/2+1-iy(i)+a,a);
        px1 = max(b/2+1-ix(i)+1,1);
        px2 = min(b/2+1-ix(i)+b,b);
        ta(ty1:ty2,tx1:tx2) = A(py1:py2,px1:px2,end+1-iz(i));
        X(:,i) = ta(:);
        %
        ta_x(ty1:ty2,tx1:tx2) = -Adx(py1:py2,px1:px2,end+1-iz(i));
        X(:,n+i) = ta_x(:);

        ta_y(ty1:ty2,tx1:tx2) = -Ady(py1:py2,px1:px2,end+1-iz(i));
        X(:,2*n+i) = ta_y(:);

        ta_z(ty1:ty2,tx1:tx2) = -Adz(py1:py2,px1:px2,end+1-iz(i));
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