function [u1] = ADMM_deconv(img, psf)
%% 3D deconvolution using ADMM based algorithm
% this ADMM algorithm is based on two publications:
% 1,Almeida, M. S. C. & Figueiredo, M. A. T. Deconvolving images with unknown
% boundaries using the alternating direction method of multipliers.
% IEEE Trans. Image Process. 22, 3074-3086 (2013).
% 2,Matakos, A., Ramani, S. & Fessler, J. A. Accelerated edge-preserving
% image restoration without boundary artifacts. IEEE Trans. Image Process. 22, 2019-2029 (2013).
% 
% psf should be a 3D matrix recording the 3D PSF
% img should be a 2D image from 3D measurement or simulation
% in defaul setting, # of img pixels in each dimension*3 = # of psf pixels in the same dimension
% Variables are defined according to our paper: Generalized recovery algorithm for 3D super-resolution microscopy
%% initial variable setting
[ny, nx, nz] = size(psf);
up_res = 3;% each pixel in image is break into 3X3 pixels, meaning dhpsf should be 3X3 larger than image
im_temp = zeros(ny, nx);
inx1 = round(up_res/2):up_res:ny;
inx2 = round(up_res/2):up_res:nx;
im_temp(inx1, inx2) = img;
y = im_temp(:)';
T = zeros(1, nz); T(nz) = 1;
A = psf; % use A matrix to represent the 3D PSF, fA is fast Fourier transform of A
% u0 is the 3D image space, x is the 3D recovered space, u1 is the sparse version of x
% fu0, fx, fu1 are the corresponding Fourier transform.
% Because most of the calculation is done in Fourier domain, corresponding
% variable in space domain may not be required.
fx = zeros(size(psf));
feta0 = fx; feta1 = fx;
lambda = 0.005; % lambda = 0.005 was very good
mu = min(1,50000*lambda);% those three parameters may difficult to decide, 0.5 was good
nu = mu/10/lambda;
N = 1000;% 1000 iterations
TT = (T'*T+mu*eye(nz));
fA = fftn(A);
invall = abs(fA).^2+nu;
thd = max(img(:))/max(psf(:))*0.1/5;% thresdhold to calculate u1 from x
u0_2d = TT\(T'*y);
u0 = reshape(u0_2d', ny, nx, nz);
fu0 = fftn(u0);
fu1 = zeros(size(u0));
%% ADMM iterations
% Calculations are mostly conducted in Fourier domain
for k = 1:N
    ftemp_1 = conj(fA).*(fu0-feta0);
    fx = (ftemp_1+nu*(fu1-feta1))./invall;% update x, equation (7)
    fAx = fA.*fx;
    feta0 = feta0-(fu0-fAx);% update eta0, equation (8)
    feta1 = feta1-(fu1-fx);% update eta1, equation (9)
    Ax_eta_2d = ifftn(fAx+feta0);
    if nz > 1
        Ax_eta_2d = reshape(XYZ_rot(Ax_eta_2d, [3,1,2]), nz, ny*nx);
    else
        Ax_eta_2d = Ax_eta_2d(:)';
    end
    u0_2d = TT\(T'*y+mu*Ax_eta_2d);
    u0 = reshape(u0_2d', ny, nx, nz);% update u0, equation (5)
    fu0 = fftn(u0);
    x_eta1 = ifftn(fx+feta1);
    u1 = max(x_eta1-thd, 0);% update u1, equation (6)
    fu1 = fftn(u1);
end
%% Get final output
u1 = fftshift(u1);
u1 = ifftshift(u1,3);