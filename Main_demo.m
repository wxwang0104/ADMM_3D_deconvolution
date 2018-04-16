%% main function for demo
clear all;
close all;
clc;

load('demo_data');
[u1] = ADMM_deconv(im, A);
paratrj = LS_solver(u1, A, Adx ,Ady ,Adz, im, 5);
[y, paras] = ML_threshold(paratrj);
paratrj = LS_solver2(paras(y==1,:), A, Adx, Ady, Adz, im, 10);
recov_pos = paratrj(:,end-2:end);

[v,h] = size(im);
figure
imagesc(2:3:3*v, 2:3:3*h, im)
colormap gray
hold on
plot(real_pos(:,1),real_pos(:,2),'go','markersize',8,'linewidth',2)
plot(recov_pos(:,1),recov_pos(:,2),'b+','markersize',8,'linewidth',2)
legend('true positions','recovered')
hold off