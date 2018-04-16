function [y, paras] = ML_threshold(paratrj)
%% Use learned threshold to get rid of overfittings
theta = importdata('Theta_of_LS.mat');
I = paratrj(:,21);
paratrj(I==0,:) = [];
paras = paratrj(:,1:4);
for i = 1:4:20
    paratrj(:,i+1:i+3) = abs(paratrj(:,i+1:i+3)-paratrj(:,i+5:i+7));
    I = paratrj(:,i);
    paratrj(:,i) = paratrj(:,i)/mean(I);
end
I = paratrj(:,21);
paratrj(:,21) = paratrj(:,21)/mean(I);
paratrj(:,22:end) = [];
n = size(paratrj,1);
paratrj = [ones(n,1), paratrj];
y = 1./(1+exp(-paratrj*theta))>0.5;