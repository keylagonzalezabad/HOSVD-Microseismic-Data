%% HOSVD 
close all;
clear all;
clc;
%% Load Data
%Load data tensor
str=load('tensor_test.mat','tensor');
str_xyz=load('xyztensor.mat','tensor_xyz');

%Convert to tensor 
clear tensor
a=str.tensor; %Microseismic Data
b=str_xyz.tensor_xyz; %XYZ Microseismic Data
X=tensor(b);
X_doub=double(X);
xs=size(X_doub);

eps = 0.0002; %Tolerance
%% HOSVD estimation
T = hosvd(X,2*sqrt(eps));
T_doub=double(T);
fprintf('Time for 10 runs of ST-HOSVD_T:\n');
tic, for i =1:10, T = hosvd(X,2*sqrt(eps),'verbosity',0); end; toc

%Factor matrices (U,V,W) and core tensor
core=T.core;
G=double(core);
U=T.U{1};
V=T.U{2};
W=T.U{3};

%% Compression ratio
s_core=size(core);
s_X=size(X);

wh_X=whos('bytes','X');
siz_X=wh_X.bytes*10e-06;
wh_T=whos('bytes','T');
siz_T=wh_T.bytes*10e-06;
CR=(1-(wh_T.bytes/wh_X.bytes))*100;

%% Singular values (Frobenius Norm)
sv_1 = zeros(size(G,1),1);
for i = 1:size(G,1)
    sv_1(i)=sqrt(sum(sum((G(i,:,:).^2))));
end
sv_1=sv_1';
sv_1=(sv_1./max(sv_1));

sv_2 = zeros(size(G,2),1);
for i = 1:size(G,2)
    sv_2(i)=sqrt(sum(sum((G(:,i,:).^2))));
end
sv_2=sv_2';
sv_2=(sv_2./max(sv_2));

sv_3 = zeros(size(G,3),1);
for i = 1:size(G,3)
    sv_3(i)=sqrt(sum(sum((G(:,:,i).^2))));
end
sv_3=sv_3';
sv_3=(sv_3./max(sv_3));

sv = {sv_1 sv_2 sv_3};

%% Plots
%% Normalized singular values 
figure(1)
for n = 1:3
subplot(1,3,n)
plot(sv{n});
title(sprintf('Mode %d',n));
end
set(gcf,'color','w');
sgtitle('Normalized singular values')

%% Comparison between Recovered and Original
k=20;

figure(2)
subplot(1,2,1)
imagesc(T_doub(:,:,k))
title('Recovered');

subplot(1,2,2)
imagesc(X_doub(:,:,k))
title('Original');
set(gcf,'color','w')
sgtitle(sprintf('Reconstruction (X vs Recovered(X))\n Tol: %.4g - CR: %.4g',eps,CR));

%% Number of variables
figure(3)
subplot(1,2,1)
V=100-CR;
pi=[CR V];
pie(pi);
title('Variables');

subplot(1,2,2)
r=(s_X(1)*s_X(2)*s_X(3))*(1-(CR/100));
p=(s_X(1)*s_X(2)*s_X(3));
bar(1, p, 'FaceColor',[0 0.4470 0.7410])
hold on
bar(2, r, 'FaceColor',[0.9290 0.6940 0.1250])
legend('Tensor','HOSVD')
title('Number of variables');
set(gca, 'xticklabels', []);
set(gcf,'color','w');

sgtitle('Number of variables')
