%% 4D matrix research
close all;
clear all;
clc;
%% Load Data
% load synthetic data tensor
load('SynthGoUp.mat');
% get 5 different sources (go up 1 meter every time)
source1=zeros(1000,10,3);
source2=zeros(1000,10,3);
source3=zeros(1000,10,3);
source4=zeros(1000,10,3);
source5=zeros(1000,10,3);
source1(:,:,:)=trc(1,:,:,:);
source2(:,:,:)=trc(2,:,:,:);
source3(:,:,:)=trc(3,:,:,:);
source4(:,:,:)=trc(4,:,:,:);
source5(:,:,:)=trc(5,:,:,:);
% visualize 3D

% create 4D tensor out of the data
% ni - i-th source is at position ni of the 4Dtensor
n2=2;
n3=4;
n4=4;
n5=5;
nend=n5;
nochange=false;
nothird=false;
nofourth=true;
nofifth=true;
% which one to pick if no change in time?
source=source2;
tensor4D=zeros(1000,10,3,nend);
if nochange==false
    for i=1:n2-1
        tensor4D(:,:,:,i)=source1;
    end
    for i=n2:n3-1
        tensor4D(:,:,:,i)=source2;
    end
    if nothird==false
        for i=n3:n4-1
            tensor4D(:,:,:,i)=source3;
        end
        if nofourth==false
            for i=n4:n5-1
                tensor4D(:,:,:,i)=source4;
            end
            if nofifth==false
                for i=n5:nend
                    tensor4D(:,:,:,i)=source5;
                end
            else
                for i=n5:nend
                    tensor4D(:,:,:,i)=source4;
                end
            end
        else
            for i=n4:nend
                tensor4D(:,:,:,i)=source3;
            end
        end 
    else
        for i=n3:nend
            tensor4D(:,:,:,i)=source2;
        end
    end
else
    for i=1:nend
        tensor4D(:,:,:,i)=source;
    end
end

save('4Dsyntheticdata.mat','tensor4D');
%% HOSVD
eps = 0.079; %Tolerance

X=tensor(tensor4D);
X_doub=double(X);
%T = hosvd(X,2*sqrt(eps));
T = hosvd(X,2*sqrt(eps),'ranks',[75 10 2 1]);
T_doub=double(T);
fprintf('Time for 10 runs of ST-HOSVD_T:\n');
tic, for i =1:10, T = hosvd(X,2*sqrt(eps),'verbosity',0); end; toc

%Factor matrices (U,V,W,F) and core tensor
core=T.core;
G=double(core);
U=T.U{1};
V=T.U{2};
W=T.U{3};
F=T.U{4};

%% Compression ratio
s_core=size(core);
s_X=size(X);

wh_X=whos('bytes','X');
siz_X=wh_X.bytes*10e-06;
wh_T=whos('bytes','T');
siz_T=wh_T.bytes*10e-06;
CR=(1-(wh_T.bytes/wh_X.bytes))*100;

%% MLSVD
[Ue,Se,sve]=mlsvd(tensor4D,[core.size]);
%T_doub_mlsvd=lmlragen(Ue,Se);

%% Plots
%% Normalized singular values MLSVD
figure(2)

sve{1}=sve{1}/max(sve{1});
sve{2}=sve{2}/max(sve{2});
sve{3}=sve{3}/max(sve{3});
sve{4}=sve{4}/max(sve{4});

for n = 1:4
subplot(1,4,n)
plot(sve{n});
title(sprintf('Mode %d',n));
end
set(gcf,'color','w');
sgtitle('Normalized singular values - MLSVD')

%% Comparison between Recovered and Original HOSVD
k=5;
figure(3)
subplot(1,2,1)
clims = [-0.004 0.004];
imagesc(T_doub(:,:,1,k))
title('Recovered');

subplot(1,2,2)
imagesc(X_doub(:,:,1,k))
title('Original');
set(gcf,'color','w')
%colorbar
sgtitle(sprintf('Reconstruction (X vs Recovered(X))\n Tol: %.4g - CR: %.4g',eps,CR));

%% Comparison between Recovered and Original HOSVD XYZ
figure(4)
j=1;
subplot(1,3,1)
imagesc(T_doub(:,:,j,k))
title('Recovered X');

subplot(1,3,2)
imagesc(T_doub(:,:,j+1,k))
title('Recovered Y');

subplot(1,3,3)
imagesc(T_doub(:,:,j+2,k))
title('Recovered Z');
set(gcf,'color','w')
%% General Structure: Microseismic traces X comp
figure(5)
subplot(1,2,1)
s1=stackedplot(X_doub(:,:,1,k),'Title','Structure of microseismic data');
s1.DisplayLabels = {'Rec 1','Rec 2','Rec 3','Rec 4','Rec 5','Rec 6','Rec 7','Rec 8','Rec 9','Rec 10'};

subplot(1,2,2)
s2=stackedplot(T_doub(:,:,1,k),'Title','Reduced Structure of microseismic data');
s2.DisplayLabels = {'Rec 1','Rec 2','Rec 3','Rec 4','Rec 5','Rec 6','Rec 7','Rec 8','Rec 9','Rec 10'};

set(gcf,'color','w');
set(gcf,'Position',[456.2,163.4,624.8,516.8])

%% General Structure: Microseismic traces -> t1,t2,t3,t4,t5
figure(6)
subplot(1,2,1)
s1=stackedplot(X_doub(:,3,1,:),'Title','Structure of microseismic data (Rec 9)');
s1.DisplayLabels = {'T1','T2','T3','T4','T5'};
subplot(1,2,2)

s2=stackedplot(T_doub(:,3,1,:),'Title','Reduced Structure of microseismic data (Rec 9)');
s2.DisplayLabels = {'T1','T2','T3','T4','T5'};
set(gcf,'color','w')

figure(7)
subplot(1,2,1)
s1=stackedplot(X_doub(:,10,1,:),'Title','Structure of microseismic data (Rec 10)');
s1.DisplayLabels = {'T1','T2','T3','T4','T5'};

subplot(1,2,2)
s2=stackedplot(T_doub(:,10,1,:),'Title','Reduced Structure of microseismic data (Rec 10)');
s2.DisplayLabels = {'T1','T2','T3','T4','T5'};
set(gcf,'color','w')

%% Difference between t1 and t5
rec=5;
dif=(T_doub(:,:,1,5)*-1)-T_doub(:,:,1,1);
dif_or=(X_doub(:,:,1,5)*-1)-X_doub(:,:,1,1);

figure(8)
subplot(1,2,1)
plot(T_doub(:,rec,1,1),'Color',[0 0.4470 0.7410]);
hold on 
plot(T_doub(:,rec,1,5),'Color',[0.8500 0.3250 0.0980]);
title(sprintf('Reduced microseismic trace for T1 and T5 (Rec: %.4g)',rec));
legend('Time 1','Time 5')
ylim([-0.0008 0.0008])

subplot(1,2,2)
plot(dif(:,rec));
title('Reduced microseismic trace for diff[T5-T1]')
legend('Time 5-Time 1')
ylim([-0.0008 0.0008])
set(gcf,'color','w')

figure(9)
subplot(1,2,1)
imagesc(dif(:,:))
title('Recovered [T5-T1]');

subplot(1,2,2)
imagesc(dif_or(:,:))
title('Original [T5-T1]');
set(gcf,'color','w')
sgtitle(sprintf('Reconstruction (X vs Recovered(X))'))

figure(10)
subplot(1,2,1)
s1=stackedplot(dif_or(:,:),'Title','Original microseismic trace for diff[T5-T1]');
s1.DisplayLabels = {'Rec 1','Rec 2','Rec 3','Rec 4','Rec 5','Rec 6','Rec 7','Rec 8','Rec 9','Rec 10'};

subplot(1,2,2)
s2=stackedplot(dif(:,:),'Title','Reduced microseismic trace for diff[T5-T1]');
s2.DisplayLabels = {'Rec 1','Rec 2','Rec 3','Rec 4','Rec 5','Rec 6','Rec 7','Rec 8','Rec 9','Rec 10'};

set(gcf,'color','w');
set(gcf,'Position',[456.2,163.4,624.8,516.8])

figure(11)
subplot(1,2,1)
plot(T_doub(:,rec,1,1),'Color',[0 0.4470 0.7410]);
hold on 
plot((T_doub(:,rec,1,4)),'Color',[0.8500 0.3250 0.0980]);
title(sprintf('Reduced microseismic trace for T1 and T4 (Rec: %.4g)',rec));
legend('Time 1','Time 4')
%ylim([-0.0008 0.0008])

subplot(1,2,2)
plot(T_doub(:,rec,1,1),'Color',[0 0.4470 0.7410]);
hold on 
plot((T_doub(:,rec,1,5)),'Color',[0.8500 0.3250 0.0980]);
title(sprintf('Reduced microseismic trace for T1 and T5 (Rec: %.4g)',rec));
legend('Time 1','Time 5')

set(gcf,'color','w');
%ylim([-0.0008 0.0008])

%% Difference between t1 and t5 (Reverse polarity *-1)
rec=5;
dif=(T_doub(:,:,1,5)*-1)-T_doub(:,:,1,1);

figure(12)
subplot(1,2,1)
plot(T_doub(:,rec,1,1),'Color',[0 0.4470 0.7410]);
hold on 
plot((T_doub(:,rec,1,5)*-1),'Color',[0.8500 0.3250 0.0980]);
title(sprintf('Reduced microseismic trace for T1 and T5 (Rec: %.4g)',rec));
legend('Time 1','Time 5')
ylim([-0.0006 0.0006])

subplot(1,2,2)
plot(dif(:,rec));
title('Reduced microseismic trace for diff[T5-T1]')
legend('Time 5-Time 1')
ylim([-0.0006 0.0006])
set(gcf,'color','w')