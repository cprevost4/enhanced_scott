clear all
close all
clc
rng(107)

%% Load data

load('Isabella_lake_preproc_subim1')
 clear MSI
SRI = HSI; clear HSI; %SRI = SRI(1:96,1:96,:);
P3 = SRF; clear SRF; 
d1 = 2; d2 = 2; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
MSI = tmprod(SRI,P3,3);

SNRh = 30; SNRm = 30;
HSI = awgn(HSI,SNRh,'measured');
MSI = awgn(MSI,SNRm,'measured');

sigma_h = 10^(-SNRh/10); sigma_m = 10^(-SNRm/10);
opts.lambda = (sigma_h^2)./(sigma_m^2);

%% SCOTT, for comparison

R = [40 40 6];
tic;
[SRI11, info] = scott(HSI, MSI, P1, P2, P3, R);
t11 = toc;
err11 = compute_metrics(SRI,SRI11,d1,d2)

%% ESCOTT

%Rmax = floor([size(HSI,1)/opts.Nblocks(1) size(HSI,2)/opts.Nblocks(2) size(MSI,3)]);

opts.Nblocks = [4 4];
R = [11 11 3];
tic;
SRI31 = escott(HSI, MSI, P1, P2, P3, R, opts);
t31 = toc;
err31 = compute_metrics(SRI,SRI31,d1,d2)

opts.Nblocks = [2 2];
R = [22 22 4];
tic;
SRI32 = escott(HSI, MSI, P1, P2, P3, R, opts);
t32 = toc;
err32 = compute_metrics(SRI,SRI32,d1,d2)

%% Figures

[W, ~, ~] = svds(tens2mat(SRI,3,[]),3);

SRI = tmprod(SRI,W',3);
rgb(:,:,1) = imadjust(rescale(SRI(:,:,1),1));
rgb(:,:,2) = imadjust(rescale(SRI(:,:,2),1));
rgb(:,:,3) = imadjust(rescale(SRI(:,:,3),1));

SRI11 = tmprod(SRI11,W',3);      
rgb4(:,:,1) = imadjust(rescale(SRI11(:,:,1),1));
rgb4(:,:,2) = imadjust(rescale(SRI11(:,:,2),1));
rgb4(:,:,3) = imadjust(rescale(SRI11(:,:,3),1));

SRI31 = tmprod(SRI31,W',3);  
rgb6(:,:,1) = imadjust(rescale(SRI31(:,:,1),1));
rgb6(:,:,2) = imadjust(rescale(SRI31(:,:,2),1));
rgb6(:,:,3) = imadjust(rescale(SRI31(:,:,3),1));  

SRI32 = tmprod(SRI32,W',3);  
rgb5(:,:,1) = imadjust(rescale(SRI32(:,:,1),1));
rgb5(:,:,2) = imadjust(rescale(SRI32(:,:,2),1));
rgb5(:,:,3) = imadjust(rescale(SRI32(:,:,3),1));   
           

figure(1)
subaxis(2,2,1); imshow(rgb); title('Ref.')
set(gca,'FontName','Times','FontSize',16);
subaxis(2,2,2); imshow(rgb4); title('SCOTT')
set(gca,'FontName','Times','FontSize',16);
subaxis(2,2,3); imshow(rgb5); title('Alg. 3 [4 4]')
set(gca,'FontName','Times','FontSize',16);
subaxis(2,2,4); imshow(rgb6); title('Alg. 3 [2 2]')
set(gca,'FontName','Times','FontSize',16);
      

      


%% Tables



table1 = ["Algorithm" "R-SNR" "CC" "SAM" "ERGAS" "Time (sec)";
    "Best" "Infty" "1" "0" "0" "0"
    "SCOTT" err11{:} t11; 
    "Alg. 3 [4 4]" err31{:} t31; 
     "Alg. 3 [2 2]" err32{:} t32;]
      

