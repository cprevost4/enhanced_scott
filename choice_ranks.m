clear all
close all
clc
%rng(107)

%% Lockwood

load('Lockwood_preproc_subim1')
% clear MSI
SRI = HSI; clear HSI; SRI = SRI(1:88,1:88,:);
P3 = SRF; clear SRF; 

d1 = 2; d2 = 2; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
MSI = tmprod(SRI,P3,3);

HSI = awgn(HSI,40,'measured');
MSI = awgn(MSI,10,'measured');
sigma_h = 10^(-40/10); sigma_m = 10^(-10/10);
opts.lambda = (sigma_h^2)./(sigma_m^2);

opts.Nblocks = [2,2];
Rmax = floor([size(HSI,1)/opts.Nblocks(1) size(HSI,2)/opts.Nblocks(2) size(MSI,3)]);

for R1 = 1:Rmax(1)
    for R3 = 1:Rmax(3)
        
        R = [R1 R1 R3]
        SRI_hat = escott(HSI, MSI, P1, P2, P3, R, opts);
        error(R1,R3) = r_snr(SRI,SRI_hat);
        
    end
end

figure
imagesc(error); colorbar; xlabel('R_1'); ylabel('R_3'); title('R-SNR as a function of ranks')

%% Isabella Lake

load('Isabella_lake_preproc_subim1')
 clear MSI
SRI = HSI; clear HSI; %SRI = SRI(1:96,1:96,:);
P3 = SRF; clear SRF; 
d1 = 2; d2 = 2; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
MSI = tmprod(SRI,P3,3);

HSI = awgn(HSI,30,'measured');
MSI = awgn(MSI,30,'measured');
opts.lambda = 1;

opts.Nblocks = [2 2];
Rmax = floor([size(HSI,1)/opts.Nblocks(1) size(HSI,2)/opts.Nblocks(2) size(MSI,3)]);

for R1 = 1:Rmax(1)
    for R3 = 1:Rmax(3)
        
        R = [R1 R1 R3]
        SRI_hat = escott(HSI, MSI, P1, P2, P3, R, opts);
        error(R1,R3) = r_snr(SRI,SRI_hat);
        
    end
end


figure
imagesc(error); colorbar; xlabel('R_1'); ylabel('R_3'); title('R-SNR as a function of ranks')

%%

amoy = 0;

opts.Nblocks = [4,4];

range_MSI = [size(MSI,1),size(MSI,2)]; 
range_HSI = [size(HSI,1),size(HSI,2)];
step_MSI = ceil(range_MSI ./ opts.Nblocks); 
step_HSI = ceil(range_HSI ./ opts.Nblocks);

for i1=1:opts.Nblocks(1)
  for i2=1:opts.Nblocks(2)
      
    %Update steps
    M_ind_min = [i1-1,i2-1].*step_MSI + 1;
    M_ind_max = min([i1,i2].*step_MSI, range_MSI);
    H_ind_min = [i1-1,i2-1].*step_HSI + 1;
    H_ind_max = min([i1,i2].*step_HSI, range_HSI);
    %Range depending on blocks
    if i1==1
        ind_MSI{1} = M_ind_min(1):M_ind_max(1)+(q-1)/2;
        ind_HSI{1} = H_ind_min(1):H_ind_max(1)+(q-1)/2;
    elseif i1==opts.Nblocks(1)    
        ind_MSI{1} = M_ind_min(1)-(q-1)/2:M_ind_max(1);
        ind_HSI{1} = H_ind_min(1)-(q-1)/2:H_ind_max(1);
    else
        ind_MSI{1} = M_ind_min(1)-(q-1)/2:M_ind_max(1)+(q-1)/2;
        ind_HSI{1} = H_ind_min(1)-(q-1)/2:H_ind_max(1)+(q-1)/2;
    end
    if i2==1
        ind_MSI{2} = M_ind_min(2):M_ind_max(2)+(q-1)/2;
        ind_HSI{2} = H_ind_min(2):H_ind_max(2)+(q-1)/2;
    elseif i2==opts.Nblocks(2)
        ind_MSI{2} = M_ind_min(2)-(q-1)/2:M_ind_max(2);
        ind_HSI{2} = H_ind_min(2)-(q-1)/2:H_ind_max(2);
    else
        ind_MSI{2} = M_ind_min(2)-(q-1)/2:M_ind_max(2)+(q-1)/2;
        ind_HSI{2} = H_ind_min(2)-(q-1)/2:H_ind_max(2)+(q-1)/2;
    end
    
    
    mat = HSI(ind_HSI{1},ind_HSI{2}, :);
    amoy = amoy + svd(tens2mat(mat,[],3));
     
  end
end

test = amoy/(i1*i2);

    figure(1)
     subplot(1,2,1); semilogy(test,'r.','MarkerSize',10,'LineWidth',1); hold on
     xlabel('Index of spectral bin','interpreter','latex'); title('$[4,4]$ pattern','interpreter','latex');
     xlim([1 30])
set(gca,'FontName','Times','FontSize',16); 