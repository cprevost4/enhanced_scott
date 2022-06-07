function SRI_hat = escott(HSI, MSI, P1, P2, Pm, R, opts)

q = 9; l = opts.lambda;

P1 = sparse(P1); P2 = sparse(P2); Pm = sparse(Pm);

if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'Nblocks') || isempty(opts.Nblocks)
    opts.Nblocks = [1,1];
end

range_MSI = [size(MSI,1),size(MSI,2)]; 
range_HSI = [size(HSI,1),size(HSI,2)];
step_MSI = ceil(range_MSI ./ opts.Nblocks); 
step_HSI = ceil(range_HSI ./ opts.Nblocks);
SRI_hat = zeros(size(MSI,1), size(MSI,2), size(HSI,3));

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
    
%      [fact1] = mlsvd(HSI(ind_HSI{1},ind_HSI{2}, :), R); % hosvd
%       Uh = cell2mat(fact1(1)); Vh = cell2mat(fact1(2)); Wh = cell2mat(fact1(3));
% %      
%      [fact2] = mlsvd(MSI(ind_MSI{1},ind_MSI{2}, :), R); % hosvd
%       Um = cell2mat(fact2(1)); Vm = cell2mat(fact2(2)); Wm = cell2mat(fact2(3));  

[Um, ~, ~] = svds(tens2mat(MSI(ind_MSI{1},ind_MSI{2}, :),1,[]),R(1));
[Vm, ~, ~] = svds(tens2mat(MSI(ind_MSI{1},ind_MSI{2}, :),2,[]),R(2));
[Wm, ~, ~] = svds(tens2mat(MSI(ind_MSI{1},ind_MSI{2}, :),3,[]), R(3));
[Uh, ~, ~] = svds(tens2mat(HSI(ind_HSI{1},ind_HSI{2}, :),1,[]),R(1));
[Vh, ~, ~] = svds(tens2mat(HSI(ind_HSI{1},ind_HSI{2}, :),2,[]),R(2));
[Wh, ~, ~] = svds(tens2mat(HSI(ind_HSI{1},ind_HSI{2}, :),3,[]), R(3));
     
%     tempU = (P1(ind_HSI{1},ind_MSI{1}) * Um) \ Uh;
%     U = Um*tempU;
    A = Um*Um'*P1(ind_HSI{1},ind_MSI{1})';
    U = (Uh'*A' + (1/l)*Um')*inv(A*A'+(1/l)*eye(size(A*A',1))); U = U';
    %U = (P1(ind_HSI{1},ind_MSI{1}) * Um * Um') \ Uh;
    Utilde = P1(ind_HSI{1},ind_MSI{1}) * U;

%     tempV = (P2(ind_HSI{2},ind_MSI{2}) * Vm) \ Vh;
%     V = Vm*tempV;
%     %V = (P2(ind_HSI{2},ind_MSI{2}) * Vm * Vm') \ Vh;

    A = Vm*Vm'*P2(ind_HSI{2},ind_MSI{2})';
    V = (Vh'*A' + (1/l)*Vm')*inv(A*A'+(1/l)*eye(size(A*A',1))); V = V';
    Vtilde = P2(ind_HSI{2},ind_MSI{2}) * V;
    
    
    A = Wh*Wh'*Pm';
    W = (Wm'*A' + l*Wh')*inv(A*A'+l*eye(size(A*A',1))); W = W';
    Wtilde = Pm*W;
    
    left = l*kron(W'*W,kron(Vtilde'*Vtilde,Utilde'*Utilde)) + kron(Wtilde'*Wtilde,kron(V'*V,U'*U));
    right = l*lmlragen({Utilde',Vtilde',W'},HSI(ind_HSI{1},ind_HSI{2}, :)) + lmlragen({U',V',Wtilde'},MSI(ind_MSI{1},ind_MSI{2}, :));
    S = left\right(:);
    %[S, flag] = pcg(left,right(:));
    S = reshape(S,R); 
    



    SRI_hat(ind_MSI{1},ind_MSI{2}, :) = lmlragen({U,V,W},S); 
    
    
    
  end
end

end

