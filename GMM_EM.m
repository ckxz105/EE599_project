function [clusterResult] = GMM_EM(X, cluster_amount, threshold)  
% ===============================================================================================  
% GMM-EM model clustering method
% Inputs:
%   X              : The input matric, each line is an instance and each column is an attribute
%   cluster_amount : The amount of clusters
%   threshold      : The condition for stop interation
%  
% Outputs:
%   clusterResult  : A n * cluster_amount dimensional matrix
%                    Each element in the clusterResult indicate the possibility of the point belong to the cluster
% =============================================================================================== 


    %% Generate Initial Conditions
    [N, D] = size(X);
    k = cluster_amount;
    % randomly choose k points as the initial centroid
    %centroid_point_number = randperm(N, k); 
    % fixed initial condition
    centroid_point_number = [100,1111,5555,7021];
    % centroid_point_number = [4545,6444,955,4696,6820,7195,6979,6999];
    % centroid_point_number = [1450,5600,10,950,5100,4800,6700,7195,6900,6990];
    centroids = X(centroid_point_number, :);
    val_miu = centroids;  
    val_phi = zeros(1, k);
    val_sigma = zeros(D, D, k);
    % calculate the distance from every point to the centroid and choose
    % the smallest distance one as the clsuter result
    dist = zeros(1,k);
    labels = zeros(N,1);
    for i = 1:N
        for j = 1:k
            dist(j) = norm(X(i,:) - centroids(j,:));
        end
        [~, temp] = min(dist);
        labels(i) = temp;
    end

    % update the parameters of phi and sigma
    for i = 1:k  
        Xi = X(labels == i, :); 
        val_phi(i) = size(Xi, 1)/N; 
        val_sigma(:, :, i) = cov(Xi); % use covariance matrix to find sigma
    end
    Lprev = -inf;  
    
    %% E/M steps
    while true % repeat until convergence  
        % step E, calculate the posteriori probability
        Px = zeros(N, k);
        for i = 1:k  
            Xshift = X - repmat(val_miu(i, :), N, 1); 
            inv_pSigma = inv(val_sigma(:, :, i));
            tmp = sum((Xshift*inv_pSigma) .* Xshift, 2);   
            coef = (2*pi)^(-D/2) * sqrt(det(inv_pSigma));  
            Px(:, i) = coef * exp(-0.5*tmp);
        end 
           
        % step M, update new parameters for miu, phi and sigma
         % new value for pGamma  
        val_gamma = Px .* repmat(val_phi, N, 1); 
        val_gamma = val_gamma ./ repmat(sum(val_gamma, 2), 1, k); 
        Ni = sum(val_gamma, 1); 
        val_miu = diag(1./Ni) * val_gamma' * X;  
        val_phi = Ni/N; 
        for ii = 1:k 
            Xshift = X-repmat(val_miu(ii, :), N, 1);
            val_sigma(:, :, ii) = (Xshift' * (diag(val_gamma(:, ii)) * Xshift)) / Ni(ii);  
        end  

        % check whether convergence or not base on likelihood function
        L = sum(log(Px*val_phi'));
        % disp(L-Lprev);
        
        % if convergence, break the loop
        disp(abs(L-Lprev));
        if abs(L-Lprev) < threshold
             break; 
        end
        Lprev = L; % update L
    end 
    
    %% return result
     clusterResult = Px;  

 end  