function [clusterResult] = myKmeans(X, cluster_amount, threshold)  
% ===============================================================================================  
% K-means clustering method
% Inputs:
%   X              : The input matric, each line is an instance and each column is an attribute
%   cluster_amount : The amount of clusters
%   threshold      : The condition for stop interation
%  
% Outputs:
%   clusterResult  : A n * 1 dimensional matrix
%                    Each line in the clusterResult indicate the clustering result
% =============================================================================================== 

    %% Generate Initial Conditions
    [N, D] = size(X);
    k = cluster_amount;
    % randomly choose k points as the initial centroid
    % centroid_point_number = randperm(N, k); 
    % fixed initial condition
    centroid_point_number = [10,998,5004,7021];
    % centroid_point_number = [4545,6444,955,4696,6820,7195,6979,6999];
    % centroid_point_number = [1450,5600,10,950,5100,4800,6700,7195,6900,6990];
    centroids = X(centroid_point_number, :);
    
    while true
        new_centroids = zeros(k,D);
        dist = zeros(1,k);
        labels = zeros(N,1);
        num = zeros(1,k);
        % generate cluster results
        for i = 1:N
            for j = 1:k
                dist(j) = norm(X(i,:) - centroids(j,:));
            end
            [~, temp] = min(dist);
            labels(i) = temp;
        end
        % calculate the new centroid position
        convergence_count = 0;
        for i = 1:k
            for j = 1:N
                if labels(j) == i
                    new_centroids(i,:) = new_centroids(i,:) + X(j,:);
                    num(i) = num(i) + 1;
                end
            end
            new_centroids(i,:) = new_centroids(i,:) ./ num(1,i);
            % check whether in convergence condition
            if norm(new_centroids(i,:) - centroids(i,:)) < threshold
                convergence_count = convergence_count + 1;
            end
        end
        
        if convergence_count == k
            break
        else
            centroids = new_centroids;
        end
    % end for while loop    
    end
    %% return results
    clusterResult = labels;
 end  