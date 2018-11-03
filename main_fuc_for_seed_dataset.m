%% EE599 Project code
% for Seed Data Set
clear;
%% Read files
path = fullfile('/Users/minazuki/Desktop/studies/master/2018Fall/EE599/EE599_project/','seeds_dataset.txt');
fidin = fopen(path,'r');
line = 0;
data = zeros(0,7);
% column_name = fgetl(fidin);
result = zeros(0,1);
while ~feof(fidin)
    tline=fgetl(fidin);
    temp_read = regexp(tline,'	','split');
    line = line + 1;
    % column 1 to 22 are the attributes
    data(line,:) = str2double(temp_read(1:7));
    % column 23 to 25 are the cluster results
     result(line,:) = str2double(temp_read(8));
end
fclose(fidin);
% clear unused variables
clear fidin path temp_read tline;
%% use GMM-EM algorithm to get cluster results
tic;
k = 3;
threshold_GMM = 1e-15;
GMM_EM_res = GMM_EM(data,k,threshold_GMM);
% in res, each column of the line with maximum value is the cluster results
[~, labels] = max(GMM_EM_res, [], 2); 
labels_adjusted1 = (labels == labels(1))* 1 + (labels == labels(119)) * 2 + (labels == labels(end)) * 3;
t1 = toc;
disp('results of GMM_EM');
Evaluate(result, labels_adjusted1)
disp('time =');
disp(t1);
%% use K-means to get cluster results
tic;
threshold_k_means = 0.1;
labels_k_means = myKmeans(data, k, threshold_k_means);
labels_adjusted2 = (labels_k_means == labels_k_means(1))* 1 + (labels_k_means == labels_k_means(119)) * 2 + (labels_k_means == labels_k_means(end)) * 3;
t2 = toc;
disp('results of K-means');
Evaluate(result, labels_adjusted2)
disp('time =');
disp(t2);