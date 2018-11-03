%% EE599 Project code
% for Seed Data Set
clear;
%% Read files
path = fullfile('/Users/minazuki/Desktop/studies/master/2018Fall/EE599/EE599_project/Anuran Calls (MFCCs)','Frogs_MFCCs.csv');
fidin = fopen(path,'r');
line = 0;
data = zeros(0,22);
column_name = fgetl(fidin);
result = strings(0,3);
while ~feof(fidin)
    tline=fgetl(fidin);
    temp_read = regexp(tline,',','split');
    line = line + 1;
    % column 1 to 22 are the attributes
    data(line,:) = str2double(temp_read(1:22));
    % column 23 to 25 are the cluster results
    result(line,:) = temp_read(23:25);
end
fclose(fidin);
% clear unused variables
clear fidin path temp_read tline;
%% use GMM-EM algorithm to get cluster results
tic;
% k = 10;
% k = 8;
k = 4;
[COEFF,SCORE,latent,tsquare] = pca(data);
data_pca = data * COEFF(:,1:13);
threshold_GMM = 1e-10;
GMM_EM_res = GMM_EM(data_pca,k,threshold_GMM);
ground_truth_k_4 = result(:,1);
ground_truth_k_4 = (ground_truth_k_4 == 'Leptodactylidae')* 1 + (ground_truth_k_4 == 'Hylidae') * 2 + (ground_truth_k_4 == 'Dendrobatidae') * 3 +  + (ground_truth_k_4 == 'Bufonidae') * 4;
ground_truth_k_8 = result(:,2);
ground_truth_k_8 = ...
     (ground_truth_k_8 == 'Adenomera')* 1 + (ground_truth_k_8 == 'Hypsiboas')* 2 ...
    +(ground_truth_k_8 == 'Ameerega')* 3 +(ground_truth_k_8 == 'Dendropsophus')* 4 ...
    +(ground_truth_k_8 == 'Leptodactylus')* 5 +(ground_truth_k_8 == 'Scinax')* 6 ...
    +(ground_truth_k_8 == 'Osteocephalus')* 7 +(ground_truth_k_8 == 'Rhinella')* 8;
ground_truth_k_10 = result(:,3);
ground_truth_k_10 = ...
     (ground_truth_k_10 == 'AdenomeraHylaedactylus')* 1 + (ground_truth_k_10 == 'HypsiboasCordobae')* 2 ...
    +(ground_truth_k_10 == 'AdenomeraAndre')* 3 +(ground_truth_k_10 == 'Ameeregatrivittata')* 4 ...
    +(ground_truth_k_10 == 'HypsiboasCinerascens')* 5 +(ground_truth_k_10 == 'HylaMinuta')* 6 ...
    +(ground_truth_k_10 == 'LeptodactylusFuscus')* 7 +(ground_truth_k_10 == 'ScinaxRuber')* 8 ...
    +(ground_truth_k_10 == 'OsteocephalusOophagus')* 9 +(ground_truth_k_10 == 'Rhinellagranulosa')* 10;
% in res, each column of the line with maximum value is the cluster results
[~, labels] = max(GMM_EM_res, [], 2);
labels_adjusted1 = zeros(7195,1);
correspond = [4,3,1,2];
for i = 1: 7195
    labels_adjusted1(i) = correspond(labels(i));
end
t1 = toc;
disp('results of GMM_EM');
% change here to get score of different values
Evaluate(ground_truth_k_4, labels_adjusted1)
disp('time =');
disp(t1);
%% use K-means to get cluster results
tic;
threshold_k_means = 0.1;
labels_k_means = myKmeans(data, k, threshold_k_means);
tabulate(labels_k_means(:))
correspond = [3,4,2,1];
labels_adjusted2 = zeros(7195,1);
for i = 1: 7195
    labels_adjusted2(i) = correspond(labels_k_means(i));
end
t2 = toc;
disp('results of K-means');
% change here to get score of different values
Evaluate(ground_truth_k_4, labels_adjusted2)
disp('time =');
disp(t2);