function [clusterIndices] = cengci(data,k,sample)
% 随机抽取2000列数据
numCols = size(data, 2);
if numCols>2000
    selectedColumns = datasample(1:numCols, 2000, 'Replace', false); % 从1到numCols中不重复地随机抽取2000个索引
    
    % 获取抽取的列数据
    X = data(:, selectedColumns);
    rSquaredMatrix = corr(X, 'type', 'Pearson').^2; % 计算相关系数并转换为 r^2
    % 获取数据维度
    
    
    % 层次聚类
    distanceMatrix = 1 - rSquaredMatrix; % 将 r^2 转换为距离
    Z = linkage(distanceMatrix, 'average'); % 使用平均连接法计算聚类链接
    
    
    numClusters = k;
    clusters = cluster(Z, 'maxclust', numClusters); % 限制聚类数量
    
    % 获取每个 SNP 的聚类分配
    snpAssignments = clusters;
    snpAssignments(end)=[];
    for i = 1:numClusters
       % clusterIndices{i} = find(snpAssignments == i);
        clusterIndices{i} = selectedColumns(find(snpAssignments == i))';
    end
    
else
    X = data;
    rSquaredMatrix = corr(X, 'type', 'Pearson').^2; % 计算相关系数并转换为 r^2
    % 获取数据维度
    
    
    % 层次聚类
    distanceMatrix = 1 - rSquaredMatrix; % 将 r^2 转换为距离
    Z = linkage(distanceMatrix, 'average'); % 使用平均连接法计算聚类链接
    
    
    numClusters = k;
    clusters = cluster(Z, 'maxclust', numClusters); % 限制聚类数量
    
    % 获取每个 SNP 的聚类分配
    snpAssignments = clusters;
    snpAssignments(end)=[];
    for i = 1:numClusters
        clusterIndices{i} = find(snpAssignments == i);
    end
    
end


end

