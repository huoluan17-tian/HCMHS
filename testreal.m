%300*195702
fileFullPath = 'C:\Users\tianjin\Desktop\experiment\real\印度冠状\数据集\转置_真实.txt';
dataTable = readtable(fileFullPath, 'Delimiter', '\t', 'ReadVariableNames', false);
data = table2array(dataTable);
for ite=1:100
    ite
epi_dim=3;
s = 2;
[m,n] = size(data);
Dim = n - 1;
HMS =2*max(100, epi_dim*min(Dim/10,100));
max_iter=1000000;
CX=Dim-epi_dim+1:Dim;
TP0 = 0.5;%全0
PAR=0.5;
F = 20;
HMCR = 0.98;
sample=300;
ln=1;
pvalue = 1/nchoosek(Dim,epi_dim);
[Task,NC,flag,Epi_Dim_FEs,ci] = HS_2021_multiTask_UnifiedCoding2021(data, epi_dim, s, HMS,max_iter,CX, TP0, PAR,F,HMCR,sample,ln);
Task;
%在这里加入G-test 检验
State=data(:,n);


%获取task中的snp组合，进行gtest评分，低于pvalue的值作为候选解，将这些组合放入文件中，争取二阶三阶，4阶都找五个
m =2;
n =4;
p =40;
result = {};


for i = 2:m
    for j = 1:n
        for k = 1:p
            Xnew=Task(i).Elite(j).X(k,1:i+1);
            tem=Gtest_score(data(:,Xnew(1,:)),State);
            if tem< pvalue&&tem>0
                result{end+1} = struct('Xnew', Xnew, 'tem', tem);
            end
            
            
        end
        
    end
    
end

unique_results = {};

% 遍历result中的每个结构体
for idx = 1:numel(result)
    % 获取当前结构体
    current_result = result{idx};
    
    % 检查当前结构体是否已经存在于unique_results中
    is_unique = true;
    for j = 1:numel(unique_results)
        % 如果Xnew和tem都相同，则当前结果不是唯一的
        if isequal(current_result.Xnew, unique_results{j}.Xnew) && current_result.tem == unique_results{j}.tem
            is_unique = false;
            break;
        end
    end
    
    % 如果当前结果是唯一的，则将其添加到unique_results中
    if is_unique
        unique_results{end+1} = current_result;
    end
end
% 将tem设置为longG格式
for idx = 1:numel(unique_results)
    format longG; % 设置输出格式为longG
    unique_results{idx}.tem = unique_results{idx}.tem;
end

tems = cellfun(@(x) x.tem, unique_results);

% 排序
[~, sorted_indices] = sort(tems);

% 根据排序后的索引重新排列unique_results
unique_results = unique_results(sorted_indices);
% 打印Xnew和tem
for idx = 1:numel(unique_results)
    disp(['Xnew: ', num2str(unique_results{idx}.Xnew)]);
    disp(['tem: ', num2str(unique_results{idx}.tem)]);
end
end
