%随机打乱治病snp之后的代码
models = {'DME -11', 'DME -12', 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'};

epi_dim=2;
s = 2;
sample=1600;
step = 0.1;



for ln=1:step:1
    for i = 1:numel(models) %每一个致病模型
        model = models{i};
        %workspace = sprintf('C:\\Users\\ASUS\\Desktop\\实验\\模拟数据集\\模拟数据集\\有边际效应\\%s', model)
        % workspace = sprintf('C:\\Users\\tianjin\\Desktop\\experiment\\模拟数据集\\无边际效应\\%s', model)
        workspace = sprintf('C:\\Users\\tianjin\\Desktop\\data\\论文\\代码\\sun\\HS-MMGKG-master\\data\\simulated_data\\%s', model)
        try
            num=0;
            filenameList = dir(workspace);
            for j = 1:numel(filenameList) %每个模型中的某个数据集
                if ~filenameList(j).isdir
                    filename = filenameList(j).name;
                    fileFullPath = fullfile(workspace, filename);
                    data=dlmread(fileFullPath,'\t',1,0);
                    [m,n] = size(data);
                    Dim = n - 1;
                    State=data(:,n);
                    % HMS =2*max(100, epi_dim*min(Dim/10,100));
                    HMS=20;
                    CX = Dim-epi_dim+1:Dim;
                    TP0 = 0.5;%全0
                    PAR=0.5;
                    F = 20;
                    HMCR = 0.98;
                    %max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000) ;
                    max_iter=3600;
                    % max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000);
                    % 选择要交换的两列
                    col_to_swap_1 = 99;
                    col_to_swap_2 = 100;
                    
                    % 从第 1000 列之前的列中随机选择两列 A 和 B
                    
                    col_A = randi([1, col_to_swap_2]);
                    col_B = randi([1, col_to_swap_2]);
                    
                    
                    % 交换第 999 列和第 1000 列与列 A 和 B 的数据
                    temp = data(:, col_to_swap_1);
                    data(:, col_to_swap_1) = data(:, col_A);
                    data(:, col_A) = temp;
                    temp = data(:, col_to_swap_2);
                    data(:, col_to_swap_2) = data(:, col_B);
                    data(:, col_B) = temp;
                    
                    [Task,NC,flag,Epi_Dim_FEs,ci] = HS_2021_multiTask_UnifiedCoding2021(data, epi_dim, s, HMS,max_iter,CX, TP0, PAR,F,HMCR,sample,ln);
                    
                    pvalue = 1/nchoosek(Dim,epi_dim);
                    %获取task中的snp组合，进行gtest评分，低于pvalue的值作为候选解，将这些组合放入文件中，争取二阶三阶，4阶都找五个
                    m =3;
                    n =4;
                    p =4;
                    result = {};
                    % 要查找的元素
                    element1 = col_A ;
                    element2 = col_B;
                    
                    
                    for i = 1:m
                        for j = 1:n
                            for k = 1:p
                                Xnew=Task(i).Elite(j).Xbest(k,1:i+1);
                                tem=Gtest_score(data(:,Xnew(1,:)),State);
                                if ismember(element1, Task(i).Elite(j).Xbest(k,1:2)) && ismember(element2, Task(i).Elite(j).Xbest(k,1:2))
                                    if tem< pvalue
                                        result{end+1} = struct('Xnew', Xnew, 'tem', tem);
                                    end
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
                    
                    unique_results
                    
                    
                end
            end
            
        catch exception
            disp(exception);
        end
    end
end




