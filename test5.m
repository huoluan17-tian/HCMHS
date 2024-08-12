%随机打乱治病snp之后的代码 
models = {'DME -11', 'DME -12', 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'};
%models = {'1_EDM-1','2_EDM-1','3_EDM-1','4_EDM-1','5_EDM-1','6_EDM-1','7_EDM-1','8_EDM-1','9_EDM-1','10_EDM-1','11_EDM-1','12_EDM-1','13_EDM-1','14_EDM-1','15_EDM-1','16_EDM-1','17_EDM-1','18_EDM-1','19_EDM-1','20_EDM-1','21_EDM-1','22_EDM-1','23_EDM-1','24_EDM-1','25_EDM-1','26_EDM-1','27_EDM-1','28_EDM-1','29_EDM-1','30_EDM-1','31_EDM-1','32_EDM-1','33_EDM-1','34_EDM-1','35_EDM-1','36_EDM-1','37_EDM-1','38_EDM-1','39_EDM-1','40_EDM-1'};
%models = {'1_EDM-1','2_EDM-1','3_EDM-1','4_EDM-1','5_EDM-1','6_EDM-1','7_EDM-1','8_EDM-1','9_EDM-1','10_EDM-1','11_EDM-1','12_EDM-1'};
%models = {'DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DME -1', 'DME -2', 'DME -3', 'DME -4'};

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
                   % HMS =2*max(100, epi_dim*min(Dim/10,100));
                    HMS=40;
                    CX = Dim-epi_dim+1:Dim;
                    TP0 = 0.5;%全0
                    PAR=0.5;
                    F = 20;
                    HMCR = 0.98;
                    %max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000) ;
                    max_iter=1600;
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

                    % 假设 Task 数组的维度为 (m, n, p)，你需要设置合适的 m、n 和 p 的值
                    m =1;
                    n =4;
                    p =4;

                    % 要查找的元素
                    element1 = col_A ;
                    element2 = col_B;
                    foundFlag = false;
                    % 遍历 Task 数组
                    for i = 1:m
                        for j = 1:n
                            for k = 1:p
                                % 检查是否同时包含 element1 和 element2
                                if ismember(element1, Task(i).Elite(j).Xbest(k,1:2)) && ismember(element2, Task(i).Elite(j).Xbest(k,1:2))
                                    num=num+1;
                                    foundFlag = true;  % 设置标志变量为 true
                                    break;
                                end
                            end
                            if foundFlag
                                break;  % 如果标志为 true，退出第二层循环
                            end
                        end
                        if foundFlag
                            break;  % 如果标志为 true，退出第三层循环
                        end
                    end
                    num
                end
            end
            % 构建包含 model、pxnew 和 psearch 参数的文件名
            % filename = sprintf('C:\\Users\\ASUS\\Desktop\\实验\\模拟数据集\\模拟数据集\\有边际效应\\聚类初始化10000次迭代补充\\%s_%f_%f.result.txt', model, pxnew, psearch);
            % filename = sprintf('C:\\Users\\tianjin\\Desktop\\experiment\\模拟数据集\\无边际效应\\ours10w统计2\\%s.result.txt', model);
             filename = sprintf('C:\\Users\\tianjin\\Desktop\\data\\论文\\代码\\sun\\HS-MMGKG-master\\data\\simulated_data\\ours测试\\%s.result.txt', model);
             fos = fopen(filename, 'w');
            fprintf(fos, '%s', num2str(num));
            fclose(fos);
        catch exception
            disp(exception);
        end
    end
end




