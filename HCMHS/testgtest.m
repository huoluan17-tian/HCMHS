%��������β�snp֮��Ĵ���
models = {'DME -11', 'DME -12', 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'};

epi_dim=2;
s = 2;
sample=1600;
step = 0.1;



for ln=1:step:1
    for i = 1:numel(models) %ÿһ���²�ģ��
        model = models{i};
        %workspace = sprintf('C:\\Users\\ASUS\\Desktop\\ʵ��\\ģ�����ݼ�\\ģ�����ݼ�\\�б߼�ЧӦ\\%s', model)
        % workspace = sprintf('C:\\Users\\tianjin\\Desktop\\experiment\\ģ�����ݼ�\\�ޱ߼�ЧӦ\\%s', model)
        workspace = sprintf('C:\\Users\\tianjin\\Desktop\\data\\����\\����\\sun\\HS-MMGKG-master\\data\\simulated_data\\%s', model)
        try
            num=0;
            filenameList = dir(workspace);
            for j = 1:numel(filenameList) %ÿ��ģ���е�ĳ�����ݼ�
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
                    TP0 = 0.5;%ȫ0
                    PAR=0.5;
                    F = 20;
                    HMCR = 0.98;
                    %max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000) ;
                    max_iter=3600;
                    % max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000);
                    % ѡ��Ҫ����������
                    col_to_swap_1 = 99;
                    col_to_swap_2 = 100;
                    
                    % �ӵ� 1000 ��֮ǰ���������ѡ������ A �� B
                    
                    col_A = randi([1, col_to_swap_2]);
                    col_B = randi([1, col_to_swap_2]);
                    
                    
                    % ������ 999 �к͵� 1000 ������ A �� B ������
                    temp = data(:, col_to_swap_1);
                    data(:, col_to_swap_1) = data(:, col_A);
                    data(:, col_A) = temp;
                    temp = data(:, col_to_swap_2);
                    data(:, col_to_swap_2) = data(:, col_B);
                    data(:, col_B) = temp;
                    
                    [Task,NC,flag,Epi_Dim_FEs,ci] = HS_2021_multiTask_UnifiedCoding2021(data, epi_dim, s, HMS,max_iter,CX, TP0, PAR,F,HMCR,sample,ln);
                    
                    pvalue = 1/nchoosek(Dim,epi_dim);
                    %��ȡtask�е�snp��ϣ�����gtest���֣�����pvalue��ֵ��Ϊ��ѡ�⣬����Щ��Ϸ����ļ��У���ȡ�������ף�4�׶������
                    m =3;
                    n =4;
                    p =4;
                    result = {};
                    % Ҫ���ҵ�Ԫ��
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
                    
                    % ����result�е�ÿ���ṹ��
                    for idx = 1:numel(result)
                        % ��ȡ��ǰ�ṹ��
                        current_result = result{idx};
                        
                        % ��鵱ǰ�ṹ���Ƿ��Ѿ�������unique_results��
                        is_unique = true;
                        for j = 1:numel(unique_results)
                            % ���Xnew��tem����ͬ����ǰ�������Ψһ��
                            if isequal(current_result.Xnew, unique_results{j}.Xnew) && current_result.tem == unique_results{j}.tem
                                is_unique = false;
                                break;
                            end
                        end
                        
                        % �����ǰ�����Ψһ�ģ�������ӵ�unique_results��
                        if is_unique
                            unique_results{end+1} = current_result;
                        end
                    end
                    % ��tem����ΪlongG��ʽ
                    for idx = 1:numel(unique_results)
                        format longG; % ���������ʽΪlongG
                        unique_results{idx}.tem = unique_results{idx}.tem;
                    end
                    
                    tems = cellfun(@(x) x.tem, unique_results);
                    
                    % ����
                    [~, sorted_indices] = sort(tems);
                    
                    % ����������������������unique_results
                    unique_results = unique_results(sorted_indices);
                    % ��ӡXnew��tem
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




