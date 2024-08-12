%��������β�snp֮��Ĵ��� 
models = {'DME -11', 'DME -12', 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'};
%models = {'1_EDM-1','2_EDM-1','3_EDM-1','4_EDM-1','5_EDM-1','6_EDM-1','7_EDM-1','8_EDM-1','9_EDM-1','10_EDM-1','11_EDM-1','12_EDM-1','13_EDM-1','14_EDM-1','15_EDM-1','16_EDM-1','17_EDM-1','18_EDM-1','19_EDM-1','20_EDM-1','21_EDM-1','22_EDM-1','23_EDM-1','24_EDM-1','25_EDM-1','26_EDM-1','27_EDM-1','28_EDM-1','29_EDM-1','30_EDM-1','31_EDM-1','32_EDM-1','33_EDM-1','34_EDM-1','35_EDM-1','36_EDM-1','37_EDM-1','38_EDM-1','39_EDM-1','40_EDM-1'};
%models = {'1_EDM-1','2_EDM-1','3_EDM-1','4_EDM-1','5_EDM-1','6_EDM-1','7_EDM-1','8_EDM-1','9_EDM-1','10_EDM-1','11_EDM-1','12_EDM-1'};
%models = {'DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DME -1', 'DME -2', 'DME -3', 'DME -4'};

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
                   % HMS =2*max(100, epi_dim*min(Dim/10,100));
                    HMS=40;
                    CX = Dim-epi_dim+1:Dim;
                    TP0 = 0.5;%ȫ0
                    PAR=0.5;
                    F = 20;
                    HMCR = 0.98;
                    %max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000) ;
                    max_iter=1600;
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

                    % ���� Task �����ά��Ϊ (m, n, p)������Ҫ���ú��ʵ� m��n �� p ��ֵ
                    m =1;
                    n =4;
                    p =4;

                    % Ҫ���ҵ�Ԫ��
                    element1 = col_A ;
                    element2 = col_B;
                    foundFlag = false;
                    % ���� Task ����
                    for i = 1:m
                        for j = 1:n
                            for k = 1:p
                                % ����Ƿ�ͬʱ���� element1 �� element2
                                if ismember(element1, Task(i).Elite(j).Xbest(k,1:2)) && ismember(element2, Task(i).Elite(j).Xbest(k,1:2))
                                    num=num+1;
                                    foundFlag = true;  % ���ñ�־����Ϊ true
                                    break;
                                end
                            end
                            if foundFlag
                                break;  % �����־Ϊ true���˳��ڶ���ѭ��
                            end
                        end
                        if foundFlag
                            break;  % �����־Ϊ true���˳�������ѭ��
                        end
                    end
                    num
                end
            end
            % �������� model��pxnew �� psearch �������ļ���
            % filename = sprintf('C:\\Users\\ASUS\\Desktop\\ʵ��\\ģ�����ݼ�\\ģ�����ݼ�\\�б߼�ЧӦ\\�����ʼ��10000�ε�������\\%s_%f_%f.result.txt', model, pxnew, psearch);
            % filename = sprintf('C:\\Users\\tianjin\\Desktop\\experiment\\ģ�����ݼ�\\�ޱ߼�ЧӦ\\ours10wͳ��2\\%s.result.txt', model);
             filename = sprintf('C:\\Users\\tianjin\\Desktop\\data\\����\\����\\sun\\HS-MMGKG-master\\data\\simulated_data\\ours����\\%s.result.txt', model);
             fos = fopen(filename, 'w');
            fprintf(fos, '%s', num2str(num));
            fclose(fos);
        catch exception
            disp(exception);
        end
    end
end




