% 这是一个测试所有目标函数性能的算法
%'DME -1',, 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'
%随机打乱治病snp之后的代码
models = {  'DME -1','DME -2', 'DME -3', 'DME -4', 'DME -5', 'DME -6', 'DME -7', 'DME -8', 'DME -9', 'DME -10', 'DME -11', 'DME -12','DNME -1', 'DNME -2', 'DNME -3', 'DNME -4', 'DNME -5', 'DNME -6', 'DNME -7', 'DNME -8', 'DNME -9', 'DNME -10', 'DNME -11', 'DNME -12', 'DNME -13', 'DNME -14'};

epi_dim=2;
s = 2;
sample=1600;
step = 0.1;



for ln=1:step:1
    for i = 1:numel(models) %每一个致病模型
        model = models{i};
        %workspace = sprintf('C:\\Users\\ASUS\\Desktop\\实验\\模拟数据集\\模拟数据集\\有边际效应\\%s', model)
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
                    snp=[];
                    hanshu=[];
                    State=data(:,n);
                    for k=1:99
                        for l=k+1:100
                            combined_data = [data(:, k), data(:, l)];
                            [JE_Score,Gtest,JS,ME,K2Score]=multi_criteriaEvaluationFuns2021(combined_data,State);
                            %[LR] = MDR_2020_2(combined_data,State);
                            % 判断 hanshu 是否为空
                            if isempty(hanshu)
                                % 如果为空，则将 K2Score 加入到 hanshu 中
                               % hanshu = K2Score;
                               hanshu =JE_Score;
                                % 将 k 和 l 合并为一个数组，放入 snp 中
                                snp = [k, l];
                            else
                                % 如果 hanshu 不为空，比较 K2Score 和 hanshu 中的值
                                if JE_Score >= hanshu
                                    % 如果 K2Score 更小，替换 hanshu 中的值
                                    hanshu =JE_Score;
                                    % 替换 snp 的值为当前的 k 和 l
                                    snp = [k, l];
                                end
                            end
                            %if k==99&&l==100
                            %    K2Score
                            %end

                        end
                    end

                    if snp==[99,100]
                        num=num+1;
                    end
                    num




                end
            end
            fos = fopen(['C:\\Users\\tianjin\\Desktop\\data\\论文\\代码\\sun\\HS-MMGKG-master\\data\\simulated_data\\目标函数\\JE检测结果\\',model,'.result.txt'], 'w');
             fprintf(fos, '%s', num2str(num));
            fclose(fos);
        catch exception
            disp(exception);
        end
    end
end




