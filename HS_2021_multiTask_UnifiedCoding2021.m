%数据集 检测阶数  增加阶数的参数  种群大小 最大迭代  TP, PAR,F,HMCR 算法的控制参数 改造成聚类初始化
function [Task,NC,flag,Epi_Dim_FEs,ci] = HS_2021_multiTask_UnifiedCoding2021(data,epi_dim,s,HMS,max_iter,CX, TP, PAR,F,HMCR,sample,ln)
% 多任务统一编码
%input--------------------------------------------------------------------
% data-----------------input dataset
% epi_dim--------------the epistasis dimension
% HMS--------------the size of harmony memory(HM)

%%-------------------------------------------------------------------------
% initial arguments
ci=0;
Epi_Dim = epi_dim;
Epi_Dim_FEs = HMS;
%F = 5;
epi_dim = epi_dim + s;
% HMCR=0.98;
% PAR=0.7;
% TP = 0.35;
fdim = length(CX);
n=size(data,2);
State=data(:,n);
K = epi_dim-1 ; % 任务数 2-->epi_dim
FitNum = 4;
bestNum = epi_dim;
%% ---------------------------------------------------------------
EliteSize =min(10*epi_dim,ceil(HMS/5)); %精英子集的大小
% CX
NC=0;
SNPs=n-1;  %% 总SNP个数
flag = -1;
maxFit = [];



clusterIndices=cengci(data,epi_dim,sample);

for k = 1:K %循环每一个任务 这里体现多任务  初始化获得task 和elite 5阶任务计算到了7阶
    %% 初始化

    dim = k+1; %计算某一个阶数
    Task(k).X=zeros(HMS,epi_dim);%初始化一个空的 Task(k).X 矩阵，用于存放任务 k 的个体
    snp=[];%创建一个空的 snp 向量，用于存放SNP的索引
    for i=1:HMS%个体生成
        l=randi([1,epi_dim]);
        if rand<=ln %类内部
            if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                snp(1)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
            else
                snp(1)=ceil(rand*SNPs);
            end

            for j=2:epi_dim %随机选择一个 SNP，然后继续选择其他 SNP，确保它们各不相同。
                if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                    snp(j)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                else
                    snp(j)=ceil(rand*SNPs);
                end

                while ismember(snp(j),snp(1:j-1)) %检查向量 snp 中第 j 个元素是否在 snp 的前 j-1 个元素中存在
                    if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                        snp(j)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                    else
                        snp(j)=ceil(rand*SNPs);
                    end

                end
            end
        else %类间初始化
            for l=1:epi_dim
                if ~isempty(clusterIndices{l})
                    snp(l)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                else
                    snp(l)=ceil(rand*SNPs);
                end
                while ismember(snp(l),snp(1:l-1)) %检查向量 snp 中第 j 个元素是否在 snp 的前 j-1 个元素中存在
                    if ~isempty(clusterIndices{l})
                        snp(l)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                    else
                        snp(l)=ceil(rand*SNPs);
                    end
                end
            end

        end
        %随机初始化了snp的取值 生成了7个
        %    snp(1)=ceil(rand*SNPs);%在每次循环迭代中，生成一个随机的 SNP 集合（snp），确保它们互不相同。
        %    for j=2:epi_dim %随机选择一个 SNP，然后继续选择其他 SNP，确保它们各不相同。
        %     snp(j)=ceil(rand*SNPs); %生成一个在区间 [1, SNPs] 内的随机整数
        %     while ismember(snp(j),snp(1:j-1)) %检查向量 snp 中第 j 个元素是否在 snp 的前 j-1 个元素中存在
        %        snp(j)=ceil(rand*SNPs);
        %  end
        %   end
        temp=snp;
        snp2=sort(snp(1:(k+1)));%将生成的 SNP 集合排序并存储在 snp2 中
        while ismember(snp,Task(k).X,'rows') %确保这个 SNP 集合没有在当前任务的 Task(k).X 中出现过
            j=ceil(rand*epi_dim);
            snp(j)=ceil(rand*SNPs);
            temp=snp;
            snp2=sort(snp(1:(k+1)));
        end

        X(i,:)=[snp2,snp(k+2 : epi_dim)];   %% 将 snp2 中的 SNP 存储在 X(i,:) 中，该矩阵用于存储有序的解
        HM(i,:)=temp;  %% 存储无序的 SNP 解。
        % 目标函数计算：
        %Fit(i,1) 到 Fit(i,4) 存储不同的目标函数值，可能是算法评价的多个指标
        [Fit(i,1),Fit(i,2),Fit(i,3),Fit(i,4)] = multi_criteriaEvaluationFuns2021(data(:,HM(i,1:dim)),State);
        NC = NC + 1;%NC 记录目标函数计算的次数
        snp=[];
    end
    maxFit = [maxFit; Fit];%将变量Fit的内容添加到maxFit的底部

    Task(k).X = X(1:HMS,:);% 在任务k中，X的前k列是有序的
    Task(k).HM = HM(1:HMS,:);% HM是无序的
    Task(k).Fit = Fit(1:HMS,:);% 存储的HM每一行的适应度值
    %这个循环迭代FitNum次，每次都将X的前EliteSize个元素分配给Elite(i).X，HM的前EliteSize个元素分配给Elite(i).HM，
    %以及Fit的前EliteSize个元素分配给Elite(i).Fit。Elite可能是一个结构体数组，用于存储优秀的解。
    for i = 1:FitNum %给精英子集赋值和任务task类似
        Elite(i).X = X(1:EliteSize,:);
        Elite(i).HM = HM(1:EliteSize,:);
        Elite(i).Fit = Fit(1:EliteSize,:);
    end
    %嵌套的循环。外层的两个循环迭代EliteSize次和FitNum次，内部的循环从EliteSize+1开始，直到Fit数组的长度。
    %在内部循环中，检查Elite(s).Fit(j,s)是否大于Fit(i,s)，如果是，就用Fit(i,:)替换Elite(s).Fit(j,:)，用X(i,:)替换Elite(s).X(j,:)，以及用HM(i,:)
    %替换Elite(s).HM(j,:)，然后通过break终止内部循环。这样做的目的是将更好的解复制到Elite结构中。
    for j = 1:EliteSize% 更新elite让其更优秀
        for s = 1:FitNum
            for i = EliteSize+1 : length(Fit(:,1))%fit的行数
                if Elite(s).Fit(j,s) > Fit(i,s)
                    Elite(s).Fit(j,:) = Fit(i,:);
                    Elite(s).X(j,:) = X(i,:);
                    Elite(s).HM(j,:) = HM(i,:);
                    break;
                end
            end
        end
    end

    Task(k).Elite = Elite;%将更新后的Elite结构分配给Task(k).Elite，将其与任务k相关联。


end







maxFit = max(maxFit);%计算变量maxFit的最大值，并将结果赋给maxFit。
% 对所有评价指标适应值进行归一化
%并且得到一些最佳的snp组合
for k = 1:K
    %对Task(k).Fit进行归一化
    Task(k).Fit = Task(k).Fit ./ maxFit;
    for i = 1:FitNum
        %对Task(k).Elite(i).Fit进行归一化
        Task(k).Elite(i).Fit = Task(k).Elite(i).Fit ./maxFit;
        %使用mink函数找到Task(k).Elite(i).Fit(:,i)中最小的bestNum个值，并将结果存储在Task(k).Elite(i).fbest中，下标存储在bestId中
        [Task(k).Elite(i).fbest, bestId] = mink(Task(k).Elite(i).Fit(:,i),bestNum);
        %根据这些最小值的索引，从Task(k).Elite(i).X中获取对应的最佳解，并存储在Task(k).Elite(i).Xbest中。
        Task(k).Elite(i).Xbest = Task(k).Elite(i).X(bestId,:);
    end
end


%先调用一个名为multi_criteriaEvaluationFuns2021的函数，传入参数data(:,CX)和State，并将返回的结果分别赋给s的前四个元素。接着，将s除以maxFit，实现对s的归一化。
[s(1),s(2),s(3),s(4)] =  multi_criteriaEvaluationFuns2021(data(:,CX),State);
s = s./maxFit;
%创建一个向量Dims，包含从2到epi_dim的整数。
Dims = [2:epi_dim];

LT=0;
%%-------------------------------------------------------------------------
tic;
while NC <= max_iter %一个嵌套的循环结构，用于实现一种优化算法
    for dim = Dims%表示不同任务
        k = dim - 1;  %% 从第k个任务群众探索     计算一个变量k，它的值是dim减去1
        Rs = rand;%生成一个随机数Rs，这个随机数的值在0到1之间。

        % if Rs < TP
        %% 迁移学习 在迁移学习的情况下，这段代码生成一个随机整数k0，其范围在1到K之间（K是一个变量）。
        %然后，它会检查k0是否等于之前计算的k，如果相等，就重新生成一个不同的k0，以确保k0不等于k。
        k0 = ceil(rand*K);
        while k0 == k
            k0 = ceil(rand*K);
        end
        %  end

        %将当前的k和dim分别存储在Ks和Ds中。
        Ks = k;
        Ds = dim;
        %分别初始化一个变量i为1，以及一个变量d，其值是一个随机整数，范围在1到FitNum之间。
        i=1;
        d = ceil(rand*FitNum);
        Xnew1=[];
        % 在这里加一个生成Xnew的算法
        
        try
            while i<= epi_dim%嵌套的while循环，它会一直执行，直到i的值小于或等于epi_dim。这个循环似乎是用于处理与优化组合相关的操作。
                %
                if Rs >= TP %% 从当前任务中进行 优化组合
                    if rand<HMCR%它检查一个随机数是否小于HMCR
                        a = ceil(rand*EliteSize);      %生成两个随机整数a和b，分别在1到EliteSize和1到epi_dim之间。然后，将Task(k).X(a, b)的值赋给Xnew(i)。
                        b = ceil(rand*epi_dim);
                        Xnew(i) = Task(k).X(a,b);
                        %  if ~isempty(clusterIndices{i})
                        %    Xnew(i) =clusterIndices{i}(ceil(rand*size(clusterIndices{i}, 1)));
                        %  end
                        if rand < PAR%检查另一个随机数是否小于PAR
                            sPar = ceil(rand*4);%生成一个随机整数sPar，其范围在1到4之间，然后生成两个随机整数b和c，分别在1到epi_dim和1到EliteSize之间。
                            b = ceil(rand*epi_dim);
                            c = ceil(rand*EliteSize);
                            while c == a%一个循环，它用于确保c不等于之前生成的a。如果c等于a，则重新生成一个不同的c
                                c = ceil(rand*EliteSize);
                            end
                            bs = ceil(rand*bestNum);%生成一个随机整数bs，其范围在1到bestNum之间。
                            switch sPar%根据sPar的值执行不同的操作。
                                case 1%如果sPar的值等于1，那么将Task(k).Elite(d).Xbest(bs, b)的值赋给Xnew(i)。
                                    Xnew(i) = Task(k).Elite(d).Xbest(bs,b);
                                case 2%如果sPar的值等于2，那么首先生成一个随机整数e，其范围在1到HMS之间。然后计算L，其值为Task(k).Elite(d).X(c, b)减去Task(k).X(e, i)。
                                    %最后，计算新的Xnew(i)，它等于原始值加上F乘以随机数乘以L，然后使用round函数四舍五入。
                                    e = ceil(rand*HMS);

                                    L = Task(k).Elite(d).X(c,b)-Task(k).X(e,i);


                                    Xnew(i) = round( Xnew(i) + F * rand * L);
                                    %检查Xnew(i)是否大于SNPs或小于1
                                    if Xnew(i) > SNPs || Xnew(i) <1
                                        Xnew(i) = ceil(rand*SNPs);
                                    end
                                case 3
                                    e = ceil(rand*HMS);%生成一个1到HMS（估计是某个数值）之间的随机整数。
                                    L = Task(k).Elite(d).Xbest(bs,b) - Task(k).X(e,b);%计算两个数组元素的差，其中Task(k).Elite(d).Xbest(bs,b)和Task(k).X(e,b)是从某个数据结构中获取的值。
                                    Xnew(i)=round(Xnew(i) + F * rand * L);%更新数组Xnew的第i个元素。F是一个常数，rand生成0到1之间的随机数。这一行代码基于L的值来更新Xnew(i)。
                                    %Xnew(i)=max(min(Xnew(i),SNPs),1);
                                    %如果Xnew(i)的值大于SNPs，则将Xnew(i)设置为SNPs减去一个随机数，该随机数的范围是0到min(10,max(L/10,1))。如果Xnew(i)的值小于1，将Xnew(i)设置为1加上一个随机数，
                                    %该随机数的范围也是0到min(10,max(L/10,1))。
                                      if Xnew(i) > SNPs || Xnew(i) <1 %检查 Xnew(i) 是否在有效范围内  在这里，如果 Xnew(i) 仍然超出范围，它会被修正。如果 Xnew(i) 大于 SNPs
                                        %，它将被减去一个随机数，否则如果 Xnew(i) 小于 1，它将被增加一个随机数。
                                        Xnew(i) = ceil(rand*SNPs);
                                    end
                                    if  Xnew(i) > SNPs
                                        Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    elseif Xnew(i) < 1
                                        Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    end
                            end
                        end
                        %上述条件不满足，就将Xnew(i)设置为1到SNPs之间的随机整数。
                    else
                        Xnew(i)=ceil(rand*SNPs);

                    end
                else %执行算法1，其他任务引导
                    if rand<HMCR %如果一个0到1之间的随机数小于HMCR，则执行以下代码块。
                        a = ceil(rand*EliteSize);%生成两个随机整数a和b，然后将Xnew(i)的值设置为从Task(k0).Elite(d).X数组中获取的一个特定值。
                        b = ceil(rand*epi_dim);

                        Xnew(i) = Task(k0).Elite(d).X(a,b);
                        if rand < PAR%一个0到1之间的随机数小于PAR，则执行以下代码块。
                            %生成三个随机整数sPar、b和c，然后b的值再次被重新赋值。while循环确保c不等于a
                            sPar = ceil(rand*3);
                            b = ceil(rand*epi_dim);
                            c = ceil(rand*EliteSize);
                            while c == a
                                c = ceil(rand*EliteSize);
                            end
                            bs = ceil(rand*bestNum);
                            switch sPar%基于sPar的值，选择不同的情况。如果sPar为1，将Xnew(i)的值设置为Task(k0).Elite(d).Xbest(bs,b)。
                                case 1
                                    Xnew(i) = Task(k0).Elite(d).Xbest(bs,b);
                                case 2
                                    %这行代码生成一个随机整数 e，范围在 1 到 HMS 之间，使用 rand 函数产生一个 0 到 1 之间的随机数，然后使用 ceil 函数向上取整。
                                    e = ceil(rand*HMS);
                                    %计算 L，这是两个数组元素之间的差，其中 Task(k0).Elite(d).X(bs,b) 和 Task(k0).X(e,i) 是数组中的特定元素。
                                    L = Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i);
                                    %这行代码更新数组 Xnew 中的第 i 个元素，将其增加一个随机数乘以 L 的结果，然后四舍五入。
                                    Xnew(i) = Xnew(i)+round(rand * (Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i)));
                                    %这行代码检查 Xnew(i) 是否大于 SNPs 或小于 1。
                                    if Xnew(i) > SNPs || Xnew(i) <1
                                        %如果 Xnew(i) 超出范围，它会被设置为在 1 和 SNPs 之间的随机整数。
                                        Xnew(i) = ceil(rand*SNPs);
                                    end


                                case 3
                                    %生成随机整数 d0，范围在 1 到 FitNum 之间。
                                    d0 = ceil(rand*FitNum);
                                    while d0 == d%这是一个循环，它确保 d0 不等于 d。
                                        d0 = ceil(rand*FitNum);
                                    end
                                    e = ceil(rand*HMS);%生成另一个随机整数 e，范围在 1 到 HMS 之间。
                                    L = Task(k0).Elite(d0).Xbest(bs,b) - Task(k0).X(e,b);%计算新的 L，类似于上面的计算。
                                    Xnew(i) = Xnew(i)+round( F * rand * L);%这行代码基于 L 的计算结果更新 Xnew(i)，并引入了一个额外的参数 F。
                                    if Xnew(i) > SNPs || Xnew(i) <1 %检查 Xnew(i) 是否在有效范围内  在这里，如果 Xnew(i) 仍然超出范围，它会被修正。如果 Xnew(i) 大于 SNPs
                                        %，它将被减去一个随机数，否则如果 Xnew(i) 小于 1，它将被增加一个随机数。
                                        Xnew(i) = ceil(rand*SNPs);
                                    end
                                    if  Xnew(i) > SNPs
                                        Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    elseif Xnew(i) < 1
                                        Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    end
                            end
                        end
                    else
                        Xnew(i)=ceil(rand*SNPs);%%如果没有满足条件的更新，那么 Xnew(i) 将被设置为一个随机整数，范围在 1 到 SNPs 之间。
                    end

                end
                %检查 i 是否等于 1 或者 Xnew(i) 是否不在数组 Xnew 的前 i-1 个元素中。如果条件为真，i 的值将递增，否则代码块结束。
                if i==1 || ~ismember(Xnew(i),Xnew(1:i-1))
                    i = i + 1;
                end

                if ( i-1 == epi_dim )
                    %将变量 Xnew 的值赋给变量 Xtemp。这创建了一个 Xtemp 的副本，以便稍后可以比较两者的值。
                    Xtemp=Xnew;
                    %将取 Xnew 中的前 dim 个元素，然后对它们进行排序，并将排序后的结果存储在 Xnew0 变量中。
                    Xnew0=sort(Xnew(1:dim));
                    %将 Xnew0 的值复制回 Xnew 中的前 dim 个元素，以确保它们按升序排列。
                    Xnew(1:dim) = Xnew0;
                    % 让xnew 不同于任务k中的snp
                    while ( ismember(Xnew(1:dim),Task(k).X(:,1:dim),'rows')  )
                        j=ceil(rand*dim);%生成一个随机整数 j，范围在1到 dim 之间。
                        r=ceil(rand*SNPs);%生成一个随机整数 r，范围在1到 SNPs 之间。
                        while ismember(r,Xnew)%检查 r 是否已经存在于 Xnew 中，如果存在，就重新生成一个随机整数 r，直到 r 不再存在于 Xnew 中。
                            r=ceil(rand*SNPs);
                        end
                        Xnew(j)=r;%将 Xnew 中的第 j 个元素设置为新生成的随机整数 r
                        Xtemp=Xnew;%再次将 Xnew 的值赋给 Xtemp，以便备份当前状态。
                        Xnew0=sort(Xnew(1:dim));%重新对 Xnew 中的前 dim 个元素进行排序。
                        Xnew(1:dim) = Xnew0;%将排序后的值赋给 Xnew 的前 dim 个元素，以确保它们按升序排列。
                    end

                    for b = 1:FitNum%这是一个 for 循环，变量 b 从1循环到 FitNum。
                        Xtemp=Xnew;%再次将 Xnew 的值赋给 Xtemp，以备份当前状态。
                        Xnew0=sort(Xnew(1:dim));%重新对 Xnew 中的前 dim 个元素进行排序。
                        Xnew(1:dim) = Xnew0;%将排序后的值赋给 Xnew 的前 dim 个元素，以确保它们按升序排列。
                        % 让xnew不同于任务k中的所有精英集合
                        while ( ismember(Xnew(1:dim),Task(k).Elite(b).X(:,1:dim),'rows')  )
                            j=ceil(rand*dim);%生成一个随机整数 j，范围在1到 dim 之间。
                            r=ceil(rand*SNPs);%生成一个随机整数 r，范围在1到 SNPs 之间。
                            while ismember(r,Xnew)%嵌套的循环，它检查 r 是否已经存在于 Xnew 中，如果存在，就重新生成
                                %一个随机整数 r，直到 r 不再存在于 Xnew 中。
                                r=ceil(rand*SNPs);
                            end
                            Xnew(j)=r;%将 Xnew 中的第 j 个元素设置为新生成的随机整数 r
                            Xtemp=Xnew;
                            Xnew0=sort(Xnew(1:dim));
                            Xnew(1:dim) = Xnew0;
                        end
                    end
                    if Rs < TP%Rs 和 TP 是可能是某些变量或参数，用于控制程序的流程。这段代码开始于一个条件判断，判断 Rs 是否小于 TP。
                        Xtemp=Xnew;%Xtemp 被赋值为 Xnew。
                        Xnew0=sort(Xnew(1:Ds));
                        Xnew(1:Ds) = Xnew0;%Xnew 中的前 Ds 个元素被排序，并赋值给 Xnew0。
                        %Xnew 的前 Ds 个元素是否存在于 Task(k0).X 中的某一行。这似乎是在检查 Xnew 的前 Ds 个元素是否与 Task(k0).X 的某些行相同。
                        while ( ismember(Xnew(1:Ds),Task(k0).X(:,1:Ds),'rows')  )
                            j=ceil(rand*Ds);%首先生成一个随机数 j，该值介于 1 和 Ds 之间
                            r=ceil(rand*SNPs);%然后生成一个随机数 r，介于 1 和 SNPs 之间。
                            while ismember(r,Xnew)%序检查 r 是否已经存在于 Xnew 中，如果是，则重新生成 r 直到找到一个不在 Xnew 中的值。
                                r=ceil(rand*SNPs);
                            end
                            %将 Xnew 中的第 j 个元素赋值为 r，然后再次排序 Xnew 的前 Ds 个元素，并将其赋值给 Xnew0。这似乎是为了确保 Xnew 的前 Ds 个元素保持排序状态。
                            Xnew(j)=r;
                            Xtemp=Xnew;
                            Xnew0=sort(Xnew(1:Ds));
                            Xnew(1:Ds) = Xnew0;
                        end
                        %使用 b 作为迭代变量，该循环重复以下操作 FitNum 次：
                        for b = 1:FitNum
                            Xtemp=Xnew;%Xtemp 被赋值为 Xnew。
                            Xnew0=sort(Xnew(1:Ds));%Xnew 中的前 Ds 个元素被排序，并赋值给 Xnew0
                            Xnew(1:Ds) = Xnew0;
                            %检查 Xnew 的前 Ds 个元素是否存在于 Task(k0).Elite(b).X 中的某一行。
                            while ( ismember(Xnew(1:Ds),Task(k0).Elite(b).X(:,1:Ds),'rows')  )
                                j=ceil(rand*Ds);%如果条件成立，再次生成随机数 j 和 r，并确保 r 不在 Xnew 中。
                                r=ceil(rand*SNPs);
                                while ismember(r,Xnew)
                                    r=ceil(rand*SNPs);
                                end
                                %将 Xnew 中的第 j 个元素赋值为 r，然后再次排序 Xnew 的前 Ds 个元素，并将其赋值给 Xnew0。
                                Xnew(j)=r;
                                Xtemp=Xnew;
                                Xnew0=sort(Xnew(1:Ds));
                                Xnew(1:Ds) = Xnew0;
                            end
                        end
                    end
                end
            end
        catch exception
            disp(exception);
        end
    

        %执行算法3
        if Rs < TP
            %这个函数的作用是对输入的数据进行多标准评估，返回四个评估分数，分别存储在 XnewScore 的前四个位置上。
            [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:Ds)),State);
            %将 XnewScore 中的每个分数除以 maxFit 的最大适应度值，用于进行归一化处理，将分数映射到 [0, 1] 范围内
            XnewScore = XnewScore ./ maxFit;
            %如果 Ds 等于 Epi_Dim，则将 Epi_Dim_FEs 增加 1。这段代码用于计算某个特定维度 Epi_Dim 下的函数评估次数。
            if Ds == Epi_Dim
                Epi_Dim_FEs = Epi_Dim_FEs + 1;
            end
            %将 NC 值加1。NC 变量可能用于追踪算法的迭代次数
            NC=NC+1;
            %对于范围在 1 到 HMS 之间的 i 进行循环操作，其中 HMS 是一个常数或变量
            for i = 1:HMS
                %查找在前 FitNum-1 个位置上，XnewScore 中的分数小于 Task(Ks).Fit(i,1:(FitNum-1)) 的索引，结果存储在 sn 中。这里 FitNum 是一个常数或变量，
                %Task(Ks).Fit(i,1:(FitNum-1)) 是一个长度为 FitNum-1 的向量。
                sn = find(XnewScore(1:(FitNum-1)) < Task(Ks).Fit(i,1:(FitNum-1)));
                %如果 sn 的长度大于 2，或者 XnewScore 的最后一个分数小于 Task(Ks).Fit(i,FitNum)，并且随机数 rand 小于 (1 - NC / max_iter)，则执行下面的操作。
               % if    length(sn) > 2 || (XnewScore(FitNum) < Task(Ks).Fit(i,FitNum) && rand < (1 - NC / max_iter))
                 if    length(sn) > 2 || (XnewScore(FitNum) < Task(Ks).Fit(i,FitNum) )
              
                    %% 将 Task(Ks) 中的第 i 行的 X 替换为 Xnew
                    Task(Ks).X(i,:) = Xnew;
                    %将 Task(Ks) 中的第 i 行的 HM 替换为 Xtemp。
                    Task(Ks).HM(i,:) = Xtemp;
                    %将 Task(Ks) 中的第 i 行的 Fit 替换为 XnewScore
                    Task(Ks).Fit(i,:) = XnewScore;
                    %                     fprintf('11\n');
                    break; % 只替换其中之一
                end
            end

            for i = 1:FitNum
                %fworst 是当前 Elite 集合中第 i 个目标函数的最差值，worstId 是该最差值在 Elite 集合中的索引
                [fworst,worstId] = max(Task(Ks).Elite(i).Fit(:,i));
                %如果当前 Elite 集合中的最差值大于 XnewScore(i)，则执行以下操作。
                if fworst > XnewScore(i)
                    %用新解 Xnew 替换 Elite 集合中最差解的决策变量。
                    Task(Ks).Elite(i).X(worstId,:) = Xnew;
                    %用新解 Xnew 替换 Elite 集合中最差解的其他相关信息
                    Task(Ks).Elite(i).HM(worstId,:) = Xtemp;
                    %用新解的评估分数 XnewScore 替换 Elite 集合中最差解的评估分数。
                    Task(Ks).Elite(i).Fit(worstId,:) = XnewScore;
                    for s = 1:bestNum%对于 bestNum 变量中的每个值，执行以下操作。
                        %如果当前目标函数的评估分数小于 Elite 集合中第 i 个目标函数的最佳值。
                        if XnewScore(i) < Task(Ks).Elite(i).fbest(s)
                            %将 Elite 集合中第 i 个目标函数的最佳值更新为当前的评估分数。
                            Task(Ks).Elite(i).fbest(s) = XnewScore(i);
                            % 将 Elite 集合中第 i 个目标函数的最佳解替换为当前解 Xnew。
                            Task(Ks).Elite(i).Xbest(s,:) = Xnew;
                            %% 扩展
                            if dim < epi_dim
                                %对新解 Xnew 的前 dim+1 个元素进行多标准评估，将结果存储在 Score 中。
                                [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State);
                                %将评估分数归一化到 [0, 1] 范围内。
                                Score = Score ./ maxFit;
                                %如果 dim+1 等于 Epi_Dim，则增加 Epi_Dim_FEs 计数器的值。
                                if dim+1 == Epi_Dim
                                    Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                end
                                %增加 NC 计数器的值。
                                NC = NC + 1;
                                %对于 FitNum 中的每个目标函数，执行以下操作。
                                for si = 1:FitNum
                                    %对于 EliteSize 中的每个值，执行以下操作。
                                    for sj = 1:EliteSize
                                        %如果当前维度 dim 下的 Elite 集合中的第 si 个目标函数的第 sj 个解在当前目标函数上的评估分数大于 Score(si)。
                                        if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                            %用排序后的新解 Xnew 替换 Elite 集合中的解。
                                            Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                            %用新解 Xnew 替换 Elite 集合中的相关信息。
                                            Task(dim).Elite(si).HM(sj,:) = Xnew;
                                            %用新的评估分数 Score 替换 Elite 集合中的评估分数。
                                            Task(dim).Elite(si).Fit(sj,:) = Score;
                                            break;
                                        end
                                    end
                                end
                            end


                            break;
                        end

                    end
                end
            end
        else
            %调用函数 multi_criteriaEvaluationFuns2021 对数据 data 的子集进行多标准评估，并将结果存储在 XnewScore 中，其中 Xnew 的前 dim 个元素被传递给该函数。
            [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim)),State);
            %将评估分数 XnewScore 归一化到 [0, 1] 范围内，通过将其除以 maxFit。
            XnewScore  = XnewScore ./ maxFit;
            %如果当前维度 dim 等于 Epi_Dim。
            if dim == Epi_Dim
                %增加 Epi_Dim_FEs 计数器的值。
                Epi_Dim_FEs = Epi_Dim_FEs + 1;
            end
            %增加 NC 计数器的值
            NC=NC+1;

            for i = 1:HMS   %对于循环索引 i 从 1 到 HMS。
                %查找满足条件的索引 sn，其中 XnewScore 的前 FitNum-1 个元素小于 Elite 集合中的第 i 个解的前 FitNum-1 个目标函数值。
                sn = find(XnewScore(1:(FitNum-1)) < Task(k).Fit(i,1:(FitNum-1)));
                %如果索引 sn 的长度大于 2 或者最后一个目标函数值 XnewScore(FitNum) 小于 Elite 集合中的第 i 个解的最后一个目标函数值，
                %同时随机数小于 (1 - NC / max_iter)。
                %if   length(sn) >2 || (XnewScore(FitNum) < Task(k).Fit(i,FitNum) &&  rand < (1 - NC / max_iter))
                 if   length(sn) >2 || (XnewScore(FitNum) < Task(k).Fit(i,FitNum) )
                
                    %% 用新解 Xnew 替换 Elite 集合中的第 i 个解的决策变量。
                    Task(k).X(i,:) = Xnew;
                    %用新解 Xtemp 替换 Elite 集合中的第 i 个解的其他相关信息。
                    Task(k).HM(i,:) = Xtemp;
                    %用新的评估分数 XnewScore 替换 Elite 集合中的第 i 个解的评估分数。
                    Task(k).Fit(i,:) = XnewScore;
                    %                     fprintf('22\n');
                    break; % 只替换其中之一
                end
            end


            %% 更新 不同评价标准 的精英集合
            for i = 1:FitNum%对于索引 i 从1到FitNum，进行循环。
                %找到 Elite 集合中的第 i 个子问题的目标函数中的最差值 fworst 和其对应的索引 worstId
                [fworst,worstId] = max(Task(k).Elite(i).Fit(:,i));
                %如果最差值 fworst 大于新解 Xnew 中的第 i 个目标函数值 XnewScore(i)
                if fworst > XnewScore(i)
                    %将新解 Xnew 替换 Elite 集合中的最差解的决策变量。
                    Task(k).Elite(i).X(worstId,:) = Xnew;
                    %将新解 Xtemp 替换 Elite 集合中的最差解的其他相关信息。
                    Task(k).Elite(i).HM(worstId,:) = Xtemp;
                    %将新的评估分数 XnewScore 替换 Elite 集合中的最差解的评估分数。
                    Task(k).Elite(i).Fit(worstId,:) = XnewScore;
                    for s = 1:bestNum%对于索引 s 从1到bestNum，进行循环。
                        %如果新解 Xnew 中的第 i 个目标函数值 XnewScore(i) 小于 Elite 集合中的第 i 个子问题的最佳值
                        if XnewScore(i) < Task(k).Elite(i).fbest(s)
                            %将最佳值 fbest(s) 更新为新的目标函数值 XnewScore(i)。
                            Task(k).Elite(i).fbest(s) = XnewScore(i);
                            %将最佳决策变量 Xbest(s, :) 更新为新解 Xnew。
                            Task(k).Elite(i).Xbest(s,:) = Xnew;
                            %% 如果当前维度 dim 小于 epi_dim。
                            if dim < epi_dim
                                %对数据的扩展子集进行多标准评估，将结果存储在 Score 中，其中 Xnew 的前 dim+1 个元素被传递给该函数。
                                [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State);
                                %将评估分数 Score 归一化到 [0, 1] 范围内，通过将其除以 maxFit。
                                Score = Score ./ maxFit;
                                %如果扩展后的维度 dim+1 等于 Epi_Dim。
                                if dim+1 == Epi_Dim
                                    %增加 Epi_Dim_FEs 计数器的值。
                                    Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                end
                                NC = NC + 1;%增加 NC 计数器的值。
                                for si = 1:FitNum
                                    %对于索引 sj 从1到EliteSize，进行循环。
                                    for sj = 1:EliteSize
                                        %如果扩展后的 Elite 集合中的第 si 个子问题的第 sj 个解的目标函数值大于扩展后的目标函数值
                                        if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                            %将扩展后的决策变量 Xnew 替换为 Elite 集合中的第 si 个子问题的第 sj 个解，并对其进行排序。
                                            Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                            %将扩展后的决策变量 Xnew 替换为 Elite 集合中的第 si 个子问题的第 sj 个解的其他相关信息。
                                            Task(dim).Elite(si).HM(sj,:) = Xnew;
                                            %将扩展后的评估分数 Score 替换为 Elite 集合中的第 si 个子问题的第 sj 个解的评估分数。
                                            Task(dim).Elite(si).Fit(sj,:) = Score;
                                            break;
                                        end
                                    end
                                end
                            end

                            break;
                        end

                    end
                end
            end

        end

        %%  The program is terminted if the Xnew is the solution.
        % cflag = 0;
        %for ci = 1:fdim

        %if ismember(CX(ci), Xnew) %检查 CX(ci) 是否存在于 Xnew 中  CX是致病模型

        %   cflag = cflag + 1;
        %end
        %end
        %if cflag == fdim
        %将 Xnew 的值赋给了一个特定的数组元素，涉及了多层的数据结构操作。
        %   Task(fdim-1).Elite(1).X(1,:) = Xnew;
        %将 XnewScore 的值赋给了一个特定的数组元素。
        %  Task(fdim-1).Elite(1).Fit(1,:) = XnewScore;
        %  flag = 1;
        %  break;
        %  end



    end
    %if flag == 1
    % break;
    %end

end

totaltime=toc;
%