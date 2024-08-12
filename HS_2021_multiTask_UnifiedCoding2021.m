%���ݼ� ������  ���ӽ����Ĳ���  ��Ⱥ��С ������  TP, PAR,F,HMCR �㷨�Ŀ��Ʋ��� ����ɾ����ʼ��
function [Task,NC,flag,Epi_Dim_FEs,ci] = HS_2021_multiTask_UnifiedCoding2021(data,epi_dim,s,HMS,max_iter,CX, TP, PAR,F,HMCR,sample,ln)
% ������ͳһ����
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
K = epi_dim-1 ; % ������ 2-->epi_dim
FitNum = 4;
bestNum = epi_dim;
%% ---------------------------------------------------------------
EliteSize =min(10*epi_dim,ceil(HMS/5)); %��Ӣ�Ӽ��Ĵ�С
% CX
NC=0;
SNPs=n-1;  %% ��SNP����
flag = -1;
maxFit = [];



clusterIndices=cengci(data,epi_dim,sample);

for k = 1:K %ѭ��ÿһ������ �������ֶ�����  ��ʼ�����task ��elite 5��������㵽��7��
    %% ��ʼ��

    dim = k+1; %����ĳһ������
    Task(k).X=zeros(HMS,epi_dim);%��ʼ��һ���յ� Task(k).X �������ڴ������ k �ĸ���
    snp=[];%����һ���յ� snp ���������ڴ��SNP������
    for i=1:HMS%��������
        l=randi([1,epi_dim]);
        if rand<=ln %���ڲ�
            if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                snp(1)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
            else
                snp(1)=ceil(rand*SNPs);
            end

            for j=2:epi_dim %���ѡ��һ�� SNP��Ȼ�����ѡ������ SNP��ȷ�����Ǹ�����ͬ��
                if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                    snp(j)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                else
                    snp(j)=ceil(rand*SNPs);
                end

                while ismember(snp(j),snp(1:j-1)) %������� snp �е� j ��Ԫ���Ƿ��� snp ��ǰ j-1 ��Ԫ���д���
                    if length(clusterIndices{l})>=epi_dim&&~isempty(clusterIndices{l})
                        snp(j)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                    else
                        snp(j)=ceil(rand*SNPs);
                    end

                end
            end
        else %����ʼ��
            for l=1:epi_dim
                if ~isempty(clusterIndices{l})
                    snp(l)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                else
                    snp(l)=ceil(rand*SNPs);
                end
                while ismember(snp(l),snp(1:l-1)) %������� snp �е� j ��Ԫ���Ƿ��� snp ��ǰ j-1 ��Ԫ���д���
                    if ~isempty(clusterIndices{l})
                        snp(l)=clusterIndices{l}(ceil(rand*size(clusterIndices{l}, 1)));
                    else
                        snp(l)=ceil(rand*SNPs);
                    end
                end
            end

        end
        %�����ʼ����snp��ȡֵ ������7��
        %    snp(1)=ceil(rand*SNPs);%��ÿ��ѭ�������У�����һ������� SNP ���ϣ�snp����ȷ�����ǻ�����ͬ��
        %    for j=2:epi_dim %���ѡ��һ�� SNP��Ȼ�����ѡ������ SNP��ȷ�����Ǹ�����ͬ��
        %     snp(j)=ceil(rand*SNPs); %����һ�������� [1, SNPs] �ڵ��������
        %     while ismember(snp(j),snp(1:j-1)) %������� snp �е� j ��Ԫ���Ƿ��� snp ��ǰ j-1 ��Ԫ���д���
        %        snp(j)=ceil(rand*SNPs);
        %  end
        %   end
        temp=snp;
        snp2=sort(snp(1:(k+1)));%�����ɵ� SNP �������򲢴洢�� snp2 ��
        while ismember(snp,Task(k).X,'rows') %ȷ����� SNP ����û���ڵ�ǰ����� Task(k).X �г��ֹ�
            j=ceil(rand*epi_dim);
            snp(j)=ceil(rand*SNPs);
            temp=snp;
            snp2=sort(snp(1:(k+1)));
        end

        X(i,:)=[snp2,snp(k+2 : epi_dim)];   %% �� snp2 �е� SNP �洢�� X(i,:) �У��þ������ڴ洢����Ľ�
        HM(i,:)=temp;  %% �洢����� SNP �⡣
        % Ŀ�꺯�����㣺
        %Fit(i,1) �� Fit(i,4) �洢��ͬ��Ŀ�꺯��ֵ���������㷨���۵Ķ��ָ��
        [Fit(i,1),Fit(i,2),Fit(i,3),Fit(i,4)] = multi_criteriaEvaluationFuns2021(data(:,HM(i,1:dim)),State);
        NC = NC + 1;%NC ��¼Ŀ�꺯������Ĵ���
        snp=[];
    end
    maxFit = [maxFit; Fit];%������Fit��������ӵ�maxFit�ĵײ�

    Task(k).X = X(1:HMS,:);% ������k�У�X��ǰk���������
    Task(k).HM = HM(1:HMS,:);% HM�������
    Task(k).Fit = Fit(1:HMS,:);% �洢��HMÿһ�е���Ӧ��ֵ
    %���ѭ������FitNum�Σ�ÿ�ζ���X��ǰEliteSize��Ԫ�ط����Elite(i).X��HM��ǰEliteSize��Ԫ�ط����Elite(i).HM��
    %�Լ�Fit��ǰEliteSize��Ԫ�ط����Elite(i).Fit��Elite������һ���ṹ�����飬���ڴ洢����Ľ⡣
    for i = 1:FitNum %����Ӣ�Ӽ���ֵ������task����
        Elite(i).X = X(1:EliteSize,:);
        Elite(i).HM = HM(1:EliteSize,:);
        Elite(i).Fit = Fit(1:EliteSize,:);
    end
    %Ƕ�׵�ѭ������������ѭ������EliteSize�κ�FitNum�Σ��ڲ���ѭ����EliteSize+1��ʼ��ֱ��Fit����ĳ��ȡ�
    %���ڲ�ѭ���У����Elite(s).Fit(j,s)�Ƿ����Fit(i,s)������ǣ�����Fit(i,:)�滻Elite(s).Fit(j,:)����X(i,:)�滻Elite(s).X(j,:)���Լ���HM(i,:)
    %�滻Elite(s).HM(j,:)��Ȼ��ͨ��break��ֹ�ڲ�ѭ������������Ŀ���ǽ����õĽ⸴�Ƶ�Elite�ṹ�С�
    for j = 1:EliteSize% ����elite���������
        for s = 1:FitNum
            for i = EliteSize+1 : length(Fit(:,1))%fit������
                if Elite(s).Fit(j,s) > Fit(i,s)
                    Elite(s).Fit(j,:) = Fit(i,:);
                    Elite(s).X(j,:) = X(i,:);
                    Elite(s).HM(j,:) = HM(i,:);
                    break;
                end
            end
        end
    end

    Task(k).Elite = Elite;%�����º��Elite�ṹ�����Task(k).Elite������������k�������


end







maxFit = max(maxFit);%�������maxFit�����ֵ�������������maxFit��
% ����������ָ����Ӧֵ���й�һ��
%���ҵõ�һЩ��ѵ�snp���
for k = 1:K
    %��Task(k).Fit���й�һ��
    Task(k).Fit = Task(k).Fit ./ maxFit;
    for i = 1:FitNum
        %��Task(k).Elite(i).Fit���й�һ��
        Task(k).Elite(i).Fit = Task(k).Elite(i).Fit ./maxFit;
        %ʹ��mink�����ҵ�Task(k).Elite(i).Fit(:,i)����С��bestNum��ֵ����������洢��Task(k).Elite(i).fbest�У��±�洢��bestId��
        [Task(k).Elite(i).fbest, bestId] = mink(Task(k).Elite(i).Fit(:,i),bestNum);
        %������Щ��Сֵ����������Task(k).Elite(i).X�л�ȡ��Ӧ����ѽ⣬���洢��Task(k).Elite(i).Xbest�С�
        Task(k).Elite(i).Xbest = Task(k).Elite(i).X(bestId,:);
    end
end


%�ȵ���һ����Ϊmulti_criteriaEvaluationFuns2021�ĺ������������data(:,CX)��State���������صĽ���ֱ𸳸�s��ǰ�ĸ�Ԫ�ء����ţ���s����maxFit��ʵ�ֶ�s�Ĺ�һ����
[s(1),s(2),s(3),s(4)] =  multi_criteriaEvaluationFuns2021(data(:,CX),State);
s = s./maxFit;
%����һ������Dims��������2��epi_dim��������
Dims = [2:epi_dim];

LT=0;
%%-------------------------------------------------------------------------
tic;
while NC <= max_iter %һ��Ƕ�׵�ѭ���ṹ������ʵ��һ���Ż��㷨
    for dim = Dims%��ʾ��ͬ����
        k = dim - 1;  %% �ӵ�k������Ⱥ��̽��     ����һ������k������ֵ��dim��ȥ1
        Rs = rand;%����һ�������Rs������������ֵ��0��1֮�䡣

        % if Rs < TP
        %% Ǩ��ѧϰ ��Ǩ��ѧϰ������£���δ�������һ���������k0���䷶Χ��1��K֮�䣨K��һ����������
        %Ȼ��������k0�Ƿ����֮ǰ�����k�������ȣ�����������һ����ͬ��k0����ȷ��k0������k��
        k0 = ceil(rand*K);
        while k0 == k
            k0 = ceil(rand*K);
        end
        %  end

        %����ǰ��k��dim�ֱ�洢��Ks��Ds�С�
        Ks = k;
        Ds = dim;
        %�ֱ��ʼ��һ������iΪ1���Լ�һ������d����ֵ��һ�������������Χ��1��FitNum֮�䡣
        i=1;
        d = ceil(rand*FitNum);
        Xnew1=[];
        % �������һ������Xnew���㷨
        
        try
            while i<= epi_dim%Ƕ�׵�whileѭ��������һֱִ�У�ֱ��i��ֵС�ڻ����epi_dim�����ѭ���ƺ������ڴ������Ż������صĲ�����
                %
                if Rs >= TP %% �ӵ�ǰ�����н��� �Ż����
                    if rand<HMCR%�����һ��������Ƿ�С��HMCR
                        a = ceil(rand*EliteSize);      %���������������a��b���ֱ���1��EliteSize��1��epi_dim֮�䡣Ȼ�󣬽�Task(k).X(a, b)��ֵ����Xnew(i)��
                        b = ceil(rand*epi_dim);
                        Xnew(i) = Task(k).X(a,b);
                        %  if ~isempty(clusterIndices{i})
                        %    Xnew(i) =clusterIndices{i}(ceil(rand*size(clusterIndices{i}, 1)));
                        %  end
                        if rand < PAR%�����һ��������Ƿ�С��PAR
                            sPar = ceil(rand*4);%����һ���������sPar���䷶Χ��1��4֮�䣬Ȼ�����������������b��c���ֱ���1��epi_dim��1��EliteSize֮�䡣
                            b = ceil(rand*epi_dim);
                            c = ceil(rand*EliteSize);
                            while c == a%һ��ѭ����������ȷ��c������֮ǰ���ɵ�a�����c����a������������һ����ͬ��c
                                c = ceil(rand*EliteSize);
                            end
                            bs = ceil(rand*bestNum);%����һ���������bs���䷶Χ��1��bestNum֮�䡣
                            switch sPar%����sPar��ִֵ�в�ͬ�Ĳ�����
                                case 1%���sPar��ֵ����1����ô��Task(k).Elite(d).Xbest(bs, b)��ֵ����Xnew(i)��
                                    Xnew(i) = Task(k).Elite(d).Xbest(bs,b);
                                case 2%���sPar��ֵ����2����ô��������һ���������e���䷶Χ��1��HMS֮�䡣Ȼ�����L����ֵΪTask(k).Elite(d).X(c, b)��ȥTask(k).X(e, i)��
                                    %��󣬼����µ�Xnew(i)��������ԭʼֵ����F�������������L��Ȼ��ʹ��round�����������롣
                                    e = ceil(rand*HMS);

                                    L = Task(k).Elite(d).X(c,b)-Task(k).X(e,i);


                                    Xnew(i) = round( Xnew(i) + F * rand * L);
                                    %���Xnew(i)�Ƿ����SNPs��С��1
                                    if Xnew(i) > SNPs || Xnew(i) <1
                                        Xnew(i) = ceil(rand*SNPs);
                                    end
                                case 3
                                    e = ceil(rand*HMS);%����һ��1��HMS��������ĳ����ֵ��֮������������
                                    L = Task(k).Elite(d).Xbest(bs,b) - Task(k).X(e,b);%������������Ԫ�صĲ����Task(k).Elite(d).Xbest(bs,b)��Task(k).X(e,b)�Ǵ�ĳ�����ݽṹ�л�ȡ��ֵ��
                                    Xnew(i)=round(Xnew(i) + F * rand * L);%��������Xnew�ĵ�i��Ԫ�ء�F��һ��������rand����0��1֮������������һ�д������L��ֵ������Xnew(i)��
                                    %Xnew(i)=max(min(Xnew(i),SNPs),1);
                                    %���Xnew(i)��ֵ����SNPs����Xnew(i)����ΪSNPs��ȥһ�����������������ķ�Χ��0��min(10,max(L/10,1))�����Xnew(i)��ֵС��1����Xnew(i)����Ϊ1����һ���������
                                    %��������ķ�ΧҲ��0��min(10,max(L/10,1))��
                                      if Xnew(i) > SNPs || Xnew(i) <1 %��� Xnew(i) �Ƿ�����Ч��Χ��  �������� Xnew(i) ��Ȼ������Χ�����ᱻ��������� Xnew(i) ���� SNPs
                                        %����������ȥһ���������������� Xnew(i) С�� 1������������һ���������
                                        Xnew(i) = ceil(rand*SNPs);
                                    end
                                    if  Xnew(i) > SNPs
                                        Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    elseif Xnew(i) < 1
                                        Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                    end
                            end
                        end
                        %�������������㣬�ͽ�Xnew(i)����Ϊ1��SNPs֮������������
                    else
                        Xnew(i)=ceil(rand*SNPs);

                    end
                else %ִ���㷨1��������������
                    if rand<HMCR %���һ��0��1֮��������С��HMCR����ִ�����´���顣
                        a = ceil(rand*EliteSize);%���������������a��b��Ȼ��Xnew(i)��ֵ����Ϊ��Task(k0).Elite(d).X�����л�ȡ��һ���ض�ֵ��
                        b = ceil(rand*epi_dim);

                        Xnew(i) = Task(k0).Elite(d).X(a,b);
                        if rand < PAR%һ��0��1֮��������С��PAR����ִ�����´���顣
                            %���������������sPar��b��c��Ȼ��b��ֵ�ٴα����¸�ֵ��whileѭ��ȷ��c������a
                            sPar = ceil(rand*3);
                            b = ceil(rand*epi_dim);
                            c = ceil(rand*EliteSize);
                            while c == a
                                c = ceil(rand*EliteSize);
                            end
                            bs = ceil(rand*bestNum);
                            switch sPar%����sPar��ֵ��ѡ��ͬ����������sParΪ1����Xnew(i)��ֵ����ΪTask(k0).Elite(d).Xbest(bs,b)��
                                case 1
                                    Xnew(i) = Task(k0).Elite(d).Xbest(bs,b);
                                case 2
                                    %���д�������һ��������� e����Χ�� 1 �� HMS ֮�䣬ʹ�� rand ��������һ�� 0 �� 1 ֮����������Ȼ��ʹ�� ceil ��������ȡ����
                                    e = ceil(rand*HMS);
                                    %���� L��������������Ԫ��֮��Ĳ���� Task(k0).Elite(d).X(bs,b) �� Task(k0).X(e,i) �������е��ض�Ԫ�ء�
                                    L = Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i);
                                    %���д���������� Xnew �еĵ� i ��Ԫ�أ���������һ����������� L �Ľ����Ȼ���������롣
                                    Xnew(i) = Xnew(i)+round(rand * (Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i)));
                                    %���д����� Xnew(i) �Ƿ���� SNPs ��С�� 1��
                                    if Xnew(i) > SNPs || Xnew(i) <1
                                        %��� Xnew(i) ������Χ�����ᱻ����Ϊ�� 1 �� SNPs ֮������������
                                        Xnew(i) = ceil(rand*SNPs);
                                    end


                                case 3
                                    %����������� d0����Χ�� 1 �� FitNum ֮�䡣
                                    d0 = ceil(rand*FitNum);
                                    while d0 == d%����һ��ѭ������ȷ�� d0 ������ d��
                                        d0 = ceil(rand*FitNum);
                                    end
                                    e = ceil(rand*HMS);%������һ��������� e����Χ�� 1 �� HMS ֮�䡣
                                    L = Task(k0).Elite(d0).Xbest(bs,b) - Task(k0).X(e,b);%�����µ� L������������ļ��㡣
                                    Xnew(i) = Xnew(i)+round( F * rand * L);%���д������ L �ļ��������� Xnew(i)����������һ������Ĳ��� F��
                                    if Xnew(i) > SNPs || Xnew(i) <1 %��� Xnew(i) �Ƿ�����Ч��Χ��  �������� Xnew(i) ��Ȼ������Χ�����ᱻ��������� Xnew(i) ���� SNPs
                                        %����������ȥһ���������������� Xnew(i) С�� 1������������һ���������
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
                        Xnew(i)=ceil(rand*SNPs);%%���û�����������ĸ��£���ô Xnew(i) ��������Ϊһ�������������Χ�� 1 �� SNPs ֮�䡣
                    end

                end
                %��� i �Ƿ���� 1 ���� Xnew(i) �Ƿ������� Xnew ��ǰ i-1 ��Ԫ���С��������Ϊ�棬i ��ֵ���������������������
                if i==1 || ~ismember(Xnew(i),Xnew(1:i-1))
                    i = i + 1;
                end

                if ( i-1 == epi_dim )
                    %������ Xnew ��ֵ�������� Xtemp���ⴴ����һ�� Xtemp �ĸ������Ա��Ժ���ԱȽ����ߵ�ֵ��
                    Xtemp=Xnew;
                    %��ȡ Xnew �е�ǰ dim ��Ԫ�أ�Ȼ������ǽ������򣬲��������Ľ���洢�� Xnew0 �����С�
                    Xnew0=sort(Xnew(1:dim));
                    %�� Xnew0 ��ֵ���ƻ� Xnew �е�ǰ dim ��Ԫ�أ���ȷ�����ǰ��������С�
                    Xnew(1:dim) = Xnew0;
                    % ��xnew ��ͬ������k�е�snp
                    while ( ismember(Xnew(1:dim),Task(k).X(:,1:dim),'rows')  )
                        j=ceil(rand*dim);%����һ��������� j����Χ��1�� dim ֮�䡣
                        r=ceil(rand*SNPs);%����һ��������� r����Χ��1�� SNPs ֮�䡣
                        while ismember(r,Xnew)%��� r �Ƿ��Ѿ������� Xnew �У�������ڣ�����������һ��������� r��ֱ�� r ���ٴ����� Xnew �С�
                            r=ceil(rand*SNPs);
                        end
                        Xnew(j)=r;%�� Xnew �еĵ� j ��Ԫ������Ϊ�����ɵ�������� r
                        Xtemp=Xnew;%�ٴν� Xnew ��ֵ���� Xtemp���Ա㱸�ݵ�ǰ״̬��
                        Xnew0=sort(Xnew(1:dim));%���¶� Xnew �е�ǰ dim ��Ԫ�ؽ�������
                        Xnew(1:dim) = Xnew0;%��������ֵ���� Xnew ��ǰ dim ��Ԫ�أ���ȷ�����ǰ��������С�
                    end

                    for b = 1:FitNum%����һ�� for ѭ�������� b ��1ѭ���� FitNum��
                        Xtemp=Xnew;%�ٴν� Xnew ��ֵ���� Xtemp���Ա��ݵ�ǰ״̬��
                        Xnew0=sort(Xnew(1:dim));%���¶� Xnew �е�ǰ dim ��Ԫ�ؽ�������
                        Xnew(1:dim) = Xnew0;%��������ֵ���� Xnew ��ǰ dim ��Ԫ�أ���ȷ�����ǰ��������С�
                        % ��xnew��ͬ������k�е����о�Ӣ����
                        while ( ismember(Xnew(1:dim),Task(k).Elite(b).X(:,1:dim),'rows')  )
                            j=ceil(rand*dim);%����һ��������� j����Χ��1�� dim ֮�䡣
                            r=ceil(rand*SNPs);%����һ��������� r����Χ��1�� SNPs ֮�䡣
                            while ismember(r,Xnew)%Ƕ�׵�ѭ��������� r �Ƿ��Ѿ������� Xnew �У�������ڣ�����������
                                %һ��������� r��ֱ�� r ���ٴ����� Xnew �С�
                                r=ceil(rand*SNPs);
                            end
                            Xnew(j)=r;%�� Xnew �еĵ� j ��Ԫ������Ϊ�����ɵ�������� r
                            Xtemp=Xnew;
                            Xnew0=sort(Xnew(1:dim));
                            Xnew(1:dim) = Xnew0;
                        end
                    end
                    if Rs < TP%Rs �� TP �ǿ�����ĳЩ��������������ڿ��Ƴ�������̡���δ��뿪ʼ��һ�������жϣ��ж� Rs �Ƿ�С�� TP��
                        Xtemp=Xnew;%Xtemp ����ֵΪ Xnew��
                        Xnew0=sort(Xnew(1:Ds));
                        Xnew(1:Ds) = Xnew0;%Xnew �е�ǰ Ds ��Ԫ�ر����򣬲���ֵ�� Xnew0��
                        %Xnew ��ǰ Ds ��Ԫ���Ƿ������ Task(k0).X �е�ĳһ�С����ƺ����ڼ�� Xnew ��ǰ Ds ��Ԫ���Ƿ��� Task(k0).X ��ĳЩ����ͬ��
                        while ( ismember(Xnew(1:Ds),Task(k0).X(:,1:Ds),'rows')  )
                            j=ceil(rand*Ds);%��������һ������� j����ֵ���� 1 �� Ds ֮��
                            r=ceil(rand*SNPs);%Ȼ������һ������� r������ 1 �� SNPs ֮�䡣
                            while ismember(r,Xnew)%���� r �Ƿ��Ѿ������� Xnew �У�����ǣ����������� r ֱ���ҵ�һ������ Xnew �е�ֵ��
                                r=ceil(rand*SNPs);
                            end
                            %�� Xnew �еĵ� j ��Ԫ�ظ�ֵΪ r��Ȼ���ٴ����� Xnew ��ǰ Ds ��Ԫ�أ������丳ֵ�� Xnew0�����ƺ���Ϊ��ȷ�� Xnew ��ǰ Ds ��Ԫ�ر�������״̬��
                            Xnew(j)=r;
                            Xtemp=Xnew;
                            Xnew0=sort(Xnew(1:Ds));
                            Xnew(1:Ds) = Xnew0;
                        end
                        %ʹ�� b ��Ϊ������������ѭ���ظ����²��� FitNum �Σ�
                        for b = 1:FitNum
                            Xtemp=Xnew;%Xtemp ����ֵΪ Xnew��
                            Xnew0=sort(Xnew(1:Ds));%Xnew �е�ǰ Ds ��Ԫ�ر����򣬲���ֵ�� Xnew0
                            Xnew(1:Ds) = Xnew0;
                            %��� Xnew ��ǰ Ds ��Ԫ���Ƿ������ Task(k0).Elite(b).X �е�ĳһ�С�
                            while ( ismember(Xnew(1:Ds),Task(k0).Elite(b).X(:,1:Ds),'rows')  )
                                j=ceil(rand*Ds);%��������������ٴ���������� j �� r����ȷ�� r ���� Xnew �С�
                                r=ceil(rand*SNPs);
                                while ismember(r,Xnew)
                                    r=ceil(rand*SNPs);
                                end
                                %�� Xnew �еĵ� j ��Ԫ�ظ�ֵΪ r��Ȼ���ٴ����� Xnew ��ǰ Ds ��Ԫ�أ������丳ֵ�� Xnew0��
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
    

        %ִ���㷨3
        if Rs < TP
            %��������������Ƕ���������ݽ��ж��׼�����������ĸ������������ֱ�洢�� XnewScore ��ǰ�ĸ�λ���ϡ�
            [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:Ds)),State);
            %�� XnewScore �е�ÿ���������� maxFit �������Ӧ��ֵ�����ڽ��й�һ������������ӳ�䵽 [0, 1] ��Χ��
            XnewScore = XnewScore ./ maxFit;
            %��� Ds ���� Epi_Dim���� Epi_Dim_FEs ���� 1����δ������ڼ���ĳ���ض�ά�� Epi_Dim �µĺ�������������
            if Ds == Epi_Dim
                Epi_Dim_FEs = Epi_Dim_FEs + 1;
            end
            %�� NC ֵ��1��NC ������������׷���㷨�ĵ�������
            NC=NC+1;
            %���ڷ�Χ�� 1 �� HMS ֮��� i ����ѭ������������ HMS ��һ�����������
            for i = 1:HMS
                %������ǰ FitNum-1 ��λ���ϣ�XnewScore �еķ���С�� Task(Ks).Fit(i,1:(FitNum-1)) ������������洢�� sn �С����� FitNum ��һ�������������
                %Task(Ks).Fit(i,1:(FitNum-1)) ��һ������Ϊ FitNum-1 ��������
                sn = find(XnewScore(1:(FitNum-1)) < Task(Ks).Fit(i,1:(FitNum-1)));
                %��� sn �ĳ��ȴ��� 2������ XnewScore �����һ������С�� Task(Ks).Fit(i,FitNum)����������� rand С�� (1 - NC / max_iter)����ִ������Ĳ�����
               % if    length(sn) > 2 || (XnewScore(FitNum) < Task(Ks).Fit(i,FitNum) && rand < (1 - NC / max_iter))
                 if    length(sn) > 2 || (XnewScore(FitNum) < Task(Ks).Fit(i,FitNum) )
              
                    %% �� Task(Ks) �еĵ� i �е� X �滻Ϊ Xnew
                    Task(Ks).X(i,:) = Xnew;
                    %�� Task(Ks) �еĵ� i �е� HM �滻Ϊ Xtemp��
                    Task(Ks).HM(i,:) = Xtemp;
                    %�� Task(Ks) �еĵ� i �е� Fit �滻Ϊ XnewScore
                    Task(Ks).Fit(i,:) = XnewScore;
                    %                     fprintf('11\n');
                    break; % ֻ�滻����֮һ
                end
            end

            for i = 1:FitNum
                %fworst �ǵ�ǰ Elite �����е� i ��Ŀ�꺯�������ֵ��worstId �Ǹ����ֵ�� Elite �����е�����
                [fworst,worstId] = max(Task(Ks).Elite(i).Fit(:,i));
                %�����ǰ Elite �����е����ֵ���� XnewScore(i)����ִ�����²�����
                if fworst > XnewScore(i)
                    %���½� Xnew �滻 Elite ����������ľ��߱�����
                    Task(Ks).Elite(i).X(worstId,:) = Xnew;
                    %���½� Xnew �滻 Elite ��������������������Ϣ
                    Task(Ks).Elite(i).HM(worstId,:) = Xtemp;
                    %���½���������� XnewScore �滻 Elite ���������������������
                    Task(Ks).Elite(i).Fit(worstId,:) = XnewScore;
                    for s = 1:bestNum%���� bestNum �����е�ÿ��ֵ��ִ�����²�����
                        %�����ǰĿ�꺯������������С�� Elite �����е� i ��Ŀ�꺯�������ֵ��
                        if XnewScore(i) < Task(Ks).Elite(i).fbest(s)
                            %�� Elite �����е� i ��Ŀ�꺯�������ֵ����Ϊ��ǰ������������
                            Task(Ks).Elite(i).fbest(s) = XnewScore(i);
                            % �� Elite �����е� i ��Ŀ�꺯������ѽ��滻Ϊ��ǰ�� Xnew��
                            Task(Ks).Elite(i).Xbest(s,:) = Xnew;
                            %% ��չ
                            if dim < epi_dim
                                %���½� Xnew ��ǰ dim+1 ��Ԫ�ؽ��ж��׼������������洢�� Score �С�
                                [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State);
                                %������������һ���� [0, 1] ��Χ�ڡ�
                                Score = Score ./ maxFit;
                                %��� dim+1 ���� Epi_Dim�������� Epi_Dim_FEs ��������ֵ��
                                if dim+1 == Epi_Dim
                                    Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                end
                                %���� NC ��������ֵ��
                                NC = NC + 1;
                                %���� FitNum �е�ÿ��Ŀ�꺯����ִ�����²�����
                                for si = 1:FitNum
                                    %���� EliteSize �е�ÿ��ֵ��ִ�����²�����
                                    for sj = 1:EliteSize
                                        %�����ǰά�� dim �µ� Elite �����еĵ� si ��Ŀ�꺯���ĵ� sj �����ڵ�ǰĿ�꺯���ϵ������������� Score(si)��
                                        if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                            %���������½� Xnew �滻 Elite �����еĽ⡣
                                            Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                            %���½� Xnew �滻 Elite �����е������Ϣ��
                                            Task(dim).Elite(si).HM(sj,:) = Xnew;
                                            %���µ��������� Score �滻 Elite �����е�����������
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
            %���ú��� multi_criteriaEvaluationFuns2021 ������ data ���Ӽ����ж��׼��������������洢�� XnewScore �У����� Xnew ��ǰ dim ��Ԫ�ر����ݸ��ú�����
            [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim)),State);
            %���������� XnewScore ��һ���� [0, 1] ��Χ�ڣ�ͨ��������� maxFit��
            XnewScore  = XnewScore ./ maxFit;
            %�����ǰά�� dim ���� Epi_Dim��
            if dim == Epi_Dim
                %���� Epi_Dim_FEs ��������ֵ��
                Epi_Dim_FEs = Epi_Dim_FEs + 1;
            end
            %���� NC ��������ֵ
            NC=NC+1;

            for i = 1:HMS   %����ѭ������ i �� 1 �� HMS��
                %������������������ sn������ XnewScore ��ǰ FitNum-1 ��Ԫ��С�� Elite �����еĵ� i �����ǰ FitNum-1 ��Ŀ�꺯��ֵ��
                sn = find(XnewScore(1:(FitNum-1)) < Task(k).Fit(i,1:(FitNum-1)));
                %������� sn �ĳ��ȴ��� 2 �������һ��Ŀ�꺯��ֵ XnewScore(FitNum) С�� Elite �����еĵ� i ��������һ��Ŀ�꺯��ֵ��
                %ͬʱ�����С�� (1 - NC / max_iter)��
                %if   length(sn) >2 || (XnewScore(FitNum) < Task(k).Fit(i,FitNum) &&  rand < (1 - NC / max_iter))
                 if   length(sn) >2 || (XnewScore(FitNum) < Task(k).Fit(i,FitNum) )
                
                    %% ���½� Xnew �滻 Elite �����еĵ� i ����ľ��߱�����
                    Task(k).X(i,:) = Xnew;
                    %���½� Xtemp �滻 Elite �����еĵ� i ��������������Ϣ��
                    Task(k).HM(i,:) = Xtemp;
                    %���µ��������� XnewScore �滻 Elite �����еĵ� i ���������������
                    Task(k).Fit(i,:) = XnewScore;
                    %                     fprintf('22\n');
                    break; % ֻ�滻����֮һ
                end
            end


            %% ���� ��ͬ���۱�׼ �ľ�Ӣ����
            for i = 1:FitNum%�������� i ��1��FitNum������ѭ����
                %�ҵ� Elite �����еĵ� i ���������Ŀ�꺯���е����ֵ fworst �����Ӧ������ worstId
                [fworst,worstId] = max(Task(k).Elite(i).Fit(:,i));
                %������ֵ fworst �����½� Xnew �еĵ� i ��Ŀ�꺯��ֵ XnewScore(i)
                if fworst > XnewScore(i)
                    %���½� Xnew �滻 Elite �����е�����ľ��߱�����
                    Task(k).Elite(i).X(worstId,:) = Xnew;
                    %���½� Xtemp �滻 Elite �����е���������������Ϣ��
                    Task(k).Elite(i).HM(worstId,:) = Xtemp;
                    %���µ��������� XnewScore �滻 Elite �����е����������������
                    Task(k).Elite(i).Fit(worstId,:) = XnewScore;
                    for s = 1:bestNum%�������� s ��1��bestNum������ѭ����
                        %����½� Xnew �еĵ� i ��Ŀ�꺯��ֵ XnewScore(i) С�� Elite �����еĵ� i ������������ֵ
                        if XnewScore(i) < Task(k).Elite(i).fbest(s)
                            %�����ֵ fbest(s) ����Ϊ�µ�Ŀ�꺯��ֵ XnewScore(i)��
                            Task(k).Elite(i).fbest(s) = XnewScore(i);
                            %����Ѿ��߱��� Xbest(s, :) ����Ϊ�½� Xnew��
                            Task(k).Elite(i).Xbest(s,:) = Xnew;
                            %% �����ǰά�� dim С�� epi_dim��
                            if dim < epi_dim
                                %�����ݵ���չ�Ӽ����ж��׼������������洢�� Score �У����� Xnew ��ǰ dim+1 ��Ԫ�ر����ݸ��ú�����
                                [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State);
                                %���������� Score ��һ���� [0, 1] ��Χ�ڣ�ͨ��������� maxFit��
                                Score = Score ./ maxFit;
                                %�����չ���ά�� dim+1 ���� Epi_Dim��
                                if dim+1 == Epi_Dim
                                    %���� Epi_Dim_FEs ��������ֵ��
                                    Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                end
                                NC = NC + 1;%���� NC ��������ֵ��
                                for si = 1:FitNum
                                    %�������� sj ��1��EliteSize������ѭ����
                                    for sj = 1:EliteSize
                                        %�����չ��� Elite �����еĵ� si ��������ĵ� sj �����Ŀ�꺯��ֵ������չ���Ŀ�꺯��ֵ
                                        if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                            %����չ��ľ��߱��� Xnew �滻Ϊ Elite �����еĵ� si ��������ĵ� sj ���⣬�������������
                                            Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                            %����չ��ľ��߱��� Xnew �滻Ϊ Elite �����еĵ� si ��������ĵ� sj ��������������Ϣ��
                                            Task(dim).Elite(si).HM(sj,:) = Xnew;
                                            %����չ����������� Score �滻Ϊ Elite �����еĵ� si ��������ĵ� sj ���������������
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

        %if ismember(CX(ci), Xnew) %��� CX(ci) �Ƿ������ Xnew ��  CX���²�ģ��

        %   cflag = cflag + 1;
        %end
        %end
        %if cflag == fdim
        %�� Xnew ��ֵ������һ���ض�������Ԫ�أ��漰�˶������ݽṹ������
        %   Task(fdim-1).Elite(1).X(1,:) = Xnew;
        %�� XnewScore ��ֵ������һ���ض�������Ԫ�ء�
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