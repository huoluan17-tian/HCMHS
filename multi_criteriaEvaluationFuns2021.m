
%function [K2Score,Gtest,JE_Score,ME,JS] = multi_criteriaEvaluationFuns2021(snp_com,state)
%function [K2Score,LR,ME,P_value,JS] = multi_criteriaEvaluationFuns2021(snp_com,state)
%function [K2Score,Gtest,JE_Score,ME,JS] = multi_criteriaEvaluationFuns2021(snp_com,state)%已记录
%function [JE_Score,Gtest,JS,ME,K2Score] = multi_criteriaEvaluationFuns2021(snp_com,state)
function [K2Score,LR,ME,JE_Score,JS] = multi_criteriaEvaluationFuns2021(snp_com,state)

%%找出状态向量 state 中不同的值，并存储在 ua 中：0或1
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0%旨在确保状态向量的值从0开始计数，而不是从一个非零值开始。
    state = state -min(ua);
end
%获取了矩阵 snp_com 的行数和列数，并将它们分别存储在 xrow 和 xcol 中
[xrow,xcol] = size(snp_com);
[Data,idx,cid]=unique(snp_com,'rows');%找出矩阵 snp_com 中唯一的基因型，并将它们存储在 Data 中，它们在snp_com中对应的位置放在idx中 cid表示snp_com中的每个与data基因型在data中的下标
[lrow,~]=size(Data);

Hcase = zeros(3,xcol);%统计两列出现的0，1，2数目 
Hcontrol = zeros(3,xcol);%统计两列出现的0，1，2数目 
sample=zeros(lrow,1);
disease=sample;%统计的疾病组9个基因型出现的数目
control=sample;%统计的对照组9个基因型出现的数目

for i=1:xrow   %% 统计每个基因型组合出现的次数 对上面进行填值
   if state(i) == 1%疾病组
       disease(cid(i)) = disease(cid(i)) + 1;
         for j = 1:xcol
              if snp_com(i,j) == 0
                  Hcase(1,j) = Hcase(1,j) + 1;
              elseif snp_com(i,j) == 1
                  Hcase(2,j) = Hcase(2,j) + 1; 
              else%
                  Hcase(3,j) = Hcase(3,j) + 1;
              end
         end
   else%对照组
       control(cid(i)) = control(cid(i)) + 1;
        for j = 1:xcol
              if snp_com(i,j) == 0
                  Hcontrol(1,j) = Hcontrol(1,j) + 1;
              elseif snp_com(i,j) == 1
                  Hcontrol(2,j) = Hcontrol(2,j) + 1; 
              else
                  Hcontrol(3,j) = Hcontrol(3,j) + 1;
              end
         end
   end
  
end
ebsol = 0.1;
A = sum((Hcase - Hcontrol).^2);%计算了 (Hcase - Hcontrol)^2 的元素和

%sample里面存储的9个基因型的数目 不分致病对照。
sample = disease + control;
%计算了两个向量disease和control之间的标准差。
 SDC = sqrt(sum(A))/sum(abs(disease-control));
 %% G-test 用到了表型
 %将disease和control两个向量按列排列在一起，形成一个2行lrow列的矩阵。然后，计算了矩阵的行和列的和，并将结果存储在F(3,1:lrow)和F(:,lrow+1)中。
     F = [disease';control'];
     F(3,1:lrow)=sum(F(1:2,1:lrow),1);
     F(:,lrow+1)=sum(F(:,1:lrow),2);
     %初始化了变量G、chi和Degree。G和chi被初始化为0，Degree被计算为(2-1)*(lrow-1)的值。
     G=0;
     chi = 0;
    % LR = 0;
     Degree=(2-1)*(lrow-1);
     %嵌套的循环结构，用于计算G-test的值。在每次迭代中，它计算了观察值（O）和期望值（E），并根据这些值更新了G和chi的值。根据F(3,j)>10的条件，如果满足条件，G和chi的值将会被更新。同时，Degree的值在每次迭代中也可能会减小。
    for i=1:2
        for j=1:lrow
            O = F(i,j)/xrow;
            E = ( ( F(i,lrow+1) * F(3,j) )/xrow )/xrow;

              if F(3,j)>10
                   if O>0
                      G = G+(xrow * O) * log(O/E);
                     % LR = LR + O * log(O/E);
                      chi = chi + (abs(O - E))^2 / E;
                   end

              elseif Degree>1 
                 Degree = Degree-0.5; 
              end
        end    
    end
    %最终计算chi和G的值。首先，chi和G的值分别乘以2，然后用于计算LR和Gtest的值。LR的计算公式是1 / G，而Gtest的计算使用了chi-square分布的累积分布函数（chi2cdf）来得到。
chi = 2 * chi;
% chiT = 1 / chi;
G=2*G;
LR = 1 / G;
Gtest = 1 - chi2cdf(G,Degree);
%1/LR;
%GtestP_value=1-chi2cdf(G,Degree);

% [disease';control';sample']
%% K2 score  用到了表型
    K2score=0;

    for i=1:lrow
        if sample(i)>0
            y=My_factorial(sample(i)+1);
            r=My_factorial(disease(i))+My_factorial(control(i));
            K2score=K2score+(r-y);
        end
    end
   K2Score =abs(K2score);



%% GINI score 用到了表型
 

        sCase = sum(disease);
        Pcase = disease./sCase;%
        Pcontrol = control./(xrow - sCase);%
        P = sample / xrow;

         Gini_Score = sCase/xrow * (1 - sum(Pcase.^2)) / ((1-sCase/xrow) * (1-sum(Pcontrol.^2)));

        
  
 %% 互熵 mutual entropy 用到了表型
     Psample = sample / xrow;
     PCC = [disease; control]/xrow;
    
     PCC = PCC(PCC>0);
     s1 = sCase/xrow; s2 = 1 - s1;
     MeY = - (s1 .* log2(s1) + s2 .* log2(s2));
   
     MeX = - sum(Psample.*log2(Psample));
     MeXY = - sum ( PCC .* log2(PCC));
     ME = 1 / ( MeY + MeX - MeXY);


%% Joint Entropy of disease genotype combinantion  用到了表型
        Psample = sample / xrow;
        Pdisease = disease / sum(disease);
        Pcontrol = control / sum(control);
        JE = 0;
        JC = 0;
        JCC = 0;
        for i = 1:lrow
           if Pdisease(i)>0
              JE = JE + Pdisease(i).*log2(Pdisease(i));
           end
            if Pcontrol(i)>0
              JC = JC + Pcontrol(i).*log2(Pcontrol(i));
            end
%             JCC = JCC + Pdisease(i).*log2(Pdisease(i)) / (Pcontrol(i).*log2(Pcontrol(i)));
            
        end


            
         %    JE_Score = SDC / (JE + JC)^2 ;% (JE + JC)^2  is quadratic sum of the joint entropy of SNP combintion
          JE_Score = (SDC / JC^2) ;
      %   JE_Score = SDC;
   
 
          
          
          %% Jensen-Shannon
 JS1 = 0;
 JS2 = 0;
 JF = 0; %Jeffrey
     for i = 1:lrow
           if Pdisease(i)>0
              JS1 = JS1 + Pdisease(i).*log2(2*Pdisease(i) / (Pdisease(i) + Pcontrol(i)) );
           end
            if Pcontrol(i)>0
              JS2 = JS2 + Pcontrol(i).*log2(2*Pcontrol(i) / (Pdisease(i) + Pcontrol(i)) );
            end
           
%             if Pdisease(i)+Pcontrol(i) > 0
%                 JF = JF + ((Pdisease(i)-Pcontrol(i))^2) / Pcontrol(i);
%             end
     end
      %JS =1- 0.5*(JS1 + JS2) ;
      JS = 1 /(JS1 + JS2);
      
%% Wasserstein Distance
%    WS = 1 / sum((Pdisease-Pcontrol).^2);
 

%    Gtest= K2Score;
%    Gini_Score = K2Score;
%    JE_Score = K2Score;


%% f is function used to calculate log form factorial
    function f=My_factorial(e)
     
        f=0;
        if e>0
            for o=1:e
                f=f+log(o);
            end
        end
    end



end


