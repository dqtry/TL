function [Ew1,Ew2,Ew3,eval]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins)
% Ew1是未完全理解JS散度的不准确的结果（a和b直接合并未平均），存在大于1的值
% Ew2和Ew3根据定义（KL散度）以及公式推导（根据分布）得到的
% eval是自定义评价指标，不是很可靠，可忽略
%% 直方图统计分布，用于计算Ew(JS散度的平方根)
% 将源域每类分成两部分进行验证哪种直方图适合于计算Ew
% A Novel Graph-Matching-Based Approach for Domain Adaptation in Classification of Remote Sensing Image Pair
% Ew(JS散度的平方根)是一个对称的距离，用于表示带权无向图中的权重，由KL散度得到
% KL散度是相对熵，即交叉熵减去理论分布的信息熵，（非负、非对称）
% Xa Xb 行表示样本，列表示特征 0～1之间的浮点型
% gt_a gt_b 表示对于的类簇号
% num_bins 直方图统计时的间隔数

%% ==============同一副图分成两部分进行测试==================
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
% % figure,hold on
% % for kkk=1:5
% % rng();
% Xa=[];Xb=[];gt_a=[];gt_b=[];
% for k1=1:length(unique(ims_gt))
%     temp = find(ims_gt==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=ims(temp,:);
%     Xa=[Xa;temp_feat(1:round(length(temp)/2),:)];
%     gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
%     Xb=[Xb;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
% num_bins=50;
% 直方图统计最好进行平滑操作
%% ==============两张图进行测试================
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% Xa=reshape(im,[],size(im,3));gt_a=im_gt+1;
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% Xb=reshape(im,[],size(im,3));gt_b=im_gt+1;
% num_bins=50;
%% 1.一维直方图，分波段，每波段用100个bin进行统计，平滑之后进行一维连接
% num_bins=[11 21 31 41 51 61 71 81 91 101 111 121 131 141 151];
% % for kt=1:length(num_bins)
edges=linspace(min(Xa(:)),max(Xa(:)),num_bins+1);%理论上应该归一化0~1
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    for k2=1:max(gt_b(:))
        tempa=find(gt_a==k1);
        tempb=find(gt_b==k2);
        feat_u=[Xa(tempa,:);Xb(tempb,:)];
        Nas=[];Nbs=[];Nus=[];
        for k3=1:size(Xa,2)
            [Na,~]=histcounts(Xa(tempa,k3),edges);
            Na=conv(Na,se_mask,'same');% 先进行平滑然后连接
            Nas=[Nas,Na+eps];
            [Nb,~]=histcounts(Xb(tempb,k3),edges);
            Nb=conv(Nb,se_mask,'same');% 先进行平滑然后连接
            Nbs=[Nbs,Nb+eps];
            [Nu,~]=histcounts(feat_u(:,k3),edges);
            Nu=conv(Nu,se_mask,'same');% 先进行平滑然后连接
            Nus=[Nus,Nu+eps];
        end
        Pa=Nas./sum(Nas);Pb=Nbs./sum(Nbs);Pu=Nus./sum(Nus);%概率
        KL_au(k1,k2)=sum(Pa.*log(Pa./Pu));
        KL_bu(k1,k2)=sum(Pb.*log(Pb./Pu));
        Ew1(k1,k2)=sqrt(0.5*KL_au(k1,k2)+0.5*KL_bu(k1,k2));
%         Ew(k1,k2)=sqrt(max(KL_au(k1,k2),KL_bu(k1,k2)));
%         Ew2(k1,k2)=sqrt(min(KL_au(k1,k2),KL_bu(k1,k2)));
        Pu2=0.5*Pa+0.5*Pb;
        KL_au2(k1,k2)=sum(Pa.*log(Pa./Pu2));
        KL_bu2(k1,k2)=sum(Pb.*log(Pb./Pu2));
        Ew2(k1,k2)=sqrt(0.5*KL_au2(k1,k2)+0.5*KL_bu2(k1,k2));
        Ew3(k1,k2)=sqrt(0.5*(sum(Pa.*log(Pa))+sum(Pb.*log(Pb)))-sum(0.5*(Pa+Pb).*log(0.5*(Pa+Pb))));
    end
end
Ewt=Ew1+100*eye(size(Ew1));%为了取得除了非对角线元素的最小值
% 非对角线元素的最小值与对角线元素的最大值之差，即类间距离的最小值与类内距离的最大值之差要大
eval(1)=min(Ewt(:))-max(diag(Ew1));
Ewt2=Ew2+100*eye(size(Ew2));
eval(2)=min(Ewt2(:))-max(diag(Ew2));
% % Ew2=Ew+100*eye(size(Ew));
% % eval(kt)=min(Ew2(:))-max(diag(Ew));
% % end
% % plot(num_bins,eval)
% % end
% % hold off