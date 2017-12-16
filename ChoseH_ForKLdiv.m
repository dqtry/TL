%% 寻找最佳直方图表示，用于计算KL散度
% 将源域每类分成两部分进行验证哪种直方图适合于计算KL散度
% A Novel Graph-Matching-Based Approach for Domain Adaptation in Classi?cation of Remote Sensing Image Pair
% KL散度是相对熵，即交叉熵减去理论分布的信息熵，（非负、非对称）
load E:\TransfLearning\area7_source\im1.mat im im_gt
ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;

feat_a=[];feat_b=[];gt_a=[];gt_b=[];
for k1=1:length(unique(ims_gt))
    temp = find(ims_gt==k1);
    temp = temp(randperm(length(temp)));
    temp_feat=ims(temp,:);
    feat_a=[feat_a;temp_feat(1:round(length(temp)/2),:)];
    gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
    feat_b=[feat_b;temp_feat(round(length(temp)/2)+1:end,:)];
    gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
end
% 直方图统计最好进行平滑操作
%% 1.一维直方图，分波段，每波段用100个bin进行统计，平滑之后进行一维连接
% % feat_b=feat_a;%% 用相同的数据进行测试
% % gt_b=gt_a;
edges1=linspace(0,1,11);
se_mask=fspecial('gaussian',[1,3]);
for k1=1:length(unique(gt_a))
    for k2=1:length(unique(gt_b))
        Nas=[];Nbs=[];
        for k3=1:size(ims,2)
            tempa=find(gt_a==k1);
            [Na,~]=histcounts(feat_a(tempa,k3),edges1);
            Na=conv(Na+eps,se_mask,'same');% 先进行平滑然后连接
            Nas=[Nas,Na];
            tempb=find(gt_b==k2);
            [Nb,~]=histcounts(feat_b(tempb,k3),edges1);
            Nb=conv(Nb+eps,se_mask,'same');% 先进行平滑然后连接
            Nbs=[Nbs,Nb];
        end
        Pa=Nas./sum(Nas);
        Pb=Nbs./sum(Nbs);
        KL_st1(k1,k2)=sum(Pa.*log(Pa./Pb));
        KL_ts1(k2,k1)=sum(Pb.*log(Pb./Pa));
    end
end
KL_stmat=cell2mat(KL_st1);
KL_tsmat=cell2mat(KL_ts1);
F1=1/(abs(sum(KL_stmat(:))-sum(KL_tsmat(:)))/numel(KL_stmat)+1);%%%%%%%%%%%+1
%% 2.二维直方图（B1+B2 B1+B3 B2+B3）统计之后得到三维数据，然后归一化
% se_mask2d=fspecial('gaussian',[3 3]);
% edges2=linspace(0,1,11);
% for k1=1:length(unique(gt_a))
%     for k2=1:length(unique(gt_b))
%         tempa=find(gt_a==k1);
%         [Na1,~,~] = histcounts2(feat_a(tempa,1),feat_a(tempa,2),edges2,edges2);
%         Na1=conv2(Na1+eps,se_mask2d,'same');% 先进行平滑然后组合
%         [Na2,~,~] = histcounts2(feat_a(tempa,1),feat_a(tempa,3),edges2,edges2);
%         Na2=conv2(Na2+eps,se_mask2d,'same');% 先进行平滑然后组合
%         [Na3,~,~] = histcounts2(feat_a(tempa,2),feat_a(tempa,3),edges2,edges2);
%         Na3=conv2(Na3+eps,se_mask2d,'same');% 先进行平滑然后组合
%         Nas=cat(3,Na1,Na2,Na3);
%         Pa=Nas./sum(Nas(:));
%         
%         tempb=find(gt_b==k2);
%         [Nb1,~,~] = histcounts2(feat_b(tempb,1),feat_b(tempb,2),edges2,edges2);
%         Nb1=conv2(Nb1+eps,se_mask2d,'same');% 先进行平滑然后组合
%         [Nb2,~,~] = histcounts2(feat_b(tempb,1),feat_b(tempb,3),edges2,edges2);
%         Nb2=conv2(Nb2+eps,se_mask2d,'same');% 先进行平滑然后组合
%         [Nb3,~,~] = histcounts2(feat_b(tempb,2),feat_b(tempb,3),edges2,edges2);
%         Nb3=conv2(Nb3+eps,se_mask2d,'same');% 先进行平滑然后组合
%         Nbs=cat(3,Nb1,Nb2,Nb3);
%         Pb=Nbs./sum(Nbs(:));
%         
%         t1=Pa.*log(Pa./Pb);
%         t2=(Pb.*log(Pb./Pa));
%         KL_st2(k1,k2)=sum(sum(sum(Pa.*log(Pa./Pb))));
%         KL_ts2(k2,k1)=sum(sum(sum(Pb.*log(Pb./Pa))));
%     end
% end      
%% 3.对图像进行PCA
% PCA,edges进行相应的修改
% [~,ims,~] = pca(ims,'Centered',false);
% ims=normcols(ims);
% % ------------------------
% % 一维直方图，分波段，每波段用100个bin进行统计，平滑之后进行一维连接
% % % feat_b=feat_a;%% 用相同的数据进行测试
% % % gt_b=gt_a;
% % ------------------------
% edges3{1,:}=linspace(0,1,101);
% edges3{2,:}=linspace(0,1,51);
% edges3{3,:}=linspace(0,1,41);
% se_mask=fspecial('gaussian',[1,3]);
% for k1=1:length(unique(gt_a))
%     for k2=1:length(unique(gt_b))
%         Nas=[];Nbs=[];
%         for k3=1:size(ims,2)
%             tempa=find(gt_a==k1);
%             [Na,~]=histcounts(feat_a(tempa,k3),edges3{k3,:});
%             Na=conv(Na+eps,se_mask,'same');% 先进行平滑然后连接
%             Nas=[Nas,Na];
%             tempb=find(gt_b==k2);
%             [Nb,~]=histcounts(feat_b(tempb,k3),edges3{k3,:});
%             Nb=conv(Nb+eps,se_mask,'same');% 先进行平滑然后连接
%             Nbs=[Nbs,Nb];
%         end
%         Pa=Nas./sum(Nas);
%         Pb=Nbs./sum(Nbs);
%         KL_st3(k1,k2)=sum(Pa.*log(Pa./Pb));
%         KL_ts3(k2,k1)=sum(Pb.*log(Pb./Pa));
%     end
% end
%% 4.高斯混合模型进行拟合(fitgmdist)，然后计算距离


        
        
        
  