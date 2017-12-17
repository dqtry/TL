%% 检查源域和目标域的分布变化情况
%% ISPRS
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
% % % for nb=1:size(ims,2)
% % %     figure,histogram(ims(:,nb));
% % % end
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% imt=reshape(im,[],size(im,3));imt_gt=im_gt+1;
% for k=1:length(unique(ims_gt))
%     temps = find(ims_gt==k);
%     tempt = find(imt_gt==k);
%     figure(1),hold on
%     scatter(ims(temps,1),ims(temps,3));
%     figure(2),hold on
%     scatter(imt(tempt,1),imt(tempt,3));
% end
%% PU和PC
midx=[52 80 23];
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
Xa=normcols(Xa);
Xb=normcols(Xb);
for kk=1:size(Xb,2)
    Xa(:,kk)=histeq(Xa(:,kk),256);
    Xb(:,kk)=histeq(Xb(:,kk),256);
end
[~, Xb] = pca(Xb,'Centered',true,'NumComponents',3);%,'NumComponents',5
% Xb=Xb(:,midx);
Xb=normcols(Xb);
[~, Xa] = pca(Xa,'Centered',true,'NumComponents',3);%,'NumComponents',5
% Xa=Xa(:,midx);
Xa=normcols(Xa);

edges=linspace(0,1,101);
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    figure;hold on
    tempa=find(gt_a==k1);
    [Na,~]=histcounts(Xa(tempa,2),edges);
    Na=conv(Na,se_mask,'same');% 先进行平滑然后连接
    tempb=find(gt_b==k1);
    [Nb,~]=histcounts(Xb(tempb,2),edges);
    Nb=conv(Nb,se_mask,'same');% 先进行平滑然后连接
    plot(edges(1:end-1),Na,edges(1:end-1),Nb);%,edges(1:end-1),Nb
end
%% 计算Ew
num_bins=50;
[Ew,Ew2,eval]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins);%% 计算Ew距离
%% 域内分成两部分进行验证JS距离有效
% rng(0);% 设置种子点，可复现
% feat_a=[];feat_b=[];gt_a=[];gt_b=[];
% for k1=1:max(gt_b(:))
%     temp = find(gt_b==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=Xb(temp,:);
%     feat_a=[feat_a;temp_feat(1:round(length(temp)/2),:)];
%     gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
%     feat_b=[feat_b;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
% [Ew,Ew2,eval]=CalculateEw(feat_a,gt_a,feat_b,gt_b,num_bins);%% 计算Ew距离
%% ---------PU和PC显示效果检测-----------
% % [~,C,~,~,midx]=kmedoids(ims',3);% midx对应RGB可能需要调整
% midx=[52 80 23];
% imchoset=reshape(Xb,size(Pc_same));
% imshow(imchoset(:,:,midx),[])
% imchoses=reshape(Xa,size(Pu_same));
% figure,imshow(imchoses(:,:,midx),[])
%% 检查源域和目标域分布情况
% figure(1)
% xlim([0 1]);ylim([0 1]);
% figure(2)
% xlim([0 1]);ylim([0 1]);
% for k=1:max(gt_a(:))% gt中含有0
%     temps = find(gt_a==k);
%     figure(1),hold on
%     scatter(Xa(temps,1),Xa(temps,2));
%     tempt = find(gt_b==k);
%     figure(2),hold on
%     scatter(Xb(tempt,1),Xb(tempt,2));
% end




