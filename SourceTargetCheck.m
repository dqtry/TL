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
ims=reshape(Pu_same,[],size(Pu_same,3));ims_gt=Ugt;
% [~, ims] = pca(ims,'Centered',true,'NumComponents',5);%,'NumComponents',5
ims=ims(:,midx);
ims=normcols(ims);

load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
imt=reshape(Pc_same(:,224:end,:),[],size(Pc_same,3));imt_gt=Cgt(:,224:end);
% [~, imt] = pca(imt,'Centered',true,'NumComponents',5);%,'NumComponents',5
imt=imt(:,midx);
imt=normcols(imt);
%% 计算Ew
% num_bins=50;
% [Ew,Ew2,eval]=CalculateEw(ims,ims_gt,imt,imt_gt,num_bins);%% 计算Ew距离
%% 域内分成两部分进行验证JS距离有效
% rng(0);% 设置种子点，可复现
% [~,~,labels]=find(imt_gt);% 避免有标记0
% feat_a=[];feat_b=[];gt_a=[];gt_b=[];
% for k1=1:length(unique(labels))
%     temp = find(imt_gt==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=imt(temp,:);
%     feat_a=[feat_a;temp_feat(1:round(length(temp)/2),:)];
%     gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
%     feat_b=[feat_b;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
% [Ew,Ew2,eval]=CalculateEw(feat_a,gt_a,feat_b,gt_b,num_bins);%% 计算Ew距离
%% ---------PU和PC显示效果检测-----------
% % [~,C,~,~,midx]=kmedoids(ims',3);% midx对应RGB可能需要调整
% midx=[52 80 23];
% imchoset=reshape(imt,size(Pc_same));
% imshow(imchoset(:,:,midx),[])
% imchoses=reshape(ims,size(Pu_same));
% figure,imshow(imchoses(:,:,midx),[])
%% 检查源域和目标域分布情况
figure(1)
xlim([0 1]);ylim([0 1]);
figure(2)
xlim([0 1]);ylim([0 1]);
for k=1:length(unique(ims_gt))-1% gt中含有0
    temps = find(ims_gt==k);
    figure(1),hold on
    scatter(ims(temps,1),ims(temps,2));
    tempt = find(imt_gt==k);
    figure(2),hold on
    scatter(imt(tempt,1),imt(tempt,2));
end




