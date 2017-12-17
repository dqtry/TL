%% 
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% gt_a=double(im_gt)+1;
% Xa=reshape(im,[],size(im,3));
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% gt_b=double(im_gt)+1;
% Xb=reshape(im,[],size(im,3));
%% 
midx=[52 80 23];
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
[~, Xa] = pca(Xa,'Centered',true,'NumComponents',3);%,'NumComponents',5
% Xa=Xa(:,midx);
Xa=normcols(Xa);
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
[~, Xb] = pca(Xb,'Centered',true,'NumComponents',3);%,'NumComponents',5
% Xb=Xb(:,midx);
Xb=normcols(Xb);
for kk=1:size(Xb,2)
    Xa(:,kk)=histeq(Xa(:,kk),100);
    Xb(:,kk)=histeq(Xb(:,kk),100);
end
%% 域内分成两部分进行验证JS距离有效
% rng(0);% 设置种子点，可复现
% imt=Xa;imt_gt=gt_a;
% Xa=[];Xb=[];gt_a=[];gt_b=[];
% for k1=1:max(imt_gt(:))
%     temp = find(imt_gt==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=imt(temp,:);
%     Xa=[Xa;temp_feat(1:round(length(temp)/2),:)];
%     gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
%     Xb=[Xb;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
%% 
% histogram
num_bins=100;
[KL_stmat1,KL_tsmat1]=CalculateKL(Xa,gt_a,Xb,gt_b,num_bins);
[Ew_hist1,Ew_hist2,eval1]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins);
% KLtree
num_clusters=100;
[KL_stmat2,KL_tsmat2]=CalculateKL_Tree(Xa,gt_a,Xb,gt_b,num_clusters);
[Ew_tree1,Ew_tree2,eval2]=CalculateEw_Tree(Xa,gt_a,Xb,gt_b,num_clusters);