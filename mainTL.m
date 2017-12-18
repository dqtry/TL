%% 
load('E:\TransfLearning\area7_source\im1.mat', 'im')
load('E:\TransfLearning\area7_source\im1.mat', 'im_gt')
ims_gt=im_gt+1;
ims=reshape(im,[],size(im,3));
load('E:\TransfLearning\area1_target\im1.mat', 'im')
load('E:\TransfLearning\area1_target\im1.mat', 'im_gt')
imt_gt=im_gt+1;
imt=reshape(im,[],size(im,3));

[KL_stmat,KL_tsmat]=CalculateKL_Tree(ims,ims_gt,imt,imt_gt,100);
[KL_stmat2,KL_tsmat2]=CalculateKL(ims,ims_gt,imt,imt_gt,50);
[Ew,Ew2,eval]=CalculateEw(ims,ims_gt,imt,imt_gt,50);

%% 
midx=[52 80 23];
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
ims=reshape(Pu_same,[],size(Pu_same,3));ims_gt=Ugt;
% [~, ims] = pca(ims,'Centered',true,'NumComponents',5);%,'NumComponents',5
% ims=ims(:,midx);
ims=normcols(ims);

load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
imt=reshape(Pc_same(:,224:end,:),[],size(Pc_same,3));imt_gt=Cgt(:,224:end);
% [~, imt] = pca(imt,'Centered',true,'NumComponents',5);%,'NumComponents',5
% imt=imt(:,midx);
imt=normcols(imt);
[KL_stmat,KL_tsmat]=CalculateKL_Tree(ims,ims_gt,imt,imt_gt,100);
Ew(k1,k2)=sqrt(0.5*KL_au+0.5*KL_bu);
Ew2(k1,k2)=sqrt(min(KL_au,KL_bu));