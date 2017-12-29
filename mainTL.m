%% 
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% gt_a=double(im_gt)+1;
% Xa=reshape(im,[],size(im,3));
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% gt_b=double(im_gt)+1;
% Xb=reshape(im,[],size(im,3));
% ims=Xa;imt=Xb;
% ims_gt=gt_a;imt_gt=gt_b;
% id_cluster = kmeans(imt,5,'MaxIter',10000,'OnlinePhase','on','Replicates',5,'Options',statset('UseParallel',1),'Display','final');%聚类
%% 
% midx=[52 80 23];
% load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
% Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
% [~, Xa] = pca(Xa,'Centered',true,'NumComponents',3);%,'NumComponents',5
% % Xa=Xa(:,midx);
% Xa=normcols(Xa);
% load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
% Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
% [~, Xb] = pca(Xb,'Centered',true,'NumComponents',3);%,'NumComponents',5
% % Xb=Xb(:,midx);
% Xb=normcols(Xb);
% for kk=1:size(Xb,2)
%     Xa(:,kk)=histeq(Xa(:,kk),100);
%     Xb(:,kk)=histeq(Xb(:,kk),100);
% end
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
% % histogram
% num_bins=100;
% [KL_stmat1,KL_tsmat1]=CalculateKL(Xa,gt_a,Xb,gt_b,num_bins);
% [Ew_hist1,Ew_hist2,eval1]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins);
% % KLtree
% num_clusters=100;
% [KL_stmat2,KL_tsmat2]=CalculateKL_Tree(Xa,gt_a,Xb,gt_b,num_clusters);
% [Ew_tree1,Ew_tree2,eval2]=CalculateEw_Tree(Xa,gt_a,Xb,gt_b,num_clusters);

%% 
n_e=4;
load(['E:\TransfLearning\PUC\CNN\PU\data',num2str(n_e),'\PcFeature.mat'], 'feat_conv','feat_fc')
% feat_fc=normcols(feat_fc);
% feat_conv=normcols(feat_conv);
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
Xa=normcols(Xa);
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Pc_same=Pc_same(:,224:end,:);
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt(:,224:end);%
Xb=normcols(Xb);
[ims,imt]=pavia_adjust(Xa,Xb,100,0.001);%灰度拉伸
num_bins=100;num_cluster=7;
inda=find(gt_a);indb=find(gt_b);%==1|gt_a==2|gt_a==6|gt_a==7 ==1|gt_b==2|gt_b==6|gt_b==7
ims=ims(inda,:);imt=imt(indb,:);
ims_gt=gt_a(inda);imt_gt=gt_b(indb);
% load(['PC_id_cluster',num2str(num_cluster)],'id_cluster');%PaviaC全图
id_cluster = kmeans([feat_fc],num_cluster,'MaxIter',10000,'OnlinePhase','on'...
    ,'Replicates',10,'Options',statset('UseParallel',1),'Distance','correlation');%聚类
[Ew1,Ew2,Ew3,eval]=CalculateEw(ims,ims_gt,imt,id_cluster,num_bins);% 利用JS散度计算边权重
pairs=iter_match(Ew2);%迭代匹配
matched_pair=cell2mat(pairs);

[~,I] = sort(matched_pair(1,:),'ascend');%% 可以用sortrows
real_matched=matched_pair(2,I);
L=4;% 对匹配的类簇上进行二次聚类匹配
for k=1:max(ims_gt(:))
    indb=find(id_cluster==real_matched(k));
    tempa=ims(ims_gt==k,:);
    tempb=imt(indb,:);
    clustera = kmeans(tempa,L,'MaxIter',10000,'OnlinePhase','on','Replicates',10,'Options',statset('UseParallel',1));%聚类,'Display','final'
    clusterb = kmeans(tempb,L,'MaxIter',10000,'OnlinePhase','on','Replicates',10,'Options',statset('UseParallel',1));%聚类,'Display','final'
    [~,~,EwL,~]=CalculateEw(tempa,clustera,tempb,clusterb,100);
    [~,m]=find(EwL==min(EwL(:)));
    all_ind{k}=indb(clusterb==m);
    acc(k)=mean(imt_gt(all_ind{k})==k);%检查二次匹配聚类效果
end
acc
% certain_pair=pairs{1,1};%certain pairs only
certain_pair=matched_pair;% all pairs
train_ind=[];y_train=[];val_ind=[];y_val=[];test_ind=[];y_test=[];
for k=1:size(certain_pair,2)
    ind_all=all_ind{certain_pair(1,k)};
    ind_all=ind_all(randperm(length(ind_all)));
    half_ind=round(0.5*length(ind_all));
    train_ind=[train_ind;ind_all(1:half_ind)];
    y_train=[y_train;certain_pair(1,k)*ones(half_ind,1)];
    val_ind=[val_ind;ind_all(half_ind+1:end)];
    y_val=[y_val;certain_pair(1,k)*ones(length(ind_all)-half_ind,1)];
end
load E:\TransfLearning\PUC\CNN\PcIMAGES.mat x_test y_test
y_train0=double(y_test(train_ind));
y_val0=double(y_test(val_ind));
% mean(y_train==y_train0)
% mean(y_val==y_val0)
x_trainT=x_test(:,:,:,train_ind);
x_valT=x_test(:,:,:,val_ind);
y_trainT=y_train;y_valT=y_val;
save(['E:\TransfLearning\PUC\CNN\PU\data',num2str(n_e),'\finetuning\DataForKerasT.mat'],...
    'x_trainT', 'y_trainT', 'x_valT', 'y_valT','-v7.3'); %matching pairs Target domain samples
% 取出PU训练样本，将PC中可靠类的样本替换
clearvars -except x_trainT x_valT y_trainT y_valT certain_pair Ew2
load(['E:\TransfLearning\PUC\CNN\PU\data',num2str(n_e),'\DataForKeras.mat'])
rm_ind_train=[];rm_ind_val=[];
for k_c=1:size(certain_pair,2)
    rm_ind_train=[rm_ind_train;find(y_train==certain_pair(1,k_c))];
    rm_ind_val=[rm_ind_val;find(y_val==certain_pair(1,k_c))];
end
x_train(:,:,:,rm_ind_train)=[];
y_train(rm_ind_train)=[];
x_val(:,:,:,rm_ind_val)=[];
y_val(rm_ind_val)=[];
x_train=cat(4,x_train,x_trainT);
y_train=[y_train;y_trainT];
x_val=cat(4,x_val,x_valT);
y_val=[y_val;y_valT];
class_weight=ones(7,1);
class_weight(certain_pair(1,:))=(length(y_train)-length(y_trainT))/length(y_trainT);
save(['E:\TransfLearning\PUC\CNN\PU\data',num2str(n_e),'\finetuning\DataForKerasTR.mat'],...
    'x_train','y_train','x_val','y_val','-v7.3') %将可靠类的PC样本替换PU样本% Target samples Replace some classes

    
