%% 利用聚类结果构建树，而后计算KL散度
% <treeKL: A distance between high dimension empirical distributions>
% clear
load('E:\TransfLearning\area7_source\im1.mat', 'im')
load('E:\TransfLearning\area7_source\im1.mat', 'im_gt')
ims_gt=im_gt+1;
ims=reshape(im,[],size(im,3));
load('E:\TransfLearning\area1_target\im1.mat', 'im')
load('E:\TransfLearning\area1_target\im1.mat', 'im_gt')
imt_gt=im_gt+1;
imt=reshape(im,[],size(im,3));
epsilon=0;se_mask=fspecial('gaussian',[1,3]);
for kkk=1:10% 检测聚类的随机性
%% 构建源域树及类统计
for k_s=1:max(ims_gt(:))
    ind= ims_gt==k_s;% tree
    Xsk=normcols(ims(ind,:));
    IDX=kmeans(Xsk,50,'MaxIter',10000);% kmeans
    model_s{k_s}=fitctree(Xsk,IDX);
    pred_s=predict(model_s{k_s},Xsk);
%     mean(IDX==pred_s);
    [Ns,edges{k_s}] = histcounts(pred_s);
    Ns = conv(Ns,se_mask);
    ps_s{k_s}=Ns./sum(Ns);
end
%% 构建目标域树及簇统计
for k_t=1:max(imt_gt(:))
    ind= imt_gt==k_t;% tree
    Xtk=normcols(imt(ind,:));
    IDX=kmeans(Xtk,50,'MaxIter',10000);% kmeans
    model_t{k_t}=fitctree(Xtk,IDX);
    pred_t=predict(model_t{k_t},Xtk);
    [Nt,edget{k_t}] = histcounts(pred_t);
    Nt = conv(Nt,se_mask);
    pt_t{k_t}=Nt./sum(Nt);
end
%% 源域树+目标域簇计算KL
for k_s=1:max(ims_gt(:))
    for k_t=1:max(imt_gt(:))
        ind=imt_gt==k_t;
        Xtk=normcols(imt(ind,:));
        pred_t=predict(model_s{k_s},Xtk);
        Nt = histcounts(pred_t,edges{k_s});
        Nt = conv(Nt,se_mask);
        ps_t{k_s,k_t}=Nt./sum(Nt);
        ind0=find(ps_t{k_s,k_t}>epsilon);% 取出非0
        KL_st{k_s,k_t}=sum((ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
    end
end
% 目标域树+源域类计算KL
for k_t=1:max(imt_gt(:))
    for k_s=1:max(ims_gt(:))
        ind=ims_gt==k_s;
        Xsk=normcols(ims(ind,:));
        pred_t=predict(model_t{k_t},Xsk);
        Ns = histcounts(pred_t,edget{k_t});
        Ns = conv(Ns,se_mask);
        pt_s{k_t,k_s}=Ns./sum(Ns);
        ind0=find(pt_s{k_t,k_s}>epsilon);
        KL_ts{k_t,k_s}=sum((pt_s{k_t,k_s}(ind0).*log(pt_s{k_t,k_s}(ind0)./pt_t{k_t}(ind0))));
    end
end
KL_stmat=cell2mat(KL_st);
KL_tsmat=cell2mat(KL_ts);
F1(kkk)=1/(abs(sum(KL_stmat(:))-sum(KL_tsmat(:)))/numel(KL_stmat)+1);
% Ew=0.5*(KL_stmat+KL_tsmat)
[a1,a2]=sort(cell2mat(KL_st));
[b1,b2]=sort(cell2mat(KL_ts));
tttt(kkk)=isequal(a2,b2);
end

% k_s=1;k_t=4;
% ind= im_gts==k_s;% tree
% Xsk=normcols(Xs(ind,:));
% IDX=kmeans(Xsk,50,'MaxIter',100);% kmeans
% modeli=fitctree(Xsk,IDX);
% pred_s=predict(modeli,Xsk);
% [Ns,edges] = histcounts(pred_s);
% ps_s=Ns./sum(Ns);
% 
% ind=im_gtt==k_t;
% Xtk=normcols(Xt(ind,:));
% pred_t=predict(modeli,Xtk);
% Nt = histcounts(pred_t,edges);
% ps_t=Nt./sum(Nt);
% KL_st=sum(ps_t.*log2(ps_t./ps_s));

