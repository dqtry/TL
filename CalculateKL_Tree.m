function [KL_stmat,KL_tsmat]=CalculateKL_Tree(Xa,gt_a,Xb,gt_b,num_clusters)
% 
se_mask=fspecial('gaussian',[1,3]);
epsilon=0;
for k_s=1:max(gt_a(:))
    ind= gt_a==k_s;% tree
    Xsk=normcols(Xa(ind,:));
    IDX=kmeans(Xsk,num_clusters,'MaxIter',10000);% kmeans
    model_s{k_s}=fitctree(Xsk,IDX);
    pred_s=predict(model_s{k_s},Xsk);
%     mean(IDX==pred_s);
    [Ns,edges{k_s}] = histcounts(pred_s);
    Ns = conv(Ns,se_mask);
    ps_s{k_s}=Ns./sum(Ns);
end
%% 构建目标域树及簇统计
for k_t=1:max(gt_b(:))
    ind= gt_b==k_t;% tree
    Xtk=normcols(Xb(ind,:));
    IDX=kmeans(Xtk,num_clusters,'MaxIter',10000);% kmeans
    model_t{k_t}=fitctree(Xtk,IDX);
    pred_t=predict(model_t{k_t},Xtk);
    [Nt,edget{k_t}] = histcounts(pred_t);
    Nt = conv(Nt,se_mask);
    pt_t{k_t}=Nt./sum(Nt);
end
%% 源域树+目标域簇计算KL
for k_s=1:max(gt_a(:))
    for k_t=1:max(gt_b(:))
        ind=gt_b==k_t;
        Xtk=normcols(Xb(ind,:));
        pred_t=predict(model_s{k_s},Xtk);
        Nt = histcounts(pred_t,edges{k_s});
        Nt = conv(Nt,se_mask);
        ps_t{k_s,k_t}=Nt./sum(Nt);
        ind0=find(ps_t{k_s,k_t}>epsilon);% 取出非0
        KL_st{k_s,k_t}=sum((ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
    end
end
% 目标域树+源域类计算KL
for k_t=1:max(gt_b(:))
    for k_s=1:max(gt_a(:))
        ind=gt_a==k_s;
        Xsk=normcols(Xa(ind,:));
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
