function [Ew,Ew2,eval]=CalculateEw_Tree(Xa,gt_a,Xb,gt_b,num_clusters)
% 
se_mask=fspecial('gaussian',[1,3]);
epsilon=0;
for k_s=1:max(gt_a(:))
    ind= gt_a==k_s;% tree
    Xsk=Xa(ind,:);
%     Xsk=normcols(Xa(ind,:));
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
    Xtk=Xb(ind,:);
%     Xtk=normcols(Xb(ind,:));
    IDX=kmeans(Xtk,num_clusters,'MaxIter',10000);% kmeans
    model_t{k_t}=fitctree(Xtk,IDX);
    pred_u=predict(model_t{k_t},Xtk);
    [Nt,edget{k_t}] = histcounts(pred_u);
    Nt = conv(Nt,se_mask);
    pt_t{k_t}=Nt./sum(Nt);
end
%% 
for k_s=1:max(gt_a(:))
    for k_t=1:max(gt_b(:))
        tempa=find(gt_a==k_s);
        tempb=find(gt_b==k_t);
        % source tree
        feat_u=[Xa(tempa,:);Xb(tempb,:)];
        pred_u=predict(model_s{k_s},feat_u);
        Nt = histcounts(pred_u,edges{k_s});
        Nt = conv(Nt,se_mask);
        ps_t{k_s,k_t}=Nt./sum(Nt);
        ind0=find(ps_t{k_s,k_t}>epsilon);% 取出非0
        KL_su{k_s,k_t}=sum((ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
        % target tree
        pred_u=predict(model_t{k_t},feat_u);
        Ns = histcounts(pred_u,edget{k_t});
        Ns = conv(Ns,se_mask);
        pt_s{k_t,k_s}=Ns./sum(Ns);
        ind0=find(pt_s{k_t,k_s}>epsilon);
        KL_tu{k_t,k_s}=sum((pt_s{k_t,k_s}(ind0).*log(pt_s{k_t,k_s}(ind0)./pt_t{k_t}(ind0))));
    end
end
KL_sumat=cell2mat(KL_su);
KL_tumat=cell2mat(KL_tu);
Ew=sqrt(0.5*KL_sumat+0.5*KL_tumat);
Ew2=sqrt(min(KL_sumat,KL_tumat));
Ewt=Ew+100*eye(size(Ew));%为了取得除了非对角线元素的最小值
% 非对角线元素的最小值与对角线元素的最大值之差，即类间距离的最小值与类内距离的最大值之差要大
eval(1)=min(Ewt(:))-max(diag(Ew));
Ewt2=Ew2+100*eye(size(Ew2));
eval(2)=min(Ewt2(:))-max(diag(Ew2));

%% 源域树+目标域簇计算KL
% for k_s=1:max(gt_a(:))
%     for k_t=1:max(gt_b(:))
%         tempa=find(gt_a==k_s);
%         tempb=find(gt_b==k_t);
%         feat_u=[Xa(tempa,:);Xb(tempb,:)];
%         pred_u=predict(model_s{k_s},feat_u);
%         Nt = histcounts(pred_u,edges{k_s});
%         Nt = conv(Nt,se_mask);
%         ps_t{k_s,k_t}=Nt./sum(Nt);
%         ind0=find(ps_t{k_s,k_t}>epsilon);% 取出非0
%         KL_su{k_s,k_t}=sum((ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
%     end
% end
% % 目标域树+源域类计算KL
% for k_t=1:max(gt_b(:))
%     for k_s=1:max(gt_a(:))
%         tempa=find(gt_a==k_s);
%         tempb=find(gt_b==k_t);
%         feat_u=[Xa(tempa,:);Xb(tempb,:)];
%         pred_u=predict(model_t{k_t},feat_u);
%         Ns = histcounts(pred_u,edget{k_t});
%         Ns = conv(Ns,se_mask);
%         pt_s{k_t,k_s}=Ns./sum(Ns);
%         ind0=find(pt_s{k_t,k_s}>epsilon);
%         KL_tu{k_t,k_s}=sum((pt_s{k_t,k_s}(ind0).*log(pt_s{k_t,k_s}(ind0)./pt_t{k_t}(ind0))));
%     end
% end