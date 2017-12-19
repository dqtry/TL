function pairs=iter_match(Ew)
% pairs=iter_match(Ew)
% 根据Ew矩阵迭代查找匹配类
% Ew: Matrix contains sqrt of JS divergence
% pairs: Cell records sequential iteration

% load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
% Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
% load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
% Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
% Xa=normcols(Xa);
% Xb=normcols(Xb);
% nbins=100;
% thresh=0.01;
% aup=zeros(1,size(Xa,2));bup=zeros(1,size(Xb,2));
% for k=1:size(Xa,2)
%     [N,edges] = histcounts(Xa(:,k),nbins);
%     cumN=cumsum(N,'reverse')./sum(N);
%     temp=find(cumN<thresh, 1, 'first');
%     aup(k)=edges(temp);%temp-1
%     Xa(:,k)=imadjust(Xa(:,k),[min(Xa(:,k)),aup(k)],[0,1]);
%     
%     [N,edges] = histcounts(Xb(:,k),nbins);
%     cumN=cumsum(N,'reverse')./sum(N);
%     temp=find(cumN<thresh, 1, 'first');
%     bup(k)=edges(temp);%temp-1
%     Xb(:,k)=imadjust(Xb(:,k),[min(Xb(:,k)),bup(k)],[0,1]);
% end
% if size(Xa,2)~=size(Xb,2)
%     error('特征维度必须一致');
% end
% num_bins=100;
% [Ew,Ew2,eval]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins);%% 计算Ew距离
inda=1:size(Ew,1);indb=1:size(Ew,2);
n_iter=1;
while ~isempty(Ew)
    ks1=[];ks2=[];matcht=[];
    for k1=1:size(Ew,1)
        for k2=1:size(Ew,2)
            if Ew(k1,k2)==min(Ew(k1,:))&& Ew(k1,k2)==min(Ew(:,k2))
                matcht=[matcht,[inda(k1);indb(k2)]];
                ks1=[ks1,k1];ks2=[ks2,k2];
            end
        end
    end
    pairs{n_iter}=matcht;
    inda(ks1)=[];indb(ks2)=[];%保留原始类别号
    Ew(ks1,:)=[];Ew(:,ks2)=[];%更新Ew
    n_iter=n_iter+1;
end
% matched_pairs=cell2mat(pairs);

