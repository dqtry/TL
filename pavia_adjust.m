function [Xa,Xb]=pavia_adjust(Xa,Xb,nbins,thresh)
% 原始PaviaU和PaviaC数据大部分分布在灰度值较小的部分，分波段进行拉伸

% load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
% Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
% Xa=normcols(Xa);
% load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
% Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
% Xb=normcols(Xb);
% nbins=100;
% thresh=0.001;
aup=zeros(1,size(Xa,2));bup=zeros(1,size(Xb,2));
for k=1:size(Xa,2)
    [N,edges] = histcounts(Xa(:,k),nbins);
    cumN=cumsum(N,'reverse')./sum(N);
    temp=find(cumN<thresh, 1, 'first');
    aup(k)=edges(temp);%temp-1
    Xa(:,k)=imadjust(Xa(:,k),[min(Xa(:,k)),aup(k)],[0,1]);
    
    [N,edges] = histcounts(Xb(:,k),nbins);
    cumN=cumsum(N,'reverse')./sum(N);
    temp=find(cumN<thresh, 1, 'first');
    bup(k)=edges(temp);%temp-1
    Xb(:,k)=imadjust(Xb(:,k),[min(Xb(:,k)),bup(k)],[0,1]);
end