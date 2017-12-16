%% 类内方差、类间方差和灰度级的总方差
% https://www.cnblogs.com/xiaomanon/p/4110006.html
load('E:\TransfLearning\area7_source\im1.mat', 'im')
X=reshape(double(im),[],size(im,3));
F2=[];
for M=4%:8
    idx = kmeans(X,M);
    intra_v=0;inter_v=0;
    N=size(X,1);mu=[];
    for k1=1:M
        temp=find(idx==k1);
        Cluster=X(temp,:);
        mu(k1,:)=mean(Cluster,1);
        n_k1=length(temp);
        omega(k1)=n_k1/N;
        intra_v=intra_v+sum((Cluster-mu(k1,:)).^2,1)./N;
    end
    mu_T=mean(X,1);% mu_T=mean(X,1);mu_T=sum(omega'.*mu);
    sigma_T=mean((X-mu_T).^2,1);
    % for k2=1:M-1
    %     for k3=k2+1:M
    %         inter_v=inter_v+omega(k2)*omega(k3)*(mu(k2,:)-mu(k3,:)).^2;
    %     end
    % end
    for k4=1:M
        inter_v=inter_v+omega(k4)*(mu(k4,:)-mu_T).^2;
    end
    F2(end+1)=sum(inter_v)+1/(sum(intra_v)+1);
    inter_v+intra_v-sigma_T%判断计算正确
end
