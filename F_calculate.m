%% ����F���ۣ��õ���ѵľ������M

load E:\TransfLearning\area1_target\im1.mat im im_gt
imt=reshape(double(im),[],size(im,3));%imt_gt=im_gt+1;
[cols,rows]=meshgrid(1:size(im_gt,2),1:size(im_gt,1));
F1=[];F2=[];F=[];
rng(0);
for M=3:8
id_cluster = kmeans(imt,M,'MaxIter',10000);
% imt=mean(imt,2);%�Ҷ�ͼ������ںʹؼ䷽��
intra_v=0;inter_v=0;
N=size(imt,1);mu=[];
for k1=1:M
    temp=find(id_cluster==k1);
    Cluster=imt(temp,:);
    mu(k1,:)=mean(Cluster,1);
    n_k1=length(temp);
    omega(k1)=n_k1/N;
    intra_v=intra_v+sum((Cluster-mu(k1,:)).^2,1)./N;
end
mu_T=mean(imt,1);% mu_T=mean(X,1);mu_T=sum(omega'.*mu);
sigma_T=mean((imt-mu_T).^2,1);
% for k2=1:M-1
%     for k3=k2+1:M
%         inter_v=inter_v+omega(k2)*omega(k3)*(mu(k2,:)-mu(k3,:)).^2;
%     end
% end
for k4=1:M
    inter_v=inter_v+omega(k4)*(mu(k4,:)-mu_T).^2;
end
F2(end+1)=sum(inter_v)+1/(sum(intra_v)+1);%%%%%%%%+1
%% 1.һάֱ��ͼ���ֲ��Σ�ÿ������100��bin����ͳ�ƣ�ƽ��֮�����һά����
load E:\TransfLearning\area7_source\im1.mat im im_gt
ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
load E:\TransfLearning\area1_target\im1.mat im %im_gt
imt=reshape(double(im),[],size(im,3));%imt_gt=im_gt+1;
% imt=ims;%% ����ͬ�����ݽ��в���
% imt_gt=ims_gt;

num_bins=50;
[KL_stmat,KL_tsmat]=CalculateKL(ims,ims_gt,imt,id_cluster,num_bins);%Դ���ࡢĿ�����
F1(end+1)=1/(abs(sum(KL_stmat(:))-sum(KL_tsmat(:)))/numel(KL_stmat)+1);%%%%%%%%%%%+1
F(end+1)=F1(end)+F2(end);
end
%% ���÷Ǽල���õ�KLɢ�ȣ��Ӷ�����F1
% load('E:\TransfLearning\area7_source\im1.mat', 'im')
% load('E:\TransfLearning\area7_source\im1.mat', 'im_gt')
% im_gts=im_gt+1;
% Xs=reshape(im,[],size(im,3));
% load('E:\TransfLearning\area1_target\im1.mat', 'im')
% Xt=reshape(double(im),[],size(im,3));
% epsilon=0;
% % ����Դ��������ͳ��
% for k_s=1:max(im_gts(:))
%     ind= im_gts==k_s;% tree
%     Xsk=normcols(Xs(ind,:));
%     IDX=kmeans(Xsk,50,'MaxIter',10000);% kmeans
%     model_s{k_s}=fitctree(Xsk,IDX);
%     pred_s=predict(model_s{k_s},Xsk);
% %     mean(IDX==pred_s)
%     [Ns,edges{k_s}] = histcounts(pred_s);
%     ps_s{k_s}=Ns./sum(Ns);
% end
% % ����Ŀ����������ͳ��
% for k_t=1:max(id_cluster(:))
%     ind= id_cluster==k_t;% tree
%     Xtk=normcols(Xt(ind,:));
%     IDX=kmeans(Xtk,50,'MaxIter',10000);% kmeans
%     model_t{k_t}=fitctree(Xtk,IDX);
%     pred_t=predict(model_t{k_t},Xtk);
%     [Nt,edget{k_t}] = histcounts(pred_t);
%     pt_t{k_t}=Nt./sum(Nt);
% end
% % Դ����+Ŀ����ؼ���KL
% for k_s=1:max(im_gts(:))
%     for k_t=1:max(id_cluster(:))
%         ind=id_cluster==k_t;
%         Xtk=normcols(Xt(ind,:));
%         pred_t=predict(model_s{k_s},Xtk);
%         Nt = histcounts(pred_t,edges{k_s});
%         ps_t{k_s,k_t}=Nt./sum(Nt);
%         ind0=find(ps_t{k_s,k_t}>epsilon);% ȡ����0
%         KL_st{k_s,k_t}=sum((ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
% %         KL_st2{k_s,k_t}=sum(abs(ps_t{k_s,k_t}(ind0).*log(ps_t{k_s,k_t}(ind0)./ps_s{k_s}(ind0))));
%     end
% end
% % Ŀ������+Դ�������KL
% for k_t=1:max(id_cluster(:))
%     for k_s=1:max(im_gts(:))
%         ind=im_gts==k_s;
%         Xsk=normcols(Xs(ind,:));
%         pred_t=predict(model_t{k_t},Xsk);
%         Ns = histcounts(pred_t,edget{k_t});
%         pt_s{k_t,k_s}=Ns./sum(Ns);
%         ind0=find(pt_s{k_t,k_s}>epsilon);
%         KL_ts{k_t,k_s}=sum((pt_s{k_t,k_s}(ind0).*log(pt_s{k_t,k_s}(ind0)./pt_t{k_t}(ind0))));
% %         KL_ts2{k_t,k_s}=sum(abs(pt_s{k_t,k_s}(ind0).*log(pt_s{k_t,k_s}(ind0)./pt_t{k_t}(ind0))));
%     end
% end
% KL_stmat=cell2mat(KL_st);
% KL_tsmat=cell2mat(KL_ts);
% F1(end+1)=1/(abs(sum(KL_stmat(:))-sum(KL_tsmat(:)))/numel(KL_stmat)+1);%%%%%%%%%%%+1

    