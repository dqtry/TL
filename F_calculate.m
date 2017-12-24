%% ����F���ۣ��õ���ѵľ������M
% clear
% rng(0);
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
Xa=normcols(Xa);
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
Xb=normcols(Xb);
[ims,imt]=pavia_adjust(Xa,Xb,100,0.001);
% [~,C,~,~,midx]=kmedoids(ims',10);
% ims=ims(:,midx);
% imt=imt(:,midx);
% load E:\TransfLearning\PUC\feats2.mat feat_a feat_b
% ims=feat_a;imt=feat_b;
% clear feat_a feat_b
% load E:\TransfLearning\PUC\Pu.mat Ugt
% gt_a=Ugt;
% load E:\TransfLearning\PUC\Pc.mat Cgt
% gt_b=Cgt;
%-----------------------------------------------------------------
inda=find(gt_a);indb=find(gt_b);%==1|gt_a==2|gt_a==6|gt_a==7 ==1|gt_b==2|gt_b==6|gt_b==7
ims=ims(inda,:);imt=imt(indb,:);
ims_gt=gt_a(inda);imt_gt=gt_b(indb);
M_vector=7*ones(5,1);%3:7 5:9
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
% load E:\TransfLearning\area1_target\im1.mat im %im_gt
% imt=reshape(double(im),[],size(im,3));%imt_gt=im_gt+1;
% % imt=ims;%% ����ͬ�����ݽ��в���
% % imt_gt=ims_gt;
% [cols,rows]=meshgrid(1:size(im_gt,2),1:size(im_gt,1));%�������������Ϣ���о���
F1=zeros(size(M_vector));F2=zeros(size(M_vector));F=zeros(size(M_vector));
num_bins=100;
for k=1:length(M_vector)
    M=M_vector(k);
    %---------------------------------------------------------------------
    id_cluster = kmeans(imt,M,'MaxIter',10000,'OnlinePhase','on','Replicates',4,'Options',statset('UseParallel',1),'Display','final');%����
    %---------------------------------------------------------------------
    %% KL divergence
    % һάֱ��ͼ���ֲ��Σ�ÿ������100��bin����ͳ�ƣ�ƽ��֮�����һά����
    [KL_stmat,KL_tsmat]=CalculateKL(ims,ims_gt,imt,id_cluster,num_bins);%Դ���ࡢĿ�����
    F1(k)=1/(abs(sum(KL_stmat(:))-sum(KL_tsmat(:)))/numel(KL_stmat)+1);
    % ------����Ew����ƥ�䣬���Ծ���Ч��
    [Ew,Ew2,Ew3,eval]=CalculateEw(ims,ims_gt,imt,id_cluster,num_bins);
%     pairs=iter_match(Ew2);
%     matched_pair=cell2mat(pairs);
%     [~,I] = sort(matched_pair(1,:));%% ������sortrows
%     real_matched=matched_pair(2,I);
%     match_id=zeros(size(id_cluster));
%     match_id2=zeros(size(id_cluster));
%     %mode(id_cluster(imt_gt==kkk))
%     for kkk=1:max(id_cluster(:))
%         match_id(id_cluster==real_matched(kkk))=kkk;
%         %real_matched��Ŀ�����Ӧ�ı�ǣ�����Ew��Դ���Ŀ����ƥ���Ͻ�����֤
%         match_id2(id_cluster==mode(id_cluster(imt_gt==kkk)))=kkk;
%     end
%     valid(k)=mean(match_id==imt_gt);
%     valid2(k)=mean(match_id2==imt_gt);
    %% ������ںʹؼ䷽��
    imt_gray=mean(imt,2);%�Ҷ�ͼ������ںʹؼ䷽��
    intra_v=0;inter_v=0;
    N=size(imt_gray,1);mu=[];
    for k1=1:M
        temp=find(id_cluster==k1);
        Cluster=imt_gray(temp,:);
        mu(k1,:)=mean(Cluster,1);
        n_k1=length(temp);
        omega(k1)=n_k1/N;
        intra_v=intra_v+sum((Cluster-mu(k1,:)).^2,1)./N;
    end
    mu_T=mean(imt_gray,1);% mu_T=mean(X,1);mu_T=sum(omega'.*mu);
    sigma_T=mean((imt_gray-mu_T).^2,1);
    % for k2=1:M-1
    %     for k3=k2+1:M
    %         inter_v=inter_v+omega(k2)*omega(k3)*(mu(k2,:)-mu(k3,:)).^2;
    %     end
    % end
    for k4=1:M
        inter_v=inter_v+omega(k4)*(mu(k4,:)-mu_T).^2;
    end
    F2(k)=sum(inter_v)+1/(sum(intra_v)+1);%+1?????????????
    F(k)=F1(k)+F2(k);
end
% plot(M_vector,F,'*-');
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

