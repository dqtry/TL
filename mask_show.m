%% 单类查看，检查类别是否对应
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
ims=reshape(Pu_same,[],size(Pu_same,3));ims_gt=Ugt;
ims=normcols(ims);

load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
imt=reshape(Pc_same,[],size(Pc_same,3));imt_gt=Cgt;
imt=normcols(imt);

midx=[52 80 23];% kmedoids得到（顺序调整）的三个显示效果较好的三个波段
a=ims(:,midx);
b=imt(:,midx);

stat=[];%统计各类个数以及总数
for id=1:7
    ind1=find(ims_gt==id);
    stat(id,1)=length(ind1);
    ind2=find(imt_gt==id);
    stat(id,2)=length(ind2);
    aa=ones(size(a));
    bb=ones(size(b));
    aa(ind1,:)=a(ind1,:);
    bb(ind2,:)=b(ind2,:);
    figure(id);
    subplot(1,2,1);imshow(reshape(aa,[size(ims_gt),3]),[])
    subplot(1,2,2);imshow(reshape(bb,[size(imt_gt),3]),[])
end
stat(id+1,:)=sum(stat);
%% 聚类效果测试（用了真实标记）
match_id=zeros(size(id_cluster));
for kkk=1:max(id_cluster(:))
    match_id(id_cluster==mode(id_cluster(imt_gt==kkk)))=kkk;
end
mean(match_id==imt_gt)
    