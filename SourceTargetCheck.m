%% 检查源域和目标域的分布变化情况
%% ISPRS
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
% % % for nb=1:size(ims,2)
% % %     figure,histogram(ims(:,nb));
% % % end
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% imt=reshape(im,[],size(im,3));imt_gt=im_gt+1;
% for k=1:length(unique(ims_gt))
%     temps = find(ims_gt==k);
%     tempt = find(imt_gt==k);
%     figure(1),hold on
%     scatter(ims(temps,1),ims(temps,3));
%     figure(2),hold on
%     scatter(imt(tempt,1),imt(tempt,3));
% end
%% PU和PC
% midx=[52 80 23];
% load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
% Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
% load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
% Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
% Xa=normcols(Xa);
% Xb=normcols(Xb);
% % for kk=1:size(Xb,2)
% %     Xa(:,kk)=histeq(Xa(:,kk),256);%均衡化效果不好
% %     Xb(:,kk)=histeq(Xb(:,kk),256);
% % end
% [~, Xa] = pca(Xa,'Centered',true,'NumComponents',3);%,'NumComponents',5
% %PCA效果也不好
% % Xa=Xa(:,midx);
% Xa=normcols(Xa);
% [~, Xb] = pca(Xb,'Centered',true,'NumComponents',3);%,'NumComponents',5
% % Xb=Xb(:,midx);
% Xb=normcols(Xb);
% 
% edges=linspace(0,1,101);
% se_mask=fspecial('gaussian',[1,3]);
% figure;hold on
% for k1=1:max(gt_a(:))
%     tempa=find(gt_a==k1);
%     [Na,~]=histcounts(Xa(tempa,1),edges);
%     Na=conv(Na,se_mask,'same');% 先进行平滑然后连接
%     tempb=find(gt_b==k1);
%     [Nb,~]=histcounts(Xb(tempb,1),edges);
%     Nb=conv(Nb,se_mask,'same');% 先进行平滑然后连接
%     plot(edges(1:end-1),Na,'r-',edges(1:end-1),Nb,'g-');%,edges(1:end-1),Nb
%     pause(2);
% end
%% 计算Ew
% num_bins=50;
% [Ew,Ew2,eval]=CalculateEw(ims,ims_gt,imt,imt_gt,num_bins);%% 计算Ew距离
%% 域内分成两部分进行验证JS距离有效
% rng(0);% 设置种子点，可复现
% feat1=[];feat2=[];gt1=[];gt2=[];
% for k1=1:max(gt_b(:))
%     temp = find(gt_b==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=Xb(temp,:);
%     feat1=[feat1;temp_feat(1:round(length(temp)/2),:)];
%     gt1=[gt1;k1*ones(round(length(temp)/2),1)];
%     feat2=[feat2;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt2=[gt2;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
% [Ew,Ew2,eval]=CalculateEw(feat1,gt1,feat2,gt2,100);%% 计算Ew距离
%% ---------PU和PC显示效果检测-----------
% % [~,C,~,~,midx]=kmedoids(ims',3);% midx对应RGB可能需要调整
% midx=[52 80 23];
% imchoset=reshape(Xb,size(Pc_same));
% imshow(imchoset(:,:,midx),[])
% imchoses=reshape(Xa,size(Pu_same));
% figure,imshow(imchoses(:,:,midx),[])
%% 检查源域和目标域分布情况
% figure(1)
% xlim([0 1]);ylim([0 1]);
% figure(2)
% xlim([0 1]);ylim([0 1]);
% for k=1:max(gt_a(:))% gt中含有0
%     temps = find(gt_a==k);
%     figure(1),hold on
%     scatter(Xa(temps,1),Xa(temps,2));
%     tempt = find(gt_b==k);
%     figure(2),hold on
%     scatter(Xb(tempt,1),Xb(tempt,2));
% end
%% PC PU统计拉伸
% midx=[52 80 23];
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
Xa=normcols(Xa);
Xb=normcols(Xb);
% % 灰度拉伸
nbins=100;
thresh=0.001;
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
% kmedoids取波段
[~,~,~,~,midx]=kmedoids(Xa',10);% midx对应RGB可能需要调整
Xa=normcols(Xa(:,midx));
Xb=normcols(Xb(:,midx));
% imshow(reshape(Xa(:,midx),[size(gt_a),3]))

if size(Xa,2)~=size(Xb,2)
    error('特征维度必须一致');
end
% [COEFF, Xa, LATENT] = pca(Xa,'Centered',false,'NumComponents',10);%,'NumComponents',5
% Xa=normcols(Xa);
% [~, Xb] = pca(Xb,'Centered',false,'NumComponents',10);%,'NumComponents',5
% Xb=normcols(Xb);
% [~,C,~,~,midx]=kmedoids(Xa',3);% midx对应RGB可能需要调整
num_bins=100;
% [KL_stmat1,KL_tsmat1]=CalculateKL(Xa,gt_a,Xb,gt_b,num_bins);
[Ew,Ew2,eval]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins);%% 计算Ew距离
%% % 自编码特征
load E:\TransfLearning\PUC\Pu.mat Ugt
gt_a=Ugt;
load E:\TransfLearning\PUC\Pc.mat Cgt
gt_b=Cgt;feata=[];featb=[];
% load E:\TransfLearning\PUC\feats21.mat feat_a feat_b
% feata=feat_a;featb=feat_b;
% load E:\TransfLearning\PUC\feats22.mat feat_a feat_b
% feata=[feata,feat_a];featb=[featb,feat_b];
load E:\TransfLearning\PUC\feats23.mat feat_a feat_b
feata=[feata,feat_a];featb=[featb,feat_b];
num_bins=100;
[Ew,Ew2,eval]=CalculateEw(feata,gt_a,featb,gt_b,num_bins);%% 计算Ew距离
pairs=iter_match(Ew);%迭代匹配
matched_pairs=cell2mat(pairs);
%% 每类各波段均值绘制折线，查看类别可分性
load E:\TransfLearning\PUC\Pu.mat Pu_same Ugt
Xa=reshape(Pu_same,[],size(Pu_same,3));gt_a=Ugt;
load E:\TransfLearning\PUC\Pc.mat Pc_same Cgt
Xb=reshape(Pc_same,[],size(Pc_same,3));gt_b=Cgt;
Xa=normcols(Xa);
Xb=normcols(Xb);
figure(1),hold on
figure(2),hold on
figure(3),hold on
linecolor=[0.9961,0,0;
    0,0.5430,0;
    0,0,0.9961;
    0,0.9961,0.9961;
    0.9961,0,0.9961;
    0.9961,0.9961,0;
    0.542977035172045,0.109376668955520,0.382818341344320];
for kk=1:max(gt_a(:))
    tempa = find(gt_a==kk);
    amean(kk,:)=mean(Xa(tempa,:),1);
    tempb = find(gt_b==kk);
    bmean(kk,:)=mean(Xb(tempb,:),1);
    cmean(kk,:)=mean([Xa(tempa,:);Xb(tempb,:)],1);
    figure(1);
    plot(1:102,amean(kk,:),'color',linecolor(kk,:));
    figure(2)
    plot(1:102,bmean(kk,:),'color',linecolor(kk,:));
    figure(3)
    plot(1:102,cmean(kk,:),'color',linecolor(kk,:));
end
% figure(1);legend('树','沥青路','自封砖','沥青浇筑','草地','裸土','阴影');title('PU');
% figure(2);legend('树','沥青路','自封砖','沥青浇筑','草地','裸土','阴影');title('PC');
% figure(3);legend('树','沥青路','自封砖','沥青浇筑','草地','裸土','阴影');title('PU+PC');
% 
% amean=amean./sum(amean,2);
% bmean=bmean./sum(bmean,2);
% cmean=cmean./sum(cmean,2);
% se_mask=fspecial('gaussian',[1,3]);
% for k1=1:7
%     for k2=1:7
%         Pa=amean(k1,:);Pb=bmean(k2,:);
%         Pa=conv(Pa,se_mask,'same');% 先进行平滑然后连接
%         Pb=conv(Pb,se_mask,'same');% 先进行平滑然后连接
%         Pu=0.5*(Pa+Pb);
%         Pu=conv(Pu,se_mask,'same');% 先进行平滑然后连接
%         KL_au=sum(Pa.*log(Pa./Pu));
%         KL_bu=sum(Pb.*log(Pb./Pu));
%         Ew(k1,k2)=sqrt(0.5*KL_au+0.5*KL_bu);
%     end
% end



