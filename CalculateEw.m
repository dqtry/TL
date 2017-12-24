function [Ew1,Ew2,Ew3,eval]=CalculateEw(Xa,gt_a,Xb,gt_b,num_bins)
% Ew1��δ��ȫ���JSɢ�ȵĲ�׼ȷ�Ľ����a��bֱ�Ӻϲ�δƽ���������ڴ���1��ֵ
% Ew2��Ew3���ݶ��壨KLɢ�ȣ��Լ���ʽ�Ƶ������ݷֲ����õ���
% eval���Զ�������ָ�꣬���Ǻܿɿ����ɺ���
%% ֱ��ͼͳ�Ʒֲ������ڼ���Ew(JSɢ�ȵ�ƽ����)
% ��Դ��ÿ��ֳ������ֽ�����֤����ֱ��ͼ�ʺ��ڼ���Ew
% A Novel Graph-Matching-Based Approach for Domain Adaptation in Classification of Remote Sensing Image Pair
% Ew(JSɢ�ȵ�ƽ����)��һ���ԳƵľ��룬���ڱ�ʾ��Ȩ����ͼ�е�Ȩ�أ���KLɢ�ȵõ�
% KLɢ��������أ��������ؼ�ȥ���۷ֲ�����Ϣ�أ����Ǹ����ǶԳƣ�
% Xa Xb �б�ʾ�������б�ʾ���� 0��1֮��ĸ�����
% gt_a gt_b ��ʾ���ڵ���غ�
% num_bins ֱ��ͼͳ��ʱ�ļ����

%% ==============ͬһ��ͼ�ֳ������ֽ��в���==================
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
% % figure,hold on
% % for kkk=1:5
% % rng();
% Xa=[];Xb=[];gt_a=[];gt_b=[];
% for k1=1:length(unique(ims_gt))
%     temp = find(ims_gt==k1);
%     temp = temp(randperm(length(temp)));
%     temp_feat=ims(temp,:);
%     Xa=[Xa;temp_feat(1:round(length(temp)/2),:)];
%     gt_a=[gt_a;k1*ones(round(length(temp)/2),1)];
%     Xb=[Xb;temp_feat(round(length(temp)/2)+1:end,:)];
%     gt_b=[gt_b;k1*ones(length(temp)-round(length(temp)/2),1)];
% end
% num_bins=50;
% ֱ��ͼͳ����ý���ƽ������
%% ==============����ͼ���в���================
% load E:\TransfLearning\area7_source\im1.mat im im_gt
% Xa=reshape(im,[],size(im,3));gt_a=im_gt+1;
% load E:\TransfLearning\area1_target\im1.mat im im_gt
% Xb=reshape(im,[],size(im,3));gt_b=im_gt+1;
% num_bins=50;
%% 1.һάֱ��ͼ���ֲ��Σ�ÿ������100��bin����ͳ�ƣ�ƽ��֮�����һά����
% num_bins=[11 21 31 41 51 61 71 81 91 101 111 121 131 141 151];
% % for kt=1:length(num_bins)
edges=linspace(min(Xa(:)),max(Xa(:)),num_bins+1);%������Ӧ�ù�һ��0~1
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    for k2=1:max(gt_b(:))
        tempa=find(gt_a==k1);
        tempb=find(gt_b==k2);
        feat_u=[Xa(tempa,:);Xb(tempb,:)];
        Nas=[];Nbs=[];Nus=[];
        for k3=1:size(Xa,2)
            [Na,~]=histcounts(Xa(tempa,k3),edges);
            Na=conv(Na,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nas=[Nas,Na+eps];
            [Nb,~]=histcounts(Xb(tempb,k3),edges);
            Nb=conv(Nb,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nbs=[Nbs,Nb+eps];
            [Nu,~]=histcounts(feat_u(:,k3),edges);
            Nu=conv(Nu,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nus=[Nus,Nu+eps];
        end
        Pa=Nas./sum(Nas);Pb=Nbs./sum(Nbs);Pu=Nus./sum(Nus);%����
        KL_au(k1,k2)=sum(Pa.*log(Pa./Pu));
        KL_bu(k1,k2)=sum(Pb.*log(Pb./Pu));
        Ew1(k1,k2)=sqrt(0.5*KL_au(k1,k2)+0.5*KL_bu(k1,k2));
%         Ew(k1,k2)=sqrt(max(KL_au(k1,k2),KL_bu(k1,k2)));
%         Ew2(k1,k2)=sqrt(min(KL_au(k1,k2),KL_bu(k1,k2)));
        Pu2=0.5*Pa+0.5*Pb;
        KL_au2(k1,k2)=sum(Pa.*log(Pa./Pu2));
        KL_bu2(k1,k2)=sum(Pb.*log(Pb./Pu2));
        Ew2(k1,k2)=sqrt(0.5*KL_au2(k1,k2)+0.5*KL_bu2(k1,k2));
        Ew3(k1,k2)=sqrt(0.5*(sum(Pa.*log(Pa))+sum(Pb.*log(Pb)))-sum(0.5*(Pa+Pb).*log(0.5*(Pa+Pb))));
    end
end
Ewt=Ew1+100*eye(size(Ew1));%Ϊ��ȡ�ó��˷ǶԽ���Ԫ�ص���Сֵ
% �ǶԽ���Ԫ�ص���Сֵ��Խ���Ԫ�ص����ֵ֮������������Сֵ�����ھ�������ֵ֮��Ҫ��
eval(1)=min(Ewt(:))-max(diag(Ew1));
Ewt2=Ew2+100*eye(size(Ew2));
eval(2)=min(Ewt2(:))-max(diag(Ew2));
% % Ew2=Ew+100*eye(size(Ew));
% % eval(kt)=min(Ew2(:))-max(diag(Ew));
% % end
% % plot(num_bins,eval)
% % end
% % hold off