function [Ew,Ew2,eval]=CalculateEw2(Xa,gt_a,Xb,gt_b,num_bins)
%% ����JS�õ�Ew
edges=linspace(0,1,num_bins+1);
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    for k2=1:max(gt_b(:))
        tempa=find(gt_a==k1);
        tempb=find(gt_b==k2);
        Nas=[];Nbs=[];Nus=[];
        for k3=1:size(Xa,2)
            [Na,~]=histcounts(Xa(tempa,k3),edges);
            Na=conv(Na,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nas=[Nas,Na+eps];
            [Nb,~]=histcounts(Xb(tempb,k3),edges);
            Nb=conv(Nb,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nbs=[Nbs,Nb+eps];
            [Nu]=0.5*Na+0.5*Nb;
            Nu=conv(Nu,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nus=[Nus,Nu+eps];
        end
        Pa=Nas./sum(Nas);Pb=Nbs./sum(Nbs);Pu=Nus./sum(Nus);%����
        KL_au(k1,k2)=sum(Pa.*log(Pa./Pu));
        KL_bu(k1,k2)=sum(Pb.*log(Pb./Pu));
        Ew(k1,k2)=sqrt(0.5*KL_au(k1,k2)+0.5*KL_bu(k1,k2));
%         Ew(k1,k2)=sqrt(max(KL_au(k1,k2),KL_bu(k1,k2)));
        Ew2(k1,k2)=sqrt(min(KL_au(k1,k2),KL_bu(k1,k2)));
    end
end
Ewt=Ew+100*eye(size(Ew));%Ϊ��ȡ�ó��˷ǶԽ���Ԫ�ص���Сֵ
% �ǶԽ���Ԫ�ص���Сֵ��Խ���Ԫ�ص����ֵ֮������������Сֵ�����ھ�������ֵ֮��Ҫ��
eval(1)=min(Ewt(:))-max(diag(Ew));
Ewt2=Ew2+100*eye(size(Ew2));
eval(2)=min(Ewt2(:))-max(diag(Ew2));