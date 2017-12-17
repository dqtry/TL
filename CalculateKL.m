function [KL_stmat,KL_tsmat]=CalculateKL(Xa,gt_a,Xb,gt_b,num_bins)
% ����һάֱ��ͼ���Ʒֲ�����KLɢ�ȼ���
% Xa Xb �б�ʾ�������б�ʾ���� 0��1֮��ĸ�����
% gt_a gt_b ��ʾ���ڵ���غ�
% num_bins ֱ��ͼͳ��ʱ�ļ����

edges=linspace(0,1,num_bins+1);
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    for k2=1:max(gt_b(:))
        Nas=[];Nbs=[];
        for k3=1:size(Xa,2)
            tempa=find(gt_a==k1);
            [Na,~]=histcounts(Xa(tempa,k3),edges);
            Na=conv(Na,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nas=[Nas,Na+eps];
            tempb=find(gt_b==k2);
            [Nb,~]=histcounts(Xb(tempb,k3),edges);
            Nb=conv(Nb,se_mask,'same');% �Ƚ���ƽ��Ȼ������
            Nbs=[Nbs,Nb+eps];
        end
        Pa=Nas./sum(Nas);
        Pb=Nbs./sum(Nbs);
        KL_stmat(k1,k2)=sum(Pa.*log(Pa./Pb));
        KL_tsmat(k2,k1)=sum(Pb.*log(Pb./Pa));
    end
end

end