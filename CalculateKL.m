function [KL_stmat,KL_tsmat]=CalculateKL(Xa,gt_a,Xb,gt_b,num_bins)
% 利用一维直方图估计分布进行KL散度计算
% Xa Xb 行表示样本，列表示特征 0～1之间的浮点型
% gt_a gt_b 表示对于的类簇号
% num_bins 直方图统计时的间隔数

edges=linspace(0,1,num_bins+1);
se_mask=fspecial('gaussian',[1,3]);
for k1=1:max(gt_a(:))
    for k2=1:max(gt_b(:))
        Nas=[];Nbs=[];
        for k3=1:size(Xa,2)
            tempa=find(gt_a==k1);
            [Na,~]=histcounts(Xa(tempa,k3),edges);
            Na=conv(Na,se_mask,'same');% 先进行平滑然后连接
            Nas=[Nas,Na+eps];
            tempb=find(gt_b==k2);
            [Nb,~]=histcounts(Xb(tempb,k3),edges);
            Nb=conv(Nb,se_mask,'same');% 先进行平滑然后连接
            Nbs=[Nbs,Nb+eps];
        end
        Pa=Nas./sum(Nas);
        Pb=Nbs./sum(Nbs);
        KL_stmat(k1,k2)=sum(Pa.*log(Pa./Pb));
        KL_tsmat(k2,k1)=sum(Pb.*log(Pb./Pa));
    end
end

end