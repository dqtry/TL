%% 三维直方图统计计算KL散度
%
load E:\TransfLearning\area7_source\im1.mat im im_gt
ims=reshape(im,[],size(im,3));ims_gt=im_gt+1;
edges1=linspace(0,1,11);
se_mask=zeros(3,3,3);
for kse=1:numel(se_mask)
    [a,b,c]=ind2sub(size(se_mask),kse);
    se_mask(kse)=sqrt(sum(abs([a,b,c]-[2,2,2])));
end
se_mask=max(se_mask(:))-se_mask+eps;
se_mask=se_mask./sum(se_mask(:));
for k=1:length(unique(ims_gt))
temps = find(ims_gt==k);
X=ims(temps,:);
hist_rgb{k}=zeros(length(edges1)-1,length(edges1)-1,length(edges1)-1,'single');
for k1=1:length(edges1)-1
    for k2=1:length(edges1)-1
        for k3=1:length(edges1)-1
            if k1~=length(edges1-1)&& k2~=length(edges1-1) &&  k3~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<edges1(k3+1));
            elseif k1==length(edges1-1) && k2~=length(edges1-1) && k3~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<=edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<edges1(k3+1));
            elseif k2==length(edges1-1) && k1~=length(edges1-1) && k3~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<=edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<edges1(k3+1));
            elseif k3==length(edges1-1) && k1~=length(edges1-1) && k2~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<=edges1(k3+1));
            elseif k1==length(edges1-1) && k2==length(edges1-1) && k3~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<=edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<=edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<edges1(k3+1));
            elseif k1==length(edges1-1) && k3==length(edges1-1) && k2~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<=edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<=edges1(k3+1));
            elseif k2==length(edges1-1) && k3==length(edges1-1) && k1~=length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<=edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<=edges1(k3+1));
            elseif k1==length(edges1-1) && k2==length(edges1-1) && k3==length(edges1-1)
                hist_rgb{k}(k1,k2,k3)=sum(X(:,1)>=edges1(k1) & X(:,1)<=edges1(k1+1) & X(:,2)>=edges1(k2) & X(:,2)<=edges1(k2+1) & X(:,3)>=edges1(k3) & X(:,3)<=edges1(k3+1));
            end
        end
    end
end
hist_rgb{k}=hist_rgb{k}./sum(hist_rgb{k}(:));
temp=convn(hist_rgb{k}+eps,se_mask,'same');
distribut(k,:)=temp(:)/sum(temp(:));
end
for k11=1:size(distribut,1)
    for k22=1:size(distribut,1)
        KL_div(k11,k22)=sum(distribut(k11,:).*log(distribut(k11,:)./distribut(k22,:)));
    end
end

%% 查看统计情况
% aa=[];
% for k4=1:length(hist_rgb)
%     aa(1,k4)=length(find(hist_rgb{1,k4}>0.1));
%     aa(2,k4)=length(find(hist_rgb{1,k4}>0.01));
%     aa(3,k4)=length(find(hist_rgb{1,k4}>0.001));
%     aa(4,k4)=length(find(hist_rgb{1,k4}>0));
% end
    
