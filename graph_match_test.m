%% 
% A Novel Graph-Matching-Based Approach for Domain Adaptation in Classification of Remote Sensing Image Pair
load('E:\TransfLearning\area1_target\im1.mat', 'im')
load('E:\TransfLearning\area1_target\im1.mat', 'im_gt')
im_gt=im_gt+1;
Xt=reshape(im,[],size(im,3));

load('E:\TransfLearning\area7_source\im1.mat', 'im')
load('E:\TransfLearning\area7_source\im1.mat', 'im_gt')
im_gt=double(im_gt)+1;
Xs=reshape(im,[],size(im,3));
% % [y,model,mse] = knKmeans(double(Xs'),max(im_gt(:))+2,@knLin);
IDX = kmeans(Xs, max(im_gt(:))+2,'MaxIter',10000,'Display','final');
figure(1),hold on
LinS=['ro';'go';'bo';'yo';'mo';'ko';'co'];
for k=1:max(IDX)
    plot3(Xs(IDX==k,1),Xs(IDX==k,2),Xs(IDX==k,3),LinS(k,:));
end
hold off
pred=reshape(IDX,size(im_gt));
imtool(pred,[])

% figure(2),hold on
% for k=1:max(im_gt(:))
%     plot3(Xs(im_gt==k,1),Xs(im_gt==k,2),Xs(im_gt==k,3),LinS(k,:));
% end
% hold off
% 
% im11=im(1:floor(size(im,1)/2),:,:);
% im11_gt=im_gt(1:floor(size(im,1)/2),:);
% Xs=reshape(im11,[],size(im11,3));
% [y,model,mse] = knKmeans(double(Xs'),max(im_gt(:))+2,@knLin);

%% 把这个转成函数
P = [0.002 0.098 0.9];%训练验证测试样本比例
ind_gt= find(im_gt);
labels_gt= double(im_gt(ind_gt));
Num_classes=length(unique(labels_gt));
trainInd=[];valInd=[];testInd=[];
if length(P)==3
    if sum(P)~=1
        warning(['比例和不为1，自动把总体样本按比例划分：',...
            num2str(P(1)/sum(P)),':',num2str(P(2)/sum(P)),':',num2str(P(3)/sum(P))]);
    end
elseif length(P)==2
    warning(['按照各类训练、验证样本个数进行划分:',num2str(P(1)),':',num2str(P(2))]);
end
for k=1:Num_classes
    ind_temp = find(labels_gt==k);
    %样本划分
    if length(P)==3
        [trainInd_k,valInd_k,testInd_k] = dividerand(length(ind_temp),P(1),P(2),P(3));
        trainInd=[trainInd;ind_temp(trainInd_k)];
        valInd=[valInd;ind_temp(valInd_k)];
        testInd=[testInd;ind_temp(testInd_k)];
    elseif length(P)==2
        temp = randperm(length(ind_temp));
        ind_temp=ind_temp(temp);
        trainInd = [trainInd;ind_temp(1:P(1))];
        valInd = [valInd;ind_temp(P(1)+1:sum(P))];
        testInd= [testInd;ind_temp(sum(P)+1:end)];
    end
    %得到的Ind是标记样本的索引而非全图，可以利用标记数据块进行检验，是否越界
end
Xs=Xs(trainInd,:);
%% 
nbins=10;
for k=1:size(Xs,2)
    [Ns(k,:),edges] = histcounts(Xs(:,k),nbins);
    [Nt(k,:),edges] = histcounts(Xt(:,k),edges);
end


