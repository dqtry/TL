%% 每类取num个进行predict
load('E:\TransfLearning\data300\DataForKeras.mat', 'y_test')
load('E:\TransfLearning\data300\DataForKeras.mat', 'x_test_t')
ind=[];
num=200;
for k=1:4
    temp=find(y_test==k);
    ind=[ind;temp(randperm(length(temp),num))];
end
y_test=y_test(ind);
x_test_t=x_test_t(:,:,:,ind);
plot(1:length(y_test),y_test)
%% 利用聚类得到的点作为源域样本
clear
load E:\TransfLearning\Source\im_r.mat im im_gt
ind_gt=find(im_gt);
label_gt=im_gt(ind_gt);
im2d=reshape(double(im),[],size(im,3));
load('E:\TransfLearning\data300\TVT_ind.mat', 'TrainingPointsCell','ValPointsCell')
load('E:\TransfLearning\MatlabCode2\res_kmedoids100_gt.mat', 'midx')
ind_chosen=[];
for i_class=1:length(unique(label_gt))
    ind_temp=find(im_gt==i_class);
    data=im2d(ind_temp,:);
    ind_chosen=[ind_chosen;ind_temp(midx{i_class})];
end
%% 取出训练样本进行predict
clear
PatchSize=9;
load('E:\TransfLearning\SourceToTarget\example2\TVT_ind.mat', 'TrainingPointsCell')
load E:\TransfLearning\Target\im1.mat im im_gt
PointsCell=TrainingPointsCell;
im_t=im;im_gt_t=im_gt;
ind_gt_t= find(im_gt_t);
labels_gt_t= single(im_gt_t(ind_gt_t));
numClasses=length(PointsCell);
ind_train_t=[]; y_train_t=[];% target
for k=1:numClasses
    ind_train_t = [ind_train_t;PointsCell{k,2}];
    y_train_t = [ y_train_t;k*ones(length(PointsCell{k,2}),1)];
end
X=reshape(single(im_t),[],3);
X = bsxfun(@minus, X, min(X));
X = bsxfun(@rdivide,X,max(X));
im_t=reshape(X,size(im_t));
IMAGES_t=GeneratePatch(im_t,PatchSize);%%4D斑块生成E:\Python\MatlabCode
IMAGES_t=IMAGES_t(:,:,:,ind_gt_t);
x_test_t=IMAGES_t(:,:,:,ind_train_t);
y_test_t=labels_gt_t(ind_train_t);%检测是否一致
save E:\TransfLearning\SourceToTarget\example2\DataTrain60.mat x_test_t y_test_t
