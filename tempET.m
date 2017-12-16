%% 计算多尺度分割中的LV
clear
seg_sv=dlmread('E:\TransfLearning\results\size50.csv',';',1,1);
[num_obj,num_layer]=size(seg_sv);
seg_sv=sum(seg_sv,1)/num_obj;
sum(seg_sv)/num_layer


x=segimr2(1:2:end);
y=segimr2(2:2:end);
plot(x,y,'*-');
figure;
plot(x(2:end),diff(y)./y(1:end-1)*100,'*-')

%% 源域光谱特征进行SVM训练
