clear;
clc;

load 'mean_HS_norm.mat';
Xc = double(PCmean');
Xu = double(PUmean');
clear PCmean;
clear PUmean;
NumBand = 7;
Xcb = atand(5*(Xc(2:end,:) - Xc(1:end-1,:)));
Xub = atand(5*(Xu(2:end,:) - Xu(1:end-1,:)));
D1 = 10000*ones(NumBand,NumBand);
D2 = 10000*ones(NumBand,NumBand);
D3 = 10000*ones(NumBand,NumBand);
D4 = 10000*ones(NumBand,NumBand);
for k1 = 1:NumBand
    for k2 = 1:NumBand
        D1(k1,k2) = mean(abs(Xc(:,k1) - Xu(:,k2)));
        D2(k1,k2) = std(Xc(:,k1) - Xu(:,k2));
        D3(k1,k2) = mean(abs(Xcb(:,k1) - Xub(:,k2)));
        D4(k1,k2) = std(Xcb(:,k1) - Xub(:,k2));
    end
end
constant1 = 10000;
matchArray = zeros(NumBand,NumBand);
kk = 0;

for k1 = 1:NumBand
    for k2 = 1:NumBand
        rul1 = (min(D1(k1,:)) == D1(k1,k2))&& (min(D1(:,k2)) == D1(k1,k2)) && (D1(k1,k2)<constant1); 
        rul2 = (min(D2(k1,:)) == D2(k1,k2)) && (min(D2(:,k2)) == D2(k1,k2));
        rul3 = (min(D3(k1,:)) == D3(k1,k2))&& (min(D3(:,k2)) == D3(k1,k2));
        rul4 = (min(D4(k1,:)) == D4(k1,k2))&& (min(D4(:,k2)) == D4(k1,k2));
        if rul1 && rul2 && rul3 && rul4
            kk = kk + 1;
            matchArray(k1,k2) = kk;
            D1(k1,:) = constant1;
            D1(:,k2) = constant1;
            D2(k1,:) = constant1;
            D2(:,k2) = constant1;
            D3(k1,:) = constant1;
            D3(:,k2) = constant1;
            D4(k1,:) = constant1;
            D4(:,k2) = constant1;
        end
    end
end

for k1 = 1:NumBand
    for k2 = 1:NumBand
        rul1 = (min(D1(k1,:)) == D1(k1,k2))&& (D1(k1,k2)<constant1); 
        rul2 = (min(D2(k1,:)) == D2(k1,k2));
        rul3 = (min(D3(k1,:)) == D3(k1,k2));
        rul4 = (min(D4(k1,:)) == D4(k1,k2));
        if rul1 && rul2 && rul3 && rul4
            kk = kk + 1;
            matchArray(k1,k2) = kk;
            D1(k1,:) = constant1;
            D1(:,k2) = constant1;
            D2(k1,:) = constant1;
            D2(:,k2) = constant1;
            D3(k1,:) = constant1;
            D3(:,k2) = constant1;
            D4(k1,:) = constant1;
            D4(:,k2) = constant1;
        end
    end
end

for k1 = 1:NumBand
    for k2 = 1:NumBand
        rul1 = (min(D1(k1,:)) == D1(k1,k2))&& (D1(k1,k2)<constant1); 
        rul2 = (min(D2(k1,:)) == D2(k1,k2))&& (D1(k1,k2)<constant1);
        rul3 = (min(D3(k1,:)) == D3(k1,k2))&& (D1(k1,k2)<constant1);
        rul4 = (min(D4(k1,:)) == D4(k1,k2))&& (D1(k1,k2)<constant1);
        if rul2 && rul4
            kk = kk + 1;
            matchArray(k1,k2) = kk;
            D1(k1,:) = constant1;
            D1(:,k2) = constant1;
            D2(k1,:) = constant1;
            D2(:,k2) = constant1;
            D3(k1,:) = constant1;
            D3(:,k2) = constant1;
            D4(k1,:) = constant1;
            D4(:,k2) = constant1;
        end
    end
end


% Xcb = zeros(4,NumBand);
% Xub = zeros(4,NumBand);
% 
% Xcb(1,:) = mean(Xc);
% Xcb(2,:) = std(Xc);
% Xcb(3,:) = mean(Xc(1:20,:));
% Xcb(4,:) = mean(Xc(end-19:end,:));
% 
% Xub(1,:) = mean(Xu);
% Xub(2,:) = std(Xu);
% Xub(3,:) = mean(Xu(1:20,:));
% Xub(4,:) = mean(Xu(end-19:end,:));
% % % Xcb = bsxfun(@minus,Xc,mean(Xc));
% % % Xub = bsxfun(@minus,Xu,mean(Xu));
% % % xs = [1, 1,1];
% % % D = 10000*ones(NumBand,NumBand,3);
% % % for k1 = 1:NumBand
% % %     for k2 = 1:NumBand
% % %         for k3 = 1:3
% % %             D(k1,k2,k3) = sum(abs(Xc(:,k1) - xs(k3)*Xu(:,k2)));
% % %         end
% % %     end
% % % end
% % % D2 = min(D,[],3);
% % % 
% % % D = 10000*ones(NumBand,NumBand,3);
% % % for k1 = 1:NumBand
% % %     for k2 = 1:NumBand
% % %         for k3 = 1:3
% % %             D(k1,k2,k3) = sum(abs(xs(k3)*Xc(:,k1) - Xu(:,k2)));
% % %         end
% % %     end
% % % end
% % % D3 = min(D,[],3);

% D = 10000*ones(NumBand,NumBand);
% std_treshC = median(Xcb(2,:));
% std_treshU = median(Xub(2,:));
% temp = abs(Xcb(4,:) - Xcb(3,:));
% diff_threshC = mean(temp);
% temp = abs(Xub(4,:) - Xub(3,:));
% diff_threshU = mean(temp);
% for k1 = 1:NumBand
%     for k2 = 1:NumBand
%         if  Xcb(2,k1) <= std_treshC 
%             if Xub(2,k2) <= std_treshU
%                 D(k1,k2) = abs(Xcb(1,k1) - Xub(1,k2));
%             else
%                 D(k1,k2) = sum(abs(Xcb(:,k1) - Xub(:,k2)));
%             end
%         elseif  abs(Xcb(4,k1)-Xcb(3,k1)) <= diff_threshC 
%             if abs(Xub(4,k1)-Xub(3,k1)) <= diff_threshU 
%                 D(k1,k2) = abs(Xcb(4,k1) - Xub(4,k2)) + abs(Xcb(3,k1) - Xub(3,k2));
%             else
%                 D(k1,k2) = sum(abs(Xc(k1,:) - Xu(k2,:)));
%             end
%         else
%              D(k1,k2) = sum(abs(Xc(k1,:) - Xu(k2,:)));
%         end
%     end
% end
% D = pdist2(Xcb',Xub'); 
% % DfXc = Xcb(2:end,:) - Xcb(1:end-1,:);
% % DfXu = Xub(2:end,:) - Xub(1:end-1,:);
% % 
% % 
% % CorrE_X = zeros(NumBand,NumBand);
% % CorrE_DX = zeros(NumBand,NumBand);
% % temp = zeros(2,2);
% % for k1 = 1:NumBand
% %     for k2 = 1:NumBand
% %         temp = corrcoef(Xc(:,k1),Xu(:,k2));
% %         CorrE_X(k1,k2) = temp(2,1);
% %         
% %         temp = corrcoef(DfXc(:,k1),DfXu(:,k2));
% %         CorrE_DX(k1,k2) = temp(2,1);
% %     end
% % end