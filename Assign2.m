clear
clc

% [numObsTij,~,~] = xlsread('Assign2Data.xlsx','ObsTij');
% [numTT,~,~] = xlsread('Assign2Data.xlsx','AutoTTmin');
% [numDist,~,~] = xlsread('Assign2Data.xlsx','AutoDistkm');
% [numParking,~,~] = xlsread('Assign2Data.xlsx','Parking');
% TAZ_header = numObsTij(1,:);
% TAZ_col = numObsTij(:,1);

load('Assign2Data')

disp('import done')

% bValStart = input('b value?');
% aValStart = input('a value?');

VOT = 15; % minimum wage
VOT = VOT/60; % minimum wage per minute
valueVeh = 0.55; % values per kilometre

costTT = numTT * VOT;
costTT(1,:) = TAZ_header;
costTT(:,1) = TAZ_col;

costVeh = numDist * valueVeh;
costVeh(1,:) = TAZ_header;
costVeh(:,1) = TAZ_col;

rpowerF = zeros(1,2);
rexpF = rpowerF;
rtopF = rexpF;

Cij = zeros(size(costTT));
Cij(1,:) = TAZ_header;
Cij(:,1) = TAZ_col;

NIter = 20;

% matrix maker here
newTij = Cij;

for i = 2:size(Cij,1)
    TAZDest = Cij(1,i);
    parkingIndex = find(numParking == TAZDest);
    costPark = numParking(parkingIndex,2);
    
    Cij(2:end,i) = costPark;
    % now the cost table is full of parking costs
end

% populate the travel cost table now

Cij = Cij + costTT + costVeh;
Cij(1,:) = TAZ_header;
Cij(:,1) = TAZ_col;

TPTA = numParking;
for i = 1:length(TPTA)
    TPTA(i,2) = sum(numObsTij(i+1,2:end)); % Oi
    TPTA(i,3) = sum(numObsTij(2:end,i+1)); % Dj
    
end
% TPTA = round(TPTA);

% can begin trying to obtain gravity models now to see how off you are
% for master = 1:100001
%     for aValStart = 0:0.1:0.1
%         for bValStart = 0:0.1:100
  row_increment = 0;
  for bValStart = -3:0.01:3
      row_increment = row_increment + 1;
            balA = zeros(1,100);
            balB = ones(1,100);
            
            % Calculate A based on B
            for Iter = 1:NIter
                for i = 2:length(Cij)
                    sigma = 0;
                    for j = 2:length(Cij)
                        sigma = sigma + balB(j-1) * TPTA(j-1,3) * powerf(i,j,bValStart,Cij); %check balA here and balB here BRYAN AND TPTA HOLY
                    end
                    balA(i-1) = 1/sigma;
                end
                
                % Calculate B based on A
                for i = 2:length(Cij)
                    sigma = 0;
                    for j = 2:length(Cij)
                        sigma = sigma + balA(j-1) * TPTA(j-1,2) * powerf(i,j,bValStart,Cij); %check balA here and balB here BRYAN AND TPTA HOLY
                    end
                    balB(i-1) = 1/sigma;
                end
                
            end
            
            % append
            balA = balA';
            numObsTij(2:101,102) = balA;
            numObsTij(102,2:101) = balB;
            
            % Calculate trip interchanges
            
            for i = 2:length(Cij)
                for j = 2:length(Cij)
                    newTij(i,j) = balA(i-1) * balB(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * powerf(i,j,bValStart,Cij);
                end
            end
            
            diffTij = newTij(2:101,2:101) - numObsTij(2:101,2:101);
            diffTij = diffTij.*diffTij;
            RMSE = sqrt(sum(sum(diffTij)));
            
            rpowerF(row_increment,:) = [bValStart,RMSE];
            disp(bValStart)
  end          
%         end
%         disp('sheperds boy')
%     end
%     disp(master)
% end
scatter(rpowerF(:,1),rpowerF(:,2),0.05)

row_increment = 0;
  for bValStart = -3:0.01:3
      row_increment = row_increment + 1;
            balA = zeros(1,100);
            balB = ones(1,100);
            
            % Calculate A based on B
            for Iter = 1:NIter
                for i = 2:length(Cij)
                    sigma = 0;
                    for j = 2:length(Cij)
                        sigma = sigma + balB(j-1) * TPTA(j-1,3) * expf(i,j,bValStart,Cij); %check balA here and balB here BRYAN AND TPTA HOLY
                    end
                    balA(i-1) = 1/sigma;
                end
                
                % Calculate B based on A
                for i = 2:length(Cij)
                    sigma = 0;
                    for j = 2:length(Cij)
                        sigma = sigma + balA(j-1) * TPTA(j-1,2) * expf(i,j,bValStart,Cij); %check balA here and balB here BRYAN AND TPTA HOLY
                    end
                    balB(i-1) = 1/sigma;
                end
                
            end
            
            % append
            balA = balA';
            numObsTij(2:101,102) = balA;
            numObsTij(102,2:101) = balB;
            
            % Calculate trip interchanges
            
            for i = 2:length(Cij)
                for j = 2:length(Cij)
                    newTij(i,j) = balA(i-1) * balB(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * expf(i,j,bValStart,Cij);
                end
            end
            
            diffTij = newTij(2:101,2:101) - numObsTij(2:101,2:101);
            diffTij = diffTij.*diffTij;
            RMSE = sqrt(sum(sum(diffTij)));
            
            rexpF(row_increment,:) = [bValStart,RMSE];
            disp(bValStart)
  end          
scatter(rexpF(:,1),rexpF(:,2),0.05)

function powerCost = powerf(i,j,b,Cij)
powerCost = Cij(i,j)^-b;
end

function expCost = expf(i,j,b,Cij)
expCost = exp(-b*Cij(i,j));
end

function topexpCost = topexpf(i,j,a,b,Cij)
topexpCost = (Cij(i,j)^a)*exp(-b*Cij(i,j));
end



% numObsTij = round(numObsTij);


