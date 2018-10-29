%% Code Initialization
clear
clc

%% Run this is the .mat file in line 14 is lost
% [numObsTij,~,~] = xlsread('Assign2Data.xlsx','ObsTij');
% [numTT,~,~] = xlsread('Assign2Data.xlsx','AutoTTmin');
% [numDist,~,~] = xlsread('Assign2Data.xlsx','AutoDistkm');
% [numParking,~,~] = xlsread('Assign2Data.xlsx','Parking');
% TAZ_header = numObsTij(1,:);
% TAZ_col = numObsTij(:,1);

% Load the pre-imported data
load('Assign2Data')

disp('import done')

%% Parameters
VOT = 33.82; % Mode value of time in dollars per hour
VOT = VOT/60; % Convert to cost per minute
valueVeh = 0.2245; % Value of vehicle travel time per kilometre

% Make a new cost matrix
costTT = numTT * VOT; % Travel time in minutes multiplied by cost per minute = cost
costTT(1,:) = TAZ_header;
costTT(:,1) = TAZ_col;

% Make another cost matrix
costVeh = numDist * valueVeh; % Distance in kilometres multiplied by cost per km = cost
costVeh(1,:) = TAZ_header;
costVeh(:,1) = TAZ_col;

% Make result matrices (RMSE v. b-val v. a-val)
rpowerF = zeros(1,2);
rexpF = rpowerF;
rtopF = rexpF;

% Make an empty Cij matrix, recall that Cij is the generalized travel cost
% parameter
Cij = zeros(size(costTT));
Cij(1,:) = TAZ_header;
Cij(:,1) = TAZ_col;

% Number of iterations for balancing factor calibration
NIter = 20;

% New trip interchange matrix, just use the template set by Cij to save on
% computational power
newTij = Cij;

% Loop from 2 (avoiding headers) to the length of the cost matrix which
% should be empty at this point
for i = 2:size(Cij,1)
    % Variable TAZDest is equal to the current column number. In an origin
    % destination matrix, the column number of the zonal destination end of
    % a trip
    TAZDest = Cij(1,i);
    % Find the index parking number in the parking costs worksheet
    % provided. If the destination zone that we are interested in has an
    % associated parking cost, this will be the returned index
    parkingIndex = find(numParking == TAZDest);
    % Therefore, using the parking index, we can start filling out our
    % parking cost O/D matrix such that at the destination trip-end trip
    % end of zone X, we are expected to pay its corresponding parking rate
    costPark = numParking(parkingIndex,2);
    % Absorb this matrix's cost into the Cij matrix to save computational
    % power
    Cij(2:end,i) = costPark;
    % now the cost table is full of parking costs
end

% We can now populate our travel cost table, which is the sum of parking
% cost (already in Cij), the travel time cost in dollars, and the vehicle
% cost in dollars
Cij = Cij + costTT + costVeh;

% Make sure the zone names stay the same since the addition in Line 75 will
% have messed it up. Restore the original O/D names.
Cij(1,:) = TAZ_header;
Cij(:,1) = TAZ_col;

% Borrow the numParking variable since it already has the zone numbers laid
% out.
TPTA = numParking;

% For the length of this new variable...
for i = 1:length(TPTA)
    % Define trip productions by each zone as the sum of each column in a
    % row of the original Tij matrix
    TPTA(i,2) = sum(numObsTij(i+1,2:end)); % Oi
    % Similar for trip attractions
    TPTA(i,3) = sum(numObsTij(2:end,i+1)); % Dj
end

% Start an incrementing value such that we can begin to populate an RMSE v.
% b-value results matrix
row_increment = 0;
% Test left and right intervals for b-value here (say 3)
for bValStart = 0.33:0.001:0.4
    % With every iteration on the for loop, the row record will increment
    % by 1
    row_increment = row_increment + 1;
    % Balacing factors pre-allocation for 100 zones (B assume 1 to start)
    balA = zeros(1,100);
    balB = ones(1,100);
    
    % Calculate A based on B for NIter iterations
    for Iter = 1:NIter
        % Calculate for each row
        for i = 2:length(Cij)
            sigma = 0;
            for j = 2:length(Cij)
                sigma = sigma + balB(j-1) * TPTA(j-1,3) * powerf(i,j,bValStart,Cij);
            end
            balA(i-1) = 1/sigma;
        end
        
        % Calculate B based on A
        for i = 2:length(Cij)
            sigma = 0;
            for j = 2:length(Cij)
                sigma = sigma + balA(j-1) * TPTA(j-1,2) * powerf(i,j,bValStart,Cij);
            end
            balB(i-1) = 1/sigma;
        end
        
    end
    
    % Transpose the A balancing factor so that it's compatible with Tij
    % integration
    balA = balA';
    
    % Append to the end of the Tij matrix
    numObsTij(2:101,102) = balA;
    numObsTij(102,2:101) = balB;
    
    clear i
    clear j
    
    % Calculate trip interchanges
    for i = 2:length(Cij)
        for j = 2:length(Cij)
            newTij(i,j) = balA(i-1) * balB(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * powerf(i,j,bValStart,Cij);
        end
    end
    
    % Create a matrix that determines the difference between estimated and
    % observed
    diffTij = newTij(2:101,2:101) - numObsTij(2:101,2:101);
    % Perform array multiplication to square the value
    diffTij = diffTij.*diffTij;
    % Determine RMSE
    RMSE = sqrt(sum(sum(diffTij))/1000000);
    
    % Record the RMSE for each b-value
    rpowerF(row_increment,:) = [bValStart,RMSE];
    disp(bValStart)
end
%         end
%         disp('sheperds boy')
%     end
%     disp(master)
% end

% Plot command readied
scatter(rpowerF(:,1),rpowerF(:,2),0.05)

row_increment = 0;
for bValStart = 0.25:0.001:0.3
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
    RMSE = sqrt(sum(sum(diffTij))/1000000);
    
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


