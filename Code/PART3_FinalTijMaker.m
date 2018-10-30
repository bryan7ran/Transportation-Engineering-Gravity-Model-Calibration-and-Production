powerb = 0.373;
expb = 0.2665;
topexpa = 0.154;
topexpb = 0.354;

newTij(2:end,2:end) = [0];
finalTij_power = newTij;
finalTij_exp = finalTij_power;
finalTij_topexp = finalTij_power;

balA_power = zeros(1,100);
balA_exp = zeros(1,100);
balA_topexp = zeros(1,100);
balB_power = ones(1,100);
balB_exp = ones(1,100);
balB_topexp = ones(1,100);

for Iter = 1:1000 % OR NIter
    for i = 2:length(finalTij_power)
        sigma_power = 0;
        sigma_exp = 0;
        sigma_topexp = 0;
        for j = 2:length(finalTij_power)
            sigma_power = sigma_power + balB_power(j-1) * TPTA(j-1,3) * powerf(i,j,powerb,Cij);
            sigma_exp = sigma_exp + balB_exp(j-1) * TPTA(j-1,3) * expf(i,j,expb,Cij);
            sigma_topexp = sigma_topexp + balB_topexp(j-1) * TPTA(j-1,3) * topexpf(i,j,topexpa,topexpb,Cij);
        end
        balA_power(i-1) = 1/sigma_power;
        balA_exp(i-1) = 1/sigma_exp;
        balA_topexp(i-1) = 1/sigma_topexp;
        disp(balA_power(4))
    end
           
    for i = 2:length(finalTij_power)
        sigma_power = 0;
        sigma_exp = 0;
        sigma_topexp = 0;
        for j = 2:length(finalTij_power)
            sigma_power = sigma_power + balA_power(j-1) * TPTA(j-1,2) * powerf(i,j,powerb,Cij);
            sigma_exp = sigma_exp + balA_exp(j-1) * TPTA(j-1,2) * expf(i,j,expb,Cij);
            sigma_topexp = sigma_topexp + balA_topexp(j-1) * TPTA(j-1,2) * topexpf(i,j,topexpa,topexpb,Cij);
        end
        balB_power(i-1) = 1/sigma_power;
        balB_exp(i-1) = 1/sigma_exp;
        balB_topexp(i-1) = 1/sigma_topexp;
    end
end

for i = 2:length(finalTij_power)
    for j = 2:length(finalTij_power)
        finalTij_power(i,j) = balA_power(i-1) * balB_power(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * powerf(i,j,powerb,Cij);
        finalTij_exp(i,j) = balA_exp(i-1) * balB_exp(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * expf(i,j,expb,Cij);
        finalTij_topexp(i,j) = balA_topexp(i-1) * balB_topexp(j-1) * TPTA(i-1,2) * TPTA(j-1,3) * topexpf(i,j,topexpa,topexpb,Cij);
    end
end

powerTPTA = numParking;
expTPTA = numParking;
topexpTPTA = numParking;
% For the length of this new variable...
for i = 1:length(powerTPTA)
    % Define trip productions by each zone as the sum of each column in a
    % row of the Tij matrix
    powerTPTA(i,2) = sum(finalTij_power(i+1,2:end)); % Oi
    expTPTA(i,2) = sum(finalTij_exp(i+1,2:end));
    topexpTPTA(i,2) = sum(finalTij_topexp(i+1,2:end));
    % Similar for trip attractions
    powerTPTA(i,3) = sum(finalTij_power(2:end,i+1)); % Dj
    expTPTA(i,3) = sum(finalTij_exp(2:end,i+1)); % Dj
    topexpTPTA(i,3) = sum(finalTij_topexp(2:end,i+1)); % Dj
end


function powerCost = powerf(i,j,b,Cij)
powerCost = Cij(i,j)^-b;
end

function expCost = expf(i,j,b,Cij)
expCost = exp(-b*Cij(i,j));
end

function topexpCost = topexpf(i,j,a,b,Cij)
topexpCost = (Cij(i,j)^a)*exp(-b*Cij(i,j));
end