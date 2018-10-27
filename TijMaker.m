powerb = 0.01;
expb = 0.01;
topexpa = 0.01;
topexpb = 0.01;

newTij(2:end,2:end) = [0];
finalTij_power = newTij;
finalTij_exp = finalTij_power;
finalTij_topexp = finalTij_power;



function powerCost = powerf(i,j,b,Cij)
powerCost = Cij(i,j)^-b;
end

function expCost = expf(i,j,b,Cij)
expCost = exp(-b*Cij(i,j));
end

function topexpCost = topexpf(i,j,a,b,Cij)
topexpCost = (Cij(i,j)^a)*exp(-b*Cij(i,j));
end