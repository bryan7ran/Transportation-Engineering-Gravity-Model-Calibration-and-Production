TPTA = floor(TPTA);
powerTPTA = floor(powerTPTA);
expTPTA = floor(expTPTA);
topexpTPTA = floor(topexpTPTA);

analysisCij = floor(Cij);
% finalTij_power = round(finalTij_power);
% finalTij_exp = round(finalTij_exp);
% finalTij_topexp = round(finalTij_topexp);
obs = round(numObsTij);
table = 0;
table_increment = 0;
for i = 2:length(Cij)
    for j = 2:length(Cij)
        table_increment = table_increment + 1;
        table(table_increment,1) = analysisCij(i,j);
        table(table_increment,2) = finalTij_power(i,j);
        table(table_increment,3) = finalTij_exp(i,j);
        table(table_increment,4) = finalTij_topexp(i,j);
        table(table_increment,5) = obs(i,j);
    end
end
% xlswrite('file456.xlsx',table);
        
        powerComp = abs(TPTA-powerTPTA);
        expComp = abs(TPTA-expTPTA);
        topexpComp = abs(TPTA-topexpTPTA);
        
        powerComp(101,2) = sum(powerComp(:,2));
        powerComp(101,3) = sum(powerComp(:,3));
        
        expComp(101,2) = sum(expComp(:,2));
        expComp(101,3) = sum(expComp(:,3));
        
        topexpComp(101,2) = sum(topexpComp(:,2));
        topexpComp(101,3) = sum(topexpComp(:,3));
