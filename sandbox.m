% Prepare top exponential function analysis similar to power function and
% exponential function

row_increment = 0;
desired_b_search = [-5,5];
desired_a_search = [-5,5];
increment = 0.05;

desired_b_search = desired_b_search(1):increment:desired_b_search(2);
desired_a_search = desired_a_search(1):increment:desired_a_search(2);

RMSErecord = zeros(length(desired_b_search)+1,length(desired_a_search)+1);
RMSErecord(2:end,1) = desired_b_search';
RMSErecord(1,2:end) = desired_a_search;

for i = 2:length(desired_b_search)+1
    bval = desired_b_search(i-1);
    disp(bval)
    for j = 2:length(desired_a_search)+1
        aval = desired_a_search(j-1);
        
        balA = zeros(1,100);
        balB = ones(1,100);
        
        for Iter = 1:NIter
            for k = 2:length(Cij)
                sigma = 0;
                for l = 2:length(Cij)
                    sigma = sigma + balB(l-1) * TPTA(l-1,3) * topexpf(k,l,aval,bval,Cij);
                end
                balA(k-1) = 1/sigma;
            end
            
            for k = 2:length(Cij)
                sigma = 0;
                for l = 2:length(Cij)
                    sigma = sigma + balA(l-1) * TPTA(l-1,2) * topexpf(k,l,aval,bval,Cij);
                end
                balB(k-1) = 1/sigma;
            end
        end
        
        balA = balA';
        numObsTij(2:101,102) = balA;
        numObsTij(102,2:101) = balB;
        
        for k = 2:length(Cij)
            for l = 2:length(Cij)
                newTij(k,l) = balA(k-1) * balB(l-1) * TPTA(k-1,2) * TPTA(l-1,3) * topexpf(k,l,aval,bval,Cij);
            end
        end
        
        diffTij = newTij(2:101,2:101) - numObsTij(2:101,2:101);
        diffTij = diffTij.*diffTij;
        RMSE = sqrt(sum(sum(diffTij))/1000000);
        
        RMSErecord(i,j) = RMSE;
        
    end
end

function topexpCost = topexpf(i,j,a,b,Cij)
topexpCost = (Cij(i,j).^a)*exp(-b*Cij(i,j));
end
% RMSErecord(i,j) = rand(1);