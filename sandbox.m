row_increment = 0;
desired_b_search = [-3,3];
desired_a_search = [-3,3];
increment = 0.01;

desired_b_search = desired_b_search(1):increment:desired_b_search(2);
desired_a_search = desired_a_search(1):increment:desired_a_search(2);

RMSErecord = zeros(length(desired_b_search)+1,length(desired_a_search)+1);
RMSErecord(2:end,1) = desired_b_search';
RMSErecord(1,2:end) = desired_a_search;

for i = 2:length(desired_b_search)+1
    for j = 2:length(desired_a_search)+1
        RMSErecord(i,j) = rand(1);
    end
end


