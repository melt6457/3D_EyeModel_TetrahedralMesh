C = cell(3,3,1);
celldisp(C);
for i=1:3
    for j=1:3
        C(i,j) = {floor(10*rand(1,3))};
    end
end

C = C';
Val = cellfun(@sum, C)

%numericVec = cell2mat(C);
