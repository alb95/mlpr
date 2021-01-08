function R = create_unbalanced_tensor(n)
R = zeros(n,n^2);
    for i = 1:n^2
        row = round(n*rand);
        if row == 0
            row = 1;
        end
        R(row,i)= 1;
    end
end
