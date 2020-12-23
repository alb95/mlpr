function R = create_tensor(n)
%the function generates a stochastic tensor of dimension n x n x n, in flattened matrix form (n x n^2) 
R = zeros(n,n^2);
for i = 1:n
    for j = 1:n^2
        if 0.6 > rand
            R(i,j) = 1;
        end
    end
end
for i = 1:n^2
    if sum(R(:,i))~=0
        R(:,i) = R(:,i)/sum(R(:,i));
    end
    if sum(R(:,i))==0
        row = round(n*rand);
        if row == 0
            row = 1;
        end
        R(row,i)= 1;
    end
end
