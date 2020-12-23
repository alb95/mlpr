function K = turnaround_samples(n, number_to_test)
%useful to generate tensors with solution curve with an "S" shape
K = {};
v = ones(n,1)/n;
j = 1;
    for i = 1:number_to_test
        R = create_tensor(n);
        [x,it,did] = continuation_eulernewton_bis(0.999, v, R);
        if did == 1
            K{j} = R;
            j = j+1;
        end
    end
end    

