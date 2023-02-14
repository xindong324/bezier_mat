function Q = getQ(n_seg, n_order, ts)
    Q = [];
    reduce_order = 4;
    for k = 1:n_seg
        Q_k = [];
        Q_k = zeros(n_order+1,n_order+1);
        %#####################################################
        % STEP 1.1: calculate Q_k of the k-th segment 
        %
        %
        %
        %
        for i = reduce_order:n_order
            for l = reduce_order:n_order
                Q_k(n_order-i+1,n_order-l+1) = factorial(i)*factorial(l)*(ts(k)^(i+l-2*reduce_order+1))/(factorial(i-reduce_order)*factorial(l-reduce_order)*(i+l-2*reduce_order+1));
                
            end
        end
        
            
        Q = blkdiag(Q, Q_k);
    end
end
