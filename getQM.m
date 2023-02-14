function [Q, M] = getQM(n_seg, n_order, ts)
    Q = [];
    M = [];
    M_k = getM(n_order);
    reduce_order = (n_order + 1)/2;
    for k = 1:n_seg
        %#####################################################
        % STEP 2.1 calculate Q_k of the k-th segment 
        % NOTE: M decide the order of the coeff c
        % c = [c0,c1,c2,..cn];
        Q_k = [];
        Q_k = zeros(n_order+1,n_order+1);
        for i = reduce_order:n_order
            for l = reduce_order:n_order
                Q_k(i+1,l+1) = factorial(i)*factorial(l)/(factorial(i-reduce_order)*factorial(l-reduce_order)*(i+l-reduce_order));
                
            end
        end
        Q_k = Q_k/ts(k)^(2*reduce_order-3);
        Q = blkdiag(Q, Q_k);
        M = blkdiag(M, M_k);
    end
end
