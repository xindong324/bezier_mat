function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        M_k = [];
        M_k = zeros(8,n_order+1);
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        %
        %
        %
        %
        %p0
        M_k(1,n_order+1) = 1; %p
        M_k(2,n_order) = 1;  %v
        M_k(3,n_order-1) = 2;   %a;
        M_k(4,n_order-2) = factorial(3); %j;
        %pt;
        for i = 1:n_order+1
            M_k(5,i) = ts(k)^(n_order+1-i);
            if(i<=n_order)
                M_k(6,i) = (n_order+1-i)*ts(k)^(n_order-i);
            end
            if(i<=n_order-1)
                M_k(7,i) = (n_order+1-i)*(n_order-i)*ts(k)^(n_order-i-1);
            end
            if(i<=n_order-2)
                M_k(8,i) = (n_order+1-i)*(n_order-i)*(n_order-i-1)*ts(k)^(n_order-i-2);
            end
        end
            
        M = blkdiag(M, M_k);
    end
end