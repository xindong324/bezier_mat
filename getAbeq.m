function [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 2.1 p,v,a constraint in start 
    Aeq_start = zeros(3, n_all_poly);
    beq_start = zeros(3, 1); 

    Aeq_start(1, 1) = 1;
    Aeq_start(2,  1) = -n_order;
    Aeq_start(2,2) = n_order;
    Aeq_start(3,1) = n_order*(n_order - 1)/ts(1);
    Aeq_start(3,2) = -2*n_order*(n_order - 1)/ts(1);
    Aeq_start(3,3) = n_order*(n_order - 1)/ts(1);

    beq_start(1:3) = start_cond';

    
    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    Aeq_end = zeros(3, n_all_poly);
    beq_end = zeros(3, 1);

    Aeq_end(1,n_all_poly) = 1;
    Aeq_end(2,n_all_poly) = n_order;
    Aeq_end(2,n_all_poly-1) = -n_order;
    Aeq_end(3,n_all_poly) = n_order*(n_order - 1)/ts(n_seg);
    Aeq_end(3,n_all_poly-1) = -2*n_order*(n_order - 1)/ts(n_seg);
    Aeq_end(3,n_all_poly-2) =  n_order*(n_order - 1)/ts(n_seg);

    beq_end(1:3) = end_cond';
    
    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);

    for i=1:n_seg-1
        Aeq_con_p(i,(i)*(n_order+1)) = -1;
        Aeq_con_p(i,i*(n_order+1)+1) = 1;
    end

    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);

    for i=1:n_seg-1
        Aeq_con_v(i,(i)*(n_order+1)) = -n_order;
        Aeq_con_v(i,i*(n_order+1) - 1) = n_order;
        Aeq_con_v(i,i*(n_order+1)+1) = -n_order;
        Aeq_con_v(i, i*(n_order+1)+2) = n_order;
    end

    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);

    for i=1:n_seg-1
        Aeq_con_a(i,(i)*(n_order+1)) = -n_order*(n_order-1)/ts(i);
        Aeq_con_a(i,(i)*(n_order+1)-1) = 2*n_order*(n_order-1)/ts(i);
        Aeq_con_a(i,(i)*(n_order+1)-2) = -n_order*(n_order-1)/ts(i);
        Aeq_con_a(i,i*(n_order+1)+1) = n_order*(n_order-1)/ts(i+1);
        Aeq_con_a(i,(i)*(n_order+1)+2) = -2*n_order*(n_order-1)/ts(i+1);
        Aeq_con_a(i,(i)*(n_order+1)+3) = n_order*(n_order-1)/ts(i+1);
    end

    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start; Aeq_end; Aeq_con];
    beq = [beq_start; beq_end; beq_con];
end