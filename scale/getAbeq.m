function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1); 
    % STEP 2.1: write expression of Aeq_start and beq_start
    %
    %
    %
    %
    %pos
    Aeq_start(1,n_order+1) = 1;
    Aeq_start(2,n_order) = factorial(1)/ts(1);
    Aeq_start(3,n_order-1) = factorial(2)/ts(1)^2;
    Aeq_start(4,n_order-2) = factorial(3)/ts(1)^3;
    
    beq_start = start_cond';
     
    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    %
    %
    %

    for i = 1:n_order+1
        Aeq_end(1,i+(n_seg-1)*(n_order+1)) = 1;
        if(i<=n_order)
            Aeq_end(2,i+(n_seg-1)*(n_order+1)) = (n_order+1-i)/ts(n_seg);
        end
        if(i<=n_order-1)
            Aeq_end(3,i+(n_seg-1)*(n_order+1)) = (n_order+1-i)*(n_order-i)/ts(n_seg)^2;
        end
        if(i<=n_order-2)
            Aeq_end(4,i+(n_seg-1)*(n_order+1)) = (n_order+1-i)*(n_order-i)*(n_order-i-1)/ts(n_seg)^3;
        end
    end
    
    beq_end = end_cond';
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    %
    %
    %
    for i = 1:n_seg-1
        % tail of seg i
        for j = 1:n_order+1
            Aeq_wp(i,j+(i-1)*(n_order+1)) = 1;
        end
        % header of seg i+1
       % Aeq_wp(i+1,(n_order+1)*i) = 1;
    end
    beq_wp = waypoints(2:n_seg,1);
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    %
    %
    %
    for i = 1:n_seg-1
        % tail of seg i
        for j = 1:n_order+1
            Aeq_con_p(i,j+(i-1)*(n_order+1)) = 1;
        end
        % header of seg i+1
        Aeq_con_p(i,(n_order+1)*(i+1)) = -1;
    end
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %
    for i = 1:n_seg-1
        % tail of seg i
        for j = 1:n_order
            Aeq_con_v(i,j+(i-1)*(n_order+1)) = (n_order+1-j)/ts(i);
        end
        % header of seg i+1
        Aeq_con_v(i,(n_order+1)*(i+1)-1) = -1/ts(i+1);
    end
    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %
    for i = 1:n_seg-1
        % tail of seg i
        for j = 1:n_order-1
            Aeq_con_a(i,j+(i-1)*(n_order+1)) = (n_order+1-j)*(n_order-j)/ts(i)^(2);
        end
        % header of seg i+1
        Aeq_con_a(i,(n_order+1)*(i+1)-2) = -2/ts(i+1)^2;
    end
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    for i = 1:n_seg-1
        % tail of seg i
        for j = 1:n_order-2
            Aeq_con_j(i,j+(i-1)*(n_order+1)) = (n_order+1-j)*(n_order-j)*(n_order-j-1)*ts(i)^(3);
        end
        % header of seg i+1
        Aeq_con_j(i,(n_order+1)*(i+1)-3) = -6/ts(i+1)^3;
    end
     
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end