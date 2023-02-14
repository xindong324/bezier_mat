function Ct = getCt(n_seg, n_order)
    %##########
    %dim decided by reduce order
    %
    state_each_wp = 4;
    state_no_p = state_each_wp-1;
    state_each_seg = 2*state_each_wp;
    dim_d = state_each_seg*n_seg;
    dim_reduced_d = state_each_wp*(n_seg+1);
    dim_reduced_df = state_each_wp*2+n_seg-1;
    dim_reduced_dp = (state_each_wp-1)*(n_seg-1);
    %d = ct*reduced_d;
    Ct = zeros(dim_d,dim_reduced_d);
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    %
    %
    %
    % df = [start_pvaj mid_p end_pvaj]'
    % dp = [mid_vaj]'
    % d = [start_pvaj end_pvaj] of each seg;
    %start
    Ct(1,1) =  1;   %p
    Ct(2,2) = 1;    %v
    Ct(3,3) = 1;    %a
    Ct(4,4) = 1;    %j
    
    %end
    Ct(dim_d-3,dim_reduced_df-3) = 1;
    Ct(dim_d-2,dim_reduced_df-2) = 1;
    Ct(dim_d-1,dim_reduced_df-1) = 1;
    Ct(dim_d,dim_reduced_df) = 1;
    
    for i = 1:n_seg-1
        %each end p of seg and start p of next seg
        %end wp
        Ct(state_each_seg*(i-1)+state_each_wp+1,state_each_wp+i) = 1;
        % next start wp
        Ct(state_each_seg*(i)+1,state_each_wp+i) = 1;
        
        %continue of vaj
        %end wv
        Ct(state_each_seg*(i-1)+state_each_wp+2,dim_reduced_df+state_no_p*(i-1)+1) = 1;
        % next start wv
        Ct(state_each_seg*(i)+2,dim_reduced_df+state_no_p*(i-1)+1) = 1;
        
        %end wa
        Ct(state_each_seg*(i-1)+state_each_wp+3,dim_reduced_df+state_no_p*(i-1)+2) = 1;
        % next start wa
        Ct(state_each_seg*(i)+3,dim_reduced_df+state_no_p*(i-1)+2) = 1;
        
        %end wj
        Ct(state_each_seg*(i-1)+state_each_wp+4,dim_reduced_df+state_no_p*(i-1)+3) = 1;
        % next start wj
        Ct(state_each_seg*(i)+4,dim_reduced_df+state_no_p*(i-1)+3) = 1;
        
    end
    
end