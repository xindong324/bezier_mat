clc;
clear;
close all
v_max = 400;
a_max = 400;
color = ['r', 'b', 'm', 'g', 'k', 'c', 'c'];

%% specify the center points of the flight corridor and the region of corridor
path = [50, 50;
       100, 120;
       180, 150;
       250, 80;
       280, 0];
x_length = 100;
y_length = 100;

path2 = [8.87096774193549	19.9708454810496;
13.9400921658986	55.5393586005831;
35.1382488479263	66.6180758017493;
65.5529953917051	53.7900874635568;
88.3640552995392	42.1282798833819;
94.3548387096774	77.9883381924198];
path2 = path2;

n_order = 7;   % 8 control points
n_seg = size(path2, 1)-1;

% calculate time distribution in proportion to distance between 2 points
dist     = zeros(n_seg, 1);
dist_sum = 0;
T        = 10;
t_sum    = 0;

for i = 1:n_seg
    dist(i) = sqrt((path2(i+1, 1)-path2(i, 1))^2 + (path2(i+1, 2) - path2(i, 2))^2);
    dist_sum = dist_sum+dist(i);
end
for i = 1:n_seg-1
    ts(i) = dist(i)/dist_sum*T;
    t_sum = t_sum+ts(i);
end
ts(n_seg) = T - t_sum;

% %% specify ts for each segment
% ts = zeros(n_seg, 1);
% for i = 1:n_seg
%     ts(i,1) = 1;
% end

poly_coef_x = MinimumSnapCorridorBezierSolver(1, path2(:, 1),  ts, n_seg, n_order, v_max, a_max);
poly_coef_y = MinimumSnapCorridorBezierSolver(2, path2(:, 2),  ts, n_seg, n_order, v_max, a_max);

for i = 0:n_seg-1
    interval = ts(i+1);
    for k = 1:n_order
        vel_coef_x(i*n_order+k) = n_order*(poly_coef_x(i*n_order+k+1) - poly_coef_x(i*n_order+k));
        vel_coef_y(i*n_order+k) = n_order*(poly_coef_y(i*n_order+k+1) - poly_coef_y(i*n_order+k));
    end
    for k = 1:n_order-1
        acc_coef_x(i*n_order+k) = n_order*(n_order-1)*(poly_coef_x(i*n_order+k+2) - 2*poly_coef_x(i*n_order+k+1) + poly_coef_x(i*n_order+k))/interval;
        acc_coef_y(i*n_order+k) = n_order*(n_order-1)*(poly_coef_y(i*n_order+k+2) - 2*poly_coef_y(i*n_order+k+1) + poly_coef_y(i*n_order+k))/interval;
    end
end



figure(1)
%% display the trajectory and cooridor
plot(path2(:,1), path2(:,2), '*r'); hold on;
% for i = 1:n_seg
%     plot_rect([corridor(1,i);corridor(2,i)], corridor(3, i), corridor(4,i));hold on;
% end
% hold on;
x_pos = [];y_pos = [];
idx = 1;

%% #####################################################
% STEP 4: draw bezier curve
tspx=[];
tspy=[];
for k = 1:n_seg
    for t = 0:0.01:1
        x_pos(idx) = 0.0;
        y_pos(idx) = 0.0;
        tmpx = 0.;
        tmpy = 0.;
        
        for i = 0:n_order
            if(idx == 101&&i==n_order)
                aaa=1;
            end
            basis_p = nchoosek(n_order, i) * t^i * (1-t)^(n_order-i);
%             tmpx = tmpx + ts(k)*poly_coef_x((k-1)*(n_order+1)+i+1)*basis_p ;
%             tmpy = tmpy + ts(k)*poly_coef_y((k-1)*(n_order+1)+i+1)*basis_p ;

            x_pos(idx) = x_pos(idx) + ts(k)*poly_coef_x((k-1)*(n_order+1)+i+1)*basis_p ;
            y_pos(idx) = y_pos(idx) + ts(k)*poly_coef_y((k-1)*(n_order+1)+i+1)*basis_p ;
        end
%         x_pos(idx) = tmpx;
%         y_pos(idx) = tmpy;
        % str = ["id: "+num2str(idx)+" x: "+num2str(x_pos(idx))+" y: "+num2str(y_pos(idx))]
        % print("id: %d, x: %lf, y: %lf",idx, tmpx, tmpy);
        idx = idx + 1;
    end

    tspx = [tspx;ts(k)*poly_coef_x((k-1)*(n_order+1)+1:(k)*(n_order+1))];
    tspy = [tspy;ts(k)*poly_coef_y((k-1)*(n_order+1)+1:(k)*(n_order+1))];

end

Mmat = getM(n_order);
invMat = [];
for i = 1:n_seg
    invMat = blkdiag(invMat,Mmat);
end
coeffx = invMat*tspx;
coeffy = invMat*tspy;


% scatter(...);
scatter(ts(1)*poly_coef_x(1:(n_order+1)),ts(1)*poly_coef_y(1:(n_order+1)),100,"r");
scatter(ts(2)*poly_coef_x((n_order+1)+1:2*(n_order+1)),ts(2)*poly_coef_y((n_order+1)+1:2*(n_order+1)),100,"g");
scatter(ts(3)*poly_coef_x(2*(n_order+1)+1:3*(n_order+1)),ts(3)*poly_coef_y(2*(n_order+1)+1:3*(n_order+1)),100,"b");
scatter(ts(4)*poly_coef_x(3*(n_order+1)+1:4*(n_order+1)),ts(4)*poly_coef_y(3*(n_order+1)+1:4*(n_order+1)),100,"c");
scatter(ts(5)*poly_coef_x(4*(n_order+1)+1:5*(n_order+1)),ts(5)*poly_coef_y(4*(n_order+1)+1:5*(n_order+1)),100,"m");
plot(x_pos,y_pos);

figure(2)

hold on
scatter(vel_coef_x(1:(n_order)),vel_coef_y(1:(n_order)),100,"r");
scatter(vel_coef_x((n_order)+1:2*(n_order)),vel_coef_y((n_order)+1:2*(n_order)),100,"g");
scatter(vel_coef_x(2*(n_order)+1:3*(n_order)),vel_coef_y(2*(n_order)+1:3*(n_order)),100,"b");
scatter(vel_coef_x(3*(n_order)+1:4*(n_order)),vel_coef_y(3*(n_order)+1:4*(n_order)),100,"c");
scatter(vel_coef_x(4*(n_order)+1:5*(n_order)),vel_coef_y(4*(n_order)+1:5*(n_order)),100,"m");
idx = 1;
for k = 1:n_seg
    for t = 0:0.01:1
        x_pos(idx) = 0.0;
        y_pos(idx) = 0.0;
        for i = 0:n_order-1
            basis_p = nchoosek(n_order-1, i) * t^i * (1-t)^(n_order-i);
            x_pos(idx) = x_pos(idx) + ts(k)*poly_coef_x((k-1)*(n_order+1)+i+1)*basis_p ;
            y_pos(idx) = y_pos(idx) + ts(k)*poly_coef_y((k-1)*(n_order+1)+i+1)*basis_p ;
        end
        idx = idx + 1;
    end

end

figure(3)
hold on
scatter(acc_coef_x(1:(n_order-1)),ts(1)*acc_coef_y(1:(n_order-1)),100,"r");
scatter(acc_coef_x((n_order-1)+1:2*(n_order-1)),acc_coef_y((n_order-1)+1:2*(n_order-1)),100,"g");
scatter(acc_coef_x(2*(n_order-1)+1:3*(n_order-1)),acc_coef_y(2*(n_order-1)+1:3*(n_order-1)),100,"b");
scatter(acc_coef_x(3*(n_order-1)+1:4*(n_order-1)),acc_coef_y(3*(n_order-1)+1:4*(n_order-1)),100,"c");
scatter(acc_coef_x(4*(n_order-1)+1:5*(n_order-1)),acc_coef_y(4*(n_order-1)+1:5*(n_order-1)),100,"m");

function poly_coef = MinimumSnapCorridorBezierSolver(axis, waypoints, ts, n_seg, n_order, v_max, a_max)
    start_cond = [waypoints(1), 0, 0,0];
    end_cond   = [waypoints(end), 0, 0,0];   
    middle_wp = waypoints(2:end-1);
    
    
    %% #####################################################
    % STEP 1: compute Q_0 of c'Q_0c
    %[Q, M]  = getQM(n_seg, n_order, ts);
    %Q_0 = M'*Q*M;
    Q_0_tmp = getQ02(n_seg, n_order,ts);

    Q_0 = nearestSPD(Q_0_tmp);
    
    %% #####################################################
    % STEP 2: get Aeq and beq
    [Aeq, beq] = getAbeq_wp2(n_seg, n_order, ts, start_cond(1:3), end_cond(1:3),middle_wp);
    
    %% #####################################################
    % STEP 3: get corridor_range, Aieq and bieq 
    
    % STEP 3.1: get corridor_range of x-axis or y-axis,
    % you can define corridor_range as [p1_min, p1_max;
    %                                   p2_min, p2_max;
    %                                   ...,
    %                                   pn_min, pn_max ];

    
    % STEP 3.2: get Aieq and bieq
    %[Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max);
    
    f = zeros(size(Q_0,1),1);
    %poly_coef = quadprog(Q_0,f,Aieq,bieq, Aeq, beq);
    poly_coef = quadprog(Q_0,f,[],[], Aeq, beq);
    %poly_coef = quadprog(Q_0,f, Aeq, beq);
end

function plot_rect(center, x_r, y_r)
    p1 = center+[-x_r;-y_r];
    p2 = center+[-x_r;y_r];
    p3 = center+[x_r;y_r];
    p4 = center+[x_r;-y_r];
    plot_line(p1,p2);
    plot_line(p2,p3);
    plot_line(p3,p4);
    plot_line(p4,p1);
end

function plot_line(p1,p2)
    a = [p1(:),p2(:)];    
    plot(a(1,:),a(2,:),'b');
end

