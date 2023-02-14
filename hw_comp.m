clc;clear;close all;
%path = ginput() * 100.0;
path = [8.87096774193549	19.9708454810496;
13.9400921658986	55.5393586005831;
35.1382488479263	66.6180758017493;
65.5529953917051	53.7900874635568;
88.3640552995392	42.1282798833819;
94.3548387096774	77.9883381924198];

n_order = 7;
n_seg = size(path, 1) - 1;
n_poly_perseg = n_order + 1;

ts = zeros(n_seg, 1);
% calculate time distribution based on distance between 2 points
dist = zeros(n_seg, 1);
dist_sum = 0;
T = 25;

t_sum = 0;
for i = 1:n_seg
    dist(i) = sqrt((path(i+1, 1) - path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
    dist_sum = dist_sum + dist(i);
end
for i = 1:n_seg-1
    ts(i) = dist(i) / dist_sum * T;
    t_sum = t_sum + ts(i);
end
ts(n_seg) = T - t_sum;
% or you can simply average the time
% for i = 1:n_seg
%     ts(i) = 1.0;
% end

poly_coef_x2 = MinimumSnapCloseformSolver(path(:, 1), ts, n_seg, n_order);
poly_coef_y2 = MinimumSnapCloseformSolver(path(:, 2), ts, n_seg, n_order);

poly_coef_x1 = MinimumSnapQPSolver(path(:, 1), ts, n_seg, n_order);
poly_coef_y1 = MinimumSnapQPSolver(path(:, 2), ts, n_seg, n_order);

X_n1 = [];
Y_n1 = [];
X_n2 = [];
Y_n2 = [];
k = 1;
tstep = 0.01;

poly_coef_x_rev = [];
poly_coef_y_rev = [];

for i=0:n_seg-1
    %#####################################################
    % STEP 4: get the coefficients of i-th segment of both x-axis
    % and y-axis
    Pxi1 = poly_coef_x1(1+i*(n_order+1):(i+1)*(n_order+1));
    Pyi1 = poly_coef_y1(1+i*(n_order+1):(i+1)*(n_order+1));

    coef_tmp_x = Pxi1(end:-1:1);
    coef_tmp_y = Pyi1(end:-1:1);
    poly_coef_x_rev = [poly_coef_x_rev;coef_tmp_x];
    poly_coef_y_rev = [poly_coef_y_rev;coef_tmp_y];
    
    Pxi2 = poly_coef_x2(1+i*(n_order+1):(i+1)*(n_order+1));
    Pyi2 = poly_coef_y2(1+i*(n_order+1):(i+1)*(n_order+1));
    for t=0:tstep:ts(i+1)
        X_n1(k)  = polyval(Pxi1,t);
        Y_n1(k)  = polyval(Pyi1,t);
        X_n2(k)  = polyval(Pxi2,t);
        Y_n2(k)  = polyval(Pyi2,t);
        k = k+1;
    end
end

plot(X_n1, Y_n1 ,'Color',[0 1.0 0],'LineWidth',2);
hold on
plot(X_n2, Y_n2 ,'Color',[1.0 0 0],'LineWidth',1);
hold on
scatter(path(1:size(path,1),1),path(1:size(path,1),2));
