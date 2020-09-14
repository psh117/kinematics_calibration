d1 = load("x_out_1.txt");
d2 = load("x_out_2.txt");
d3 = load("x_out_3.txt");

close all
d1_true = [-0.65, 0.0, 0.064];
d2_true = [0.0, 0.45, 0.064];
d3_true = [0.675, 0, 0.064];

for i=1:3
    figure
    hold on
    plot (d1(:,i) - d1_true(i))
    plot (d2(:,i) - d2_true(i))
    plot (d3(:,i) - d3_true(i))
    legend('x','y','z')
end
%legend('x1','x2','x3','y1','y2','y3','z1','z2','z3')

% -0.00113419  -0.00316993  0.000600501  -0.00200384  0.000897991  -0.00657346 -2.02488e-12
%in deg:   -0.0649845    -0.181624    0.0344062    -0.114812    0.0514511    -0.376632 -1.16017e-10