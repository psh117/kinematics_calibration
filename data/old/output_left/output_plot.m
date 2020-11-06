d1 = load("x_out_1.txt");
d2 = load("x_out_2.txt");

close all
d1_true = [-0.55, 0.0, 0.064];
d2_true = [0.0, -0.45, 0.064];
figure 
hold on
for i=1:3
    plot (d1(:,i) - d1_true(i))
    plot (d2(:,i) - d2_true(i))
end
legend('x1','x2','y1','y2','z1','z2')
