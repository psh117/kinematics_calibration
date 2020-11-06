d1 = load("x_67.5_y_0_case1/panda_right_x_out.txt");
d2 = load("x_67.5_y_0_case2/panda_right_x_out.txt");

close all

figure
plot (d1(:,1))
hold on
plot (d2(:,1))

figure
plot (d1(:,2))
hold on
plot (d2(:,2))

figure
plot (d1(:,3))
hold on
plot (d2(:,3))