d1 = load("x_-0.55_y_0_1_case_1/panda_left_x_out.txt");
d2 = load("x_-0.55_y_0_1_case_2/panda_left_x_out.txt");

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