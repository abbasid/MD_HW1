function F = fourbar_eqs( theta )

R_1 = 0.1;
R_2 = 0.18;
R_3 = 0.25;
R_4 = 0.3;

theta_1 = theta(1);
theta_2 = theta(2);
theta_3 = theta(3);

F(1) = R_3*cosd(theta_3) + R_2*cosd(theta_2) - R_4 - R_1*cosd(theta_1);
F(2) = R_3*sind(theta_3) + R_2*sind(theta_2) - R_1*sind(theta_1);

end

