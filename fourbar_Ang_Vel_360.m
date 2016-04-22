function ang_vel_360 = fourbar_Ang_Vel_360( theta_dot_360 )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
R_ab = 100;
R_bc = 180;
R_dc = 250;
theta_360 = [0 -25.57 18.1];
theta_dot1 = theta_dot_360(1);
theta_dot2 = theta_dot_360(2);
ang_vel_360 = [R_ab*cosd(theta_360(1, 1)) * 120 - R_dc*cosd(theta_360(1, 3))*theta_dot1 ...
                             + R_bc*cosd(theta_360(1, 2))*theta_dot2; ...
               -R_ab*sind(theta_360(1, 1)) * 120 + R_dc*sind(theta_360(1, 3))*theta_dot1 ...
                             - R_bc*sind(theta_360(1, 2))*theta_dot2];


end

