function ang_vel_270= fourbar_Ang_Vel_270( theta_dot_270 )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
R_ab = 100;
R_bc = 180;
R_dc = 250;
theta_270 = [270 -70.56 16.2];
theta_dot1 = theta_dot_270(1);
theta_dot2 = theta_dot_270(2);

ang_vel_270 = [R_ab*cosd(theta_270(1, 1)) * 120 - R_dc*cosd(theta_270(1, 3)*theta_dot1 ...
                             - R_bc*cosd(theta_270(1, 2))*theta_dot2); ...
                R_ab*-sind(theta_270(1, 1)) * 120 + R_dc*sind(theta_270(1, 3)*theta_dot1 ...
                             - R_bc*sind(theta_270(1, 2))*theta_dot2)];


end

