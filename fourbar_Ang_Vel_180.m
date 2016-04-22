function ang_vel_180 = fourbar_Ang_Vel_180( theta_dot_180 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
R_ab = 100;
R_bc = 180;
R_dc = 250;
theta_180 = [180 -82.1 45.49];

ang_vel_180 = [R_ab*cosd(theta_180(1, 1)) * 120 - R_dc*cosd(theta_180(1, 3))*theta_dot_180(1) ...
                             - R_bc*cosd(theta_180(1, 2))*theta_dot_180(2); ...
                R_ab*-sind(theta_180(1, 1)) * 120 + R_dc*sind(theta_180(1, 3))*theta_dot_180(1) ...
                             - R_bc*sind(theta_180(1, 2))*theta_dot_180(2)];

end

