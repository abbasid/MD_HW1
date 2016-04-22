function ang_vel_90 = fourbar_Ang_Vel_90( theta_dot_90 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
R_ab = 100;
R_bc = 180;
R_dc = 250;
theta_90 = [90, -33.69, 53.07];
theta_dot1 = theta_dot_90(1);
theta_dot2 = theta_dot_90(2);
ang_vel_90 = [R_ab*cosd(theta_90(1, 1))*120 - R_dc*cosd(theta_90(1, 3))*theta_dot1 ...
                             + R_bc*cosd(theta_90(1, 2))*theta_dot2; ...
              -R_ab*sind(theta_90(1, 1))*120 + R_dc*sind(theta_90(1, 3))*theta_dot1 ...
                            - R_bc*sind(theta_90(1, 2))*theta_dot2];

end

