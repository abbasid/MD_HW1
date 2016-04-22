% calculate four-bar angular velocity 

addpath('D:\NTU_Graduate_school\Semester_104-2\Machine Dynamics\HW1\');
theta = [90 -33.69 53.07; ...
    180 -82.1 45.49; ...
    270 -70.56 16.2; ...
    0 -25.57 18.1]; % theta_1, theta_2, theta_3
%% calculate 0 ~ 360 four bar angle
theta_test = zeros(720, 3);
theta_test(1, :) = [90 -33.69 53.07];

for i = 2: 720
    theta_test(i, :) = fsolve(@fourbar_eqs, [theta_test(i-1, 1)+1 theta_test(i-1, 2) theta_test(i-1, 3)]);
    if theta_test(i, 1) > 359 + 90
        break;
    end
end
theta_test2 = theta_test((1:457), :);

theta = deg2rad(theta);
R1 = 0.1;
R2 = 0.18;
R3 = 0.25;
Rg3b = 0.43;
Rg3c = 0.250;
R43 = 0.125;
R41 = 0.125;
Rbe = 0.36;

%% calculate angular velocity
theta_2_3_dot = zeros(4, 2);
d_theta_23 = zeros(size(theta_test2, 1), 2);
theta_test2 = deg2rad(theta_test2);

% four angles (0 90 180 270) 
for i = 1: 4
   theta_2_3_dot(i, :) =  [-R3*sin(theta(i, 3)), -R2*sin(theta(i, 2));...
                          R3*cos(theta(i, 3)), R2*cos(theta(i, 2))] \ [-R1*sin(theta(i, 1))*2*2*pi, R1*cos(theta(i, 1))*2*2*pi]'; %inv(A)*B
end

% angle from 0 ~ 360 
for i = 1 : size(theta_test2, 1)
   d_theta_23(i, :) = [-R3*sin(theta_test2(i, 3)), -R2*sin(theta_test2(i, 2));...
                          R3*cos(theta_test2(i, 3)), R2*cos(theta_test2(i, 2))] \ [-R1*sin(theta_test2(i, 1))*2*2*pi, R1*cos(theta_test2(i, 1))*2*2*pi]'; 
end

%% calculate angular acceleration
theta_2_3_double_dot = zeros(4, 2);
theta_2_3_dot_squar = [[4*pi 4*pi 4*pi 4*pi]', theta_2_3_dot].^2;
dd_theta_23 = zeros(size(theta_test2, 1), 2);
d_theta_1 = ones(size(theta_test2, 1), 1)*4*pi;
d_theta_23_squar = [d_theta_1, d_theta_23].^2;

syms a b c;
syms a_2 b_2 c_2;


B = [R3*cos(c)*c_2 + R2*cos(b)*b_2 - R1*cos(a)*a_2, R3*sin(c)*c_2 + R2*sin(b)*b_2 - R1*sin(a)*a_2]';
% four angle angular acceleration( 0 90 180 270 )
for i = 1: 4
    A = [-R2*cos(theta(i, 2)), -R3*sin(theta(i, 3)); R2*cos(theta(i, 2)), R3 * cos(theta(i, 3))];
    real_B = eval(subs(B, [a, b, c, a_2, b_2, c_2], [theta(i, 1), theta(i, 2), theta(i, 3), theta_2_3_dot_squar(i, 1), theta_2_3_dot_squar(i, 2), theta_2_3_dot_squar(i, 3)]));
    theta_2_3_double_dot(i, :) = A \ real_B;
end

% angle from 0 ~ 360
for i = 1: size(theta_test2, 1)
     A = [-R2*cos(theta_test2(i, 2)), -R3*sin(theta_test2(i, 3)); R2*cos(theta_test2(i, 2)), R3 * cos(theta_test2(i, 3))];
    real_B = eval(subs(B, [a, b, c, a_2, b_2, c_2], [theta_test2(i, 1), theta_test2(i, 2), theta_test2(i, 3), d_theta_23_squar(i, 1), d_theta_23_squar(i, 2), d_theta_23_squar(i, 3)]));
    dd_theta_23(i, :) = A \ real_B;
end

%% calculate static balancing force

unknown_forces = zeros(9, 4); %F_12x, F_12y, F_23x, F_23y, F_34x, F_34y, F_41x, F_41y, T_12
unknown_forces360 = zeros(9, size(theta_test2, 1));

syms R_2y R_2x R_g3by R_g3bx R_g3cy R_g3cx R_4y R_4x 

coefficient = [1 0 1 0 0 0 0 0 0; ...
                0 1 0 1 0 0 0 0 0; ...
                0 0 R_2y R_2x 0 0 0 0 1; ...
                0 0 1 0 1 0 0 0 0; ...
                0 0 0 1 0 1 0 0 0; ...
                0 0 R_g3by R_g3bx R_g3cy R_g3cx 0 0 0; ...
                0 0 0 0 1 0 1 0 0; ...
                0 0 0 0 0 1 0 1 0; ...
                0 0 0 0 R_4y R_4x 0 0 0];
            
right_hand_side = [0 -1*9.8 0 0 -2*9.8 0 0 -0.2*9.8 -0.2*9.8*R_4x]';

% four angle static force
for i = 1: 4
   real_coefficient = eval(subs(coefficient, [R_2y, R_2x, R_g3by, R_g3bx, R_g3cy, R_g3cx, R_4y, R_4x], [R2*sind(theta(i, 1)), R2*cosd(theta(i, 1)), Rg3b*sind(theta(i, 2)), Rg3b*cosd(theta(i, 2)), ...
       Rg3c*sind(theta(i, 2)), Rg3c*cosd(theta(i, 2)), R3*sind(theta(i, 3)), R3*cosd(theta(i, 3))]));
    real_right_hand_side = eval(subs(right_hand_side, R_4x, R3*cosd(theta(i, 3))));
   unknown_forces(:, i) = real_coefficient \ real_right_hand_side;
end

% static four from 0 ~ 360
for i = 1 : size(theta_test2, 1)
     real_coefficient = eval(subs(coefficient, [R_2y, R_2x, R_g3by, R_g3bx, R_g3cy, R_g3cx, R_4y, R_4x], [R2*sind(theta_test2(i, 1)), R2*cosd(theta_test2(i, 1)), Rg3b*sind(theta_test2(i, 2)), Rg3b*cosd(theta_test2(i, 2)), ...
       Rg3c*sind(theta_test2(i, 2)), Rg3c*cosd(theta_test2(i, 2)), R3*sind(theta_test2(i, 3)), R3*cosd(theta_test2(i, 3))]));
    real_right_hand_side = eval(subs(right_hand_side, R_4x, R3*cosd(theta_test2(i, 3))));
   unknown_forces360(:, i) = real_coefficient \ real_right_hand_side;
end


%% calculate acceleration
acc = zeros(3, 2, 4); % link x (acceleration, angle) x R1 angle 
acc360 = zeros(3, 2, size(theta_test2));
omega = [[4*pi, 4*pi, 4*pi, 4*pi]', theta_2_3_dot];
omega360 = [d_theta_1, d_theta_23];
alpha = [[0, 0, 0, 0]', theta_2_3_double_dot];
alpha360 = [zeros(size(theta_test2, 1), 1), dd_theta_23];

for i = 1: 4
    acc_r1_x = R1*omega(i, 1).^2 * cos(theta(i, 1));
    acc_r1_y = R1*omega(i, 1).^2 * sin(theta(i, 1));
    acc(1, 1, i) = 0;
    acc(1, 2, i) = 0;
    
    acc_r2_2r1_x = Rbe*((omega(i, 2).^2) * cos(theta(i, 2)) + alpha(i, 2) * sin(theta(i, 2)));
    acc_r2_2r1_y = Rbe*((omega(i, 2).^2) * sin(theta(i, 3)) + alpha(i, 2) * cos(theta(i, 2)));
    acc_B_x = acc_r2_2r1_x + acc_r1_x;
    acc_B_y = acc_r2_2r1_y + acc_r1_y;
    acc(2, 1, i) = sqrt(acc_B_x.^2 + acc_B_y.^2);
    acc(2, 2, i) = atan2(acc_B_y, acc_B_x);
    
    acc_r3_x = (R3/2)*(omega(i, 3).^2 * cos(theta(i, 3)) + alpha(i, 3) * sin(theta(i, 3)));
    acc_r3_y = (R3/2)*(omega(i, 3).^2 * sin(theta(i, 3)) + alpha(i, 3) * cos(theta(i, 3)));
    acc(3, 1, i) = sqrt(acc_r3_x.^2 + acc_r3_y.^2);
    acc(3, 2, i) = atan2(acc_r3_y, acc_r3_x);

end

for i = 1 : size(theta_test2, 1)
    acc_r1_x = R1*omega360(i, 1).^2 * cos(theta_test2(i, 1));
    acc_r1_y = R1*omega360(i, 1).^2 * sin(theta_test2(i, 1));
    acc360(1, 1, i) = 0;
    acc360(1, 2, i) = 0;
    
    acc_r2_2r1_x = Rbe*((omega360(i, 2).^2) * cos(theta_test2(i, 2)) + alpha360(i, 2) * sin(theta_test2(i, 2)));
    acc_r2_2r1_y = Rbe*((omega360(i, 2).^2) * sin(theta_test2(i, 3)) + alpha360(i, 2) * cos(theta_test2(i, 2)));
    acc_B_x = acc_r2_2r1_x + acc_r1_x;
    acc_B_y = acc_r2_2r1_y + acc_r1_y;
    acc360(2, 1, i) = sqrt(acc_B_x.^2 + acc_B_y.^2);
    acc360(2, 2, i) = atan2(acc_B_y, acc_B_x);
    
    acc_r3_x = (R3/2)*(omega360(i, 3).^2 * cos(theta_test2(i, 3)) + alpha360(i, 3) * sin(theta_test2(i, 3)));
    acc_r3_y = (R3/2)*(omega360(i, 3).^2 * sin(theta_test2(i, 3)) + alpha360(i, 3) * cos(theta_test2(i, 3)));
    acc360(3, 1, i) = sqrt(acc_r3_x.^2 + acc_r3_y.^2);
    acc360(3, 2, i) = atan2(acc_r3_y, acc_r3_x);
end

%% Calculate D'Alambert forces
I2 = 0.02;
I3 = 0.06;
I4 = 0.005;
m2 = 1;
m3 = 2;
m4 = 0.2;


unknown_dynamic_forces = zeros(9, 4); %F_12x, F_12y, F_23x, F_23y, F_34x, F_34y, F_41x, F_41y, T_12
unknown_dynamic_forces360 = zeros(9, size(theta_test2, 1));
syms Igg2alpha2 M3Ag3x M3Ag3y Igg3Alpha3 M4Ag4x M4Ag4y Igg4Alpha4;
right_hand_side = [0 -1*9.8 Igg2alpha2 M3Ag3x M3Ag3y-2*9.8 Igg3Alpha3 M4Ag4x M4Ag4y-0.2*9.8 Igg4Alpha4-0.2*9.8*R_4x]';

for i = 1 : size(unknown_dynamic_forces, 2)
    real_right_hand_side = eval(subs(right_hand_side, [Igg2alpha2 M3Ag3x M3Ag3y Igg3Alpha3 M4Ag4x M4Ag4y Igg4Alpha4 R_4x], ...
                                [I2*alpha(i, 1) m2*acc(2, 1, i)*cos(acc(1, 2, i)) m2*acc(2, 1, i)*sin(acc(1, 2, i)) ...
                                 I3*alpha(i, 2) m4*acc(3, 1, i)*cos(acc(3, 2, i)) m4*acc(3, 1, i)*sin(acc(3, 2, i)) (I4 + m4*R41.^2) * alpha(i, 3) R3*cosd(theta(i, 3))]));
     real_coefficient = eval(subs(coefficient, [R_2y, R_2x, R_g3by, R_g3bx, R_g3cy, R_g3cx, R_4y, R_4x], [R2*sind(theta(i, 1)), R2*cosd(theta(i, 1)), Rg3b*sind(theta(i, 2)), Rg3b*cosd(theta(i, 2)), ...
       Rg3c*sind(theta(i, 2)), Rg3c*cosd(theta(i, 2)), R3*sind(theta(i, 3)), R3*cosd(theta(i, 3))]));  
   unknown_dynamic_forces(:, i) = real_coefficient \ real_right_hand_side;
end

for i = 1 : size(unknown_dynamic_forces360, 2)
    real_right_hand_side = eval(subs(right_hand_side, [Igg2alpha2 M3Ag3x M3Ag3y Igg3Alpha3 M4Ag4x M4Ag4y Igg4Alpha4 R_4x], ...
                                [I2*alpha360(i, 1) m2*acc360(2, 1, i)*cos(acc360(1, 2, i)) m2*acc360(2, 1, i)*sin(acc360(1, 2, i)) ...
                                 I3*alpha360(i, 2) m4*acc360(3, 1, i)*cos(acc360(3, 2, i)) m4*acc360(3, 1, i)*sin(acc360(3, 2, i)) (I4 + m4*R41.^2) * alpha360(i, 3) R3*cosd(theta_test2(i, 3))]));
     real_coefficient = eval(subs(coefficient, [R_2y, R_2x, R_g3by, R_g3bx, R_g3cy, R_g3cx, R_4y, R_4x], [R2*sin(theta_test2(i, 1)), R2*cos(theta_test2(i, 1)), Rg3b*sin(theta_test2(i, 2)), Rg3b*cos(theta_test2(i, 2)), ...
       Rg3c*sin(theta_test2(i, 2)), Rg3c*cos(theta_test2(i, 2)), R3*sin(theta_test2(i, 3)), R3*cos(theta_test2(i, 3))]));  
   unknown_dynamic_forces360(:, i) = real_coefficient \ real_right_hand_side;
end

%% Calculate shaking forces and shaking moments
shaking_forces = zeros(4, 4); % four R1 angle x (x component,  y component, magnitude, phase)
shaking_forces360 = zeros(size(theta_test2, 1), 4);
shaking_moments = zeros(4, 1);
shaking_moments360 = zeros(size(theta_test2, 1), 4); 

for i = 1: size(shaking_forces, 1)
   shaking_forces(i, 1) = -(unknown_dynamic_forces(1, i) + unknown_dynamic_forces(7, i)); 
   shaking_forces(i, 2) = -(unknown_dynamic_forces(2, i) + unknown_dynamic_forces(8, i));
   shaking_forces(i, 3) = sqrt(shaking_forces(i, 1).^2 + shaking_forces(i, 2).^2);
   shaking_forces(i, 4) = atan2(shaking_forces(i, 2), shaking_forces(i, 1));
   shaking_moments(i, 1) = -(unknown_dynamic_forces(9, i) + unknown_dynamic_forces(5, i)*R3*sin(theta(i, 3)) + unknown_dynamic_forces(6, i)* R3*cos(theta(i, 3))+ m4*9.8*R41*cos(theta(i, 3)));
end

for i = 1: size(shaking_forces360, 1)
   shaking_forces360(i, 1) = -(unknown_dynamic_forces360(1, i) + unknown_dynamic_forces360(7, i)); 
   shaking_forces360(i, 2) = -(unknown_dynamic_forces360(2, i) + unknown_dynamic_forces360(8, i));
   shaking_forces360(i, 3) = sqrt(shaking_forces360(i, 1).^2 + shaking_forces360(i, 2).^2);
   shaking_forces360(i, 4) = atan2(shaking_forces360(i, 2), shaking_forces360(i, 1));
   shaking_moments360(i, 1) = -(unknown_dynamic_forces360(9, i) + unknown_dynamic_forces360(5, i)*R3*sin(theta_test2(i, 3)) + unknown_dynamic_forces360(6, i)* R3*cos(theta_test2(i, 3))+ m4*9.8*R41*cos(theta_test2(i, 3)));
end

scatter(shaking_forces360(:, 1), shaking_forces360(:, 2));
polar(shaking_forces360(:, 4), shaking_forces360(:, 3));










