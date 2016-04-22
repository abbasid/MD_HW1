function [ psi, phi ] = fourbar_angle( theta, linkage )

theta = theta*pi/180;
k_1 = linkage(4)/linkage(1);
k_2 = linkage(4)/linkage(3);
k_3 = (linkage(1).^2 - linkage(2).^2 + linkage(3).^2 + linkage(4).^2)/(2*linkage(1)*linkage(3));
k_4 = linkage(4)/ linkage(2);
k_5 = (linkage(3).^2 - linkage(4).^2 - linkage(1).^2 - linkage(2).^2)/(2*linkage(1)*linkage(2));

A = cos(theta) - k_1 - k_2*cos(theta) + k_3;
B = -2*sin(theta);
C = k_1 - (k_2 + 1)*cos(theta) + k_3;

beta_1 = (-B + sqrt(B.^2 - 4*A*C))/(2*A);
beta_2 = (-B - sqrt(B.^2 - 4*A*C))/(2*A);
psi(1) = 2*atan(beta_1)*180/pi;
psi(2) = 2*atan(beta_2)*180/pi;

D = cos(theta) - k_1 + k_4*cos(theta) + k_5;
E = -2*sin(theta);
F = k_1 + (k_4 - 1)*cos(theta) + k_5;

beta_1 = (-E + sqrt(E.^2 - 4*D*F))/(2*D);
beta_2 = (-E - sqrt(E.^2 - 4*D*F))/(2*D);
phi(1) = 2*atan(beta_1)*180/pi;
phi(2) = 2*atan(beta_2)*180/pi;
end

