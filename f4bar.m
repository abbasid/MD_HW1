function [values,form] = f4bar(r,theta1,theta2,td2,tdd2,sigma,driver)

% program revised from Waldron's

% This function analyzes a four-bar linkage when the crank is the

% driving link.  The input values are:

% theta1,theta2 are angles in degrees

%r(1)    = length  of vector 1 (frame)

%r(2)    = length  of vector 2 (crank)

%r(3)    = length  of vector 3 (coupler)

%r(4)    = length  of vector 4 (rocker or slider offset)

%td2     = crank or coupler angular velocity (rad/sec)

%tdd2    = crank or coupler angular acceleration (rad/sec^2)

%sigma   = +1 or -1.  Identifies assembly mode

%driver = 0 for crank as driver; 1 for coupler as driver

% The results are returned in the vector "values".  The answers are

% stored in values according to the following:

 

%values (1:4,1)   = link position

%values (1:4,2)   = link angles in degrees

%values (1:4,3)   = link angular velocities

%values (1:4,4)   = link angular accelerations

%values (1,5)   = velocity of point Q

%values (2,5)   = velocity of point P

%values (3,5)   = acceleration of point Q

%values (4,5)   = acceleration of point P

%vakyes (4,6)   =absolute values of values(:,4)

%vakyes (4,7)   =angles in degrees of values(:,4)

%form         = assembly flag.  If form = 0, mechanism cannot be

%                 assembled.

 

%convert input data

 

values=zeros(4,7);

% if coupler is the driver, interchange the vetor 3 & 2

if sigma>0, sigma=1;else sigma=-1;end

if driver==1,

   r=[r(1) r(3) r(2) r(4)];

end

rr=r.*r;

fact=pi/180;

theta=zeros(4,1);

td=zeros(4,1);

tdd=zeros(4,1);

theta(1:2)=[theta1 theta2]*fact;

t1=theta(1);

tx=theta(2);

s1=sin(t1);

c1=cos(t1);

sx=sin(tx);

cx=cos(tx);

 

% position calculations

A=2*r(1)*r(4)*c1-2*r(2)*r(4)*cx;

C=rr(1)+rr(2)+rr(4)-rr(3)-2*r(1)*r(2)*(c1*cx+s1*sx);

B=2*r(1)*r(4)*s1-2*r(2)*r(4)*sx;

arg=B*B-C*C+A*A;

if (arg>=0)

   form=1;

 

% Check for the denominator equal to zero

 

   if abs(C-A)>=eps

      t4=2*atan((-B+sigma*sqrt(arg))/(C-A));

      s4=sin(t4);

      c4=cos(t4);

      t3=atan2((r(1)*s1+r(4)*s4-r(2)*sx),(r(1)*c1+r(4)*c4-r(2)*cx));

      s3=sin(t3);

      c3=cos(t3);

   elseif abs(C-A)<eps 

% If the denominator is zero, compute theta(3) first

 

      A=-2*r(1)*r(3)*c1+2*r(2)*r(3)*cx;

      B=-2*r(1)*r(3)*s1+2*r(2)*r(3)*sx;

      C=rr(1)+rr(2)+rr(3)-rr(4)-2*r(1)*r(2)*(c1*cx+s1*sx);

      arg=B*B-C*C+A*A;

      if (arg>=0)

         t3=2*atan((-B-sigma*sqrt(arg))/(C-A));

         s3=sin(t3);

         c3=cos(t3);

         t4=atan2((-r(1)*s1+r(3)*s3+r(2)*sx),(-r(1)*c1+r(3)*c3+r(2)*cx));

         s4=sin(t4);

         c4=cos(t4);

      end  

   end

   theta(3)=t3;

   theta(4)=t4;

  

 

   %velocity calculation

      td(2)=td2;

      AM=[-r(3)*s3, r(4)*s4; -r(3)*c3, r(4)*c4];

      BM=[r(2)*td(2)*sx;r(2)*td(2)*cx];

      CM=AM\BM;

      td(3)=CM(1);

      td(4)=CM(2);

  

   %acceleration calculation

  

      tdd(2)=tdd2;

      BM=[r(2)*tdd(2)*sx+r(2)*td(2)*td(2)*cx+r(3)*td(3)*td(3)*c3-r(4)*td(4)*td(4)*c4;...

       r(2)*tdd(2)*cx-r(2)*td(2)*td(2)*sx-r(3)*td(3)*td(3)*s3+r(4)*td(4)*td(4)*s4];

      CM=AM\BM;

      tdd(3)=CM(1);

      tdd(4)=CM(2);

  %store results in array values

  % coordinates of P and Q

  if driver==1,

     r=[r(1) r(3) r(2) r(4)];

     c2=c3;c3=cx;s2=s3;s3=sx;

     td=[td(1) td(3) td(2) td(4)];

     tdd=[tdd(1) tdd(3) tdd(2) tdd(4)];

     theta=[theta(1) theta(3) theta(2) theta(4)];

  else

     c2=cx;s2=sx;

  end

  values(1,1)=r(1).*exp(i*theta(1));%vector r1

  values(2,1)=r(2).*exp(i*theta(2));%vector rq=r2

  values(3,1)= values(2,1)+ r(3).*exp(i*theta(3));%vector rp=r2+r3

  values(4,1)=r(4).*exp(i*theta(4));%vector r4

  for j=1:4,

     values(j,2)=theta(j)/fact;

     values(j,3)=td(j);

     values(j,4)=tdd(j);

  end % position vectors

  values(1,5)=r(2).*exp(i*theta(2));%velocity for point Q

  values(2,5)=r(4).*exp(i*theta(4));%velocity for point P

  values(3,5)=i*r(2).*(tdd(2)-td(2).*td(2)).*exp(i*theta(2));%accel of Q

  values(4,5)=i*r(4).*(tdd(4)-td(4).*td(4)).*exp(i*theta(4));%accel of P

 

for j=1:4,

values(j,1)=r(j)*exp(i*theta(j));

values(j,2)=theta(j)/fact;

values(j,3)=td(j);

values(j,4)=tdd(j);

end % position vectors

values(1,5)=r(2)*td(2)*exp(i*theta(2));%velocity for point Q

values(2,5)=r(4)*td(4)*exp(i*theta(4));%velocity for point P

values(3,5)=r(2)*(i*tdd(2)-td(2)*td(2))*exp(i*theta(2));%accel of Q

values(4,5)=r(4)*(i*tdd(4)-td(4)*td(4))*exp(i*theta(4));%accel of P

for j=1:4,

values(j,6)=abs(values(j,5)); %absolute values for values(:,4)

values(j,7)=angle(values(j,5))/fact; %angles for values(:,4)

end

else

form=0;

if driver==1,

r=[r(1) r(3) r(2) r(4)];

for j=1:4, values(j,1)=r(j).*exp(i*theta(j));end % position vectors

end

end

 