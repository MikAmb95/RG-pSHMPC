

% Characteristics of the robot---------------------------------------------
a1=2; %length of the link #1 [m]
a2=2; %length of the link #2 [m]
a3=2; %length of the link #3 [m]
ac1=a1/2; %CoM link #1 [m]
ac2=a2/2; %CoM link #2 [m]
ac3=a3/2; %CoM link #3 [m]
m1=2; %mass link #1 [kg]
m2=2; %mass link #2 [kg]
m3=2; %mass link #3 [kg]
mb=5; %mass link #1 [kg]
mp=1; %mass EE [kg]
Izz1=1/12*m1*a1^2; %moment of inertia link #1 [kg*m^2]
Izz2=1/12*m2*a2^2; %moment of inertia link #2 [kg*m^2]
Izz3=1/12*m3*a3^2; %moment of inertia link #3 [kg*m^2]
W1=1; % (the base is modeled as a rectangular) side #1 [m]
W2=0.5; % (the base is modeled as a rectangular) side #2 [m]
Izzb=1/12*mb*(4*W1^2+4*W2^2); % moment of inertia base [kg*m^2]
ex=W1*0-0.1*0; %where the robot is connected to the base (x-direction)
ey=W2*0-0.1*0; %where the robot is connected to the base (y-direction)


x = xx;

xc=x(1,:); %x_base
yc=x(2,:); %y_base
psi=x(3,:); %orient_base
q1=x(4,:); %q1 robot
q2=x(5,:); %q2 robot
q3=x(6,:); %q3 robot


xb=xc+ex*cos(psi)+ey*sin(psi); %x_base
yb=yc-ey*cos(psi)+ex*sin(psi); %y_base

xe1=xc+a1*cos(psi+q1)+ex*cos(psi)+ey*sin(psi);
ye1=yc+a1*sin(psi+q1)-ey*cos(psi)+ex*sin(psi);

xe2=xc+a1*cos(psi+q1)+ex*cos(psi)+ey*sin(psi)+a2*cos(psi+q1+q2);
ye2=yc+a1*sin(psi+q1)-ey*cos(psi)+ex*sin(psi)+a2*sin(psi+q1+q2);

xe3=xc+a3*cos(psi+q1+q2+q3)+a1*cos(psi+q1)+ex*cos(psi)+ey*sin(psi)+a2*cos(psi+q1+q2);
ye3=yc+a3*sin(psi+q1+q2+q3)+a1*sin(psi+q1)-ey*cos(psi)+ex*sin(psi)+a2*sin(psi+q1+q2);


hold on;grid on;box on

axis equal

i = size(x,2);

i = 1;

for i = 1:size(x,2)

x11=xc(i)+W1.*cos(psi(i))-W2.*sin(psi(i));
y11=yc(i)+W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
x21=xc(i)-W1.*cos(psi(i))+W2.*cos(psi(i)+pi/2);
y21=yc(i)-W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
x31=xc(i)+W1.*cos(psi(i))+W2.*sin(psi(i));
y31=yc(i)+W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);
x41=xc(i)-W1.*cos(psi(i))-W2.*cos(psi(i)+pi/2);
y41=yc(i)-W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);
Line1 = line([xb(i);xe1(i)],[yb(i);ye1(i)],'color',[0.3,0.75,1]);
Line2 = line([xe1(i);xe2(i)],[ye1(i);ye2(i)],'color',[0.3,0.75,1]);
Line3 = line([xe2(i);xe3(i)],[ye2(i);ye3(i)],'color',[0.3,0.75,1]);
Line4 = line([x11;x31],[y11;y31],'color',[0.9,0.7,0.15]);
Line5 = line([x41;x31],[y41;y31],'color',[0.9,0.7,0.15]);
Line6 = line([x21;x41],[y21;y41],'color',[0.9,0.7,0.15]);
Line7 = line([x21;x11],[y21;y11],'color',[0.9,0.7,0.15]);


pause(0.1)

if i == size(x,2); break 
end

% delete(Line1);
% delete(Line2);
% delete(Line3);
% delete(Line4);
% delete(Line5);
% delete(Line6);
% delete(Line7);



end