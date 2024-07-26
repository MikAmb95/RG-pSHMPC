
function out = FFS_dynamic_model(x,tau)
    %run('robot_paramenters.m') %function with the dynamic paramenters of the system (robot+target)
    
    
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
    mb=20; %mass base #1 [kg]
    mp=1; %mass EE [kg]
    Izz1=1/12*m1*a1^2; %moment of inertia link #1 [kg*m^2]
    Izz2=1/12*m2*a2^2; %moment of inertia link #2 [kg*m^2]
    Izz3=1/12*m3*a3^2; %moment of inertia link #3 [kg*m^2]
    W1=1; % (the base is modeled as a rectangular) side #1 [m]
    W2=0.5; % (the base is modeled as a rectangular) side #2 [m]
    Izzb=1/12*mb*(4*W1^2+4*W2^2); % moment of inertia base [kg*m^2]
    ex=W1*0-0.1*0; %where the robot is connected to the base (x-direction)
    ey=W2*0-0.1*0; %where the robot is connected to the base (y-direction)

    % States---------------------------------------------------------------
    
    %Position
    xc=x(1); yc=x(2); psi=x(3);
    q1=x(4);    q2=x(5); q3=x(6);
    
    
    %Velocity
    dxc=x(7);    dyc=x(8);    dpsi=x(9);
    dq1=x(10);    dq2=x(11);    dq3=x(12);
    
    
   
    
    
    % Model----------------------------------------------------------------
    m11=m1+m2+m3+mb+mp;
    m12=0;
    m13=ey*m1*cos(psi)+ey*m2*cos(psi)+ey*m3*cos(psi)+ey*mp*cos(psi)-ex*m1*sin(psi)-ex*m2*sin(psi)-ex*m3*sin(psi)-ex*mp*sin(psi)-a2*m3*sin(psi+q1+q2)-ac2*m2*sin(psi+q1+q2)-a2*mp*sin(psi+q1+q2)-ac3*m3*sin(psi+q1+q2+q3)-a3*mp*sin(psi+q1+q2+q3)-a1*m2*sin(psi+q1)-a1*m3*sin(psi+q1)-ac1*m1*sin(psi+q1)-a1*mp*sin(psi+q1);
    m14=-a2*m3*sin(psi+q1+q2)-ac2*m2*sin(psi+q1+q2)-a2*mp*sin(psi+q1+q2)-ac3*m3*sin(psi+q1+q2+q3)-a3*mp*sin(psi+q1+q2+q3)-a1*m2*sin(psi+q1)-a1*m3*sin(psi+q1)-ac1*m1*sin(psi+q1)-a1*mp*sin(psi+q1);
    m15=-a2*m3*sin(psi+q1+q2)-ac2*m2*sin(psi+q1+q2)-a2*mp*sin(psi+q1+q2)-ac3*m3*sin(psi+q1+q2+q3)-a3*mp*sin(psi+q1+q2+q3);
    m16=-sin(psi+q1+q2+q3)*(ac3*m3+a3*mp);
    m22=m1+m2+m3+mb+mp;
    m23=ex*m1*cos(psi)+ex*m2*cos(psi)+ex*m3*cos(psi)+ex*mp*cos(psi)+ey*m1*sin(psi)+ey*m2*sin(psi)+ey*m3*sin(psi)+ey*mp*sin(psi)+a2*m3*cos(psi+q1+q2)+ac2*m2*cos(psi+q1+q2)+a2*mp*cos(psi+q1+q2)+ac3*m3*cos(psi+q1+q2+q3)+a3*mp*cos(psi+q1+q2+q3)+a1*m2*cos(psi+q1)+a1*m3*cos(psi+q1)+ac1*m1*cos(psi+q1)+a1*mp*cos(psi+q1);
    m24=a2*m3*cos(psi+q1+q2)+ac2*m2*cos(psi+q1+q2)+a2*mp*cos(psi+q1+q2)+ac3*m3*cos(psi+q1+q2+q3)+a3*mp*cos(psi+q1+q2+q3)+a1*m2*cos(psi+q1)+a1*m3*cos(psi+q1)+ac1*m1*cos(psi+q1)+a1*mp*cos(psi+q1);
    m25=a2*m3*cos(psi+q1+q2)+ac2*m2*cos(psi+q1+q2)+a2*mp*cos(psi+q1+q2)+ac3*m3*cos(psi+q1+q2+q3)+a3*mp*cos(psi+q1+q2+q3);
    m26=cos(psi+q1+q2+q3)*(ac3*m3+a3*mp);
    m33=Izz1+Izz2+Izz3+Izzb+a1^2*m2+a1^2*m3+a2^2*m3+ac1^2*m1+ac2^2*m2+ac3^2*m3+a1^2*mp+a2^2*mp+a3^2*mp+ex^2*m1+ex^2*m2+ex^2*m3+ey^2*m1+ey^2*m2+ey^2*m3+ex^2*mp+ey^2*mp+2*a1*ac3*m3*cos(q2+q3)+2*a1*a3*mp*cos(q2+q3)+2*a2*ex*m3*cos(q1+q2)+2*ac2*ex*m2*cos(q1+q2)+2*a2*ex*mp*cos(q1+q2)-2*a2*ey*m3*sin(q1+q2)-2*ac2*ey*m2*sin(q1+q2)-2*a2*ey*mp*sin(q1+q2)+2*a1*a2*m3*cos(q2)+2*a1*ac2*m2*cos(q2)+2*a2*ac3*m3*cos(q3)+2*a1*a2*mp*cos(q2)+2*a2*a3*mp*cos(q3)+2*a1*ex*m2*cos(q1)+2*a1*ex*m3*cos(q1)+2*ac1*ex*m1*cos(q1)+2*a1*ex*mp*cos(q1)-2*a1*ey*m2*sin(q1)-2*a1*ey*m3*sin(q1)-2*ac1*ey*m1*sin(q1)-2*a1*ey*mp*sin(q1)+2*ac3*ex*m3*cos(q1+q2+q3)+2*a3*ex*mp*cos(q1+q2+q3)-2*ac3*ey*m3*sin(q1+q2+q3)-2*a3*ey*mp*sin(q1+q2+q3);
    m34=Izz1+Izz2+Izz3+a1^2*m2+a1^2*m3+a2^2*m3+ac1^2*m1+ac2^2*m2+ac3^2*m3+a1^2*mp+a2^2*mp+a3^2*mp+2*a1*ac3*m3*cos(q2+q3)+2*a1*a3*mp*cos(q2+q3)+a2*ex*m3*cos(q1+q2)+ac2*ex*m2*cos(q1+q2)+a2*ex*mp*cos(q1+q2)-a2*ey*m3*sin(q1+q2)-ac2*ey*m2*sin(q1+q2)-a2*ey*mp*sin(q1+q2)+2*a1*a2*m3*cos(q2)+2*a1*ac2*m2*cos(q2)+2*a2*ac3*m3*cos(q3)+2*a1*a2*mp*cos(q2)+2*a2*a3*mp*cos(q3)+a1*ex*m2*cos(q1)+a1*ex*m3*cos(q1)+ac1*ex*m1*cos(q1)+a1*ex*mp*cos(q1)-a1*ey*m2*sin(q1)-a1*ey*m3*sin(q1)-ac1*ey*m1*sin(q1)-a1*ey*mp*sin(q1)+ac3*ex*m3*cos(q1+q2+q3)+a3*ex*mp*cos(q1+q2+q3)-ac3*ey*m3*sin(q1+q2+q3)-a3*ey*mp*sin(q1+q2+q3);
    m35=Izz2+Izz3+a2^2*m3+ac2^2*m2+ac3^2*m3+a2^2*mp+a3^2*mp+a1*ac3*m3*cos(q2+q3)+a1*a3*mp*cos(q2+q3)+a2*ex*m3*cos(q1+q2)+ac2*ex*m2*cos(q1+q2)+a2*ex*mp*cos(q1+q2)-a2*ey*m3*sin(q1+q2)-ac2*ey*m2*sin(q1+q2)-a2*ey*mp*sin(q1+q2)+a1*a2*m3*cos(q2)+a1*ac2*m2*cos(q2)+2*a2*ac3*m3*cos(q3)+a1*a2*mp*cos(q2)+2*a2*a3*mp*cos(q3)+ac3*ex*m3*cos(q1+q2+q3)+a3*ex*mp*cos(q1+q2+q3)-ac3*ey*m3*sin(q1+q2+q3)-a3*ey*mp*sin(q1+q2+q3);
    m36=Izz3+ac3^2*m3+a3^2*mp+a1*ac3*m3*cos(q2+q3)+a1*a3*mp*cos(q2+q3)+a2*ac3*m3*cos(q3)+a2*a3*mp*cos(q3)+ac3*ex*m3*cos(q1+q2+q3)+a3*ex*mp*cos(q1+q2+q3)-ac3*ey*m3*sin(q1+q2+q3)-a3*ey*mp*sin(q1+q2+q3);
    m44=Izz1+Izz2+Izz3+a1^2*m2+a1^2*m3+a2^2*m3+ac1^2*m1+ac2^2*m2+ac3^2*m3+a1^2*mp+a2^2*mp+a3^2*mp+2*a1*ac3*m3*cos(q2+q3)+2*a1*a3*mp*cos(q2+q3)+2*a1*a2*m3*cos(q2)+2*a1*ac2*m2*cos(q2)+2*a2*ac3*m3*cos(q3)+2*a1*a2*mp*cos(q2)+2*a2*a3*mp*cos(q3);
    m45=Izz2+Izz3+a2^2*m3+ac2^2*m2+ac3^2*m3+a2^2*mp+a3^2*mp+a1*ac3*m3*cos(q2+q3)+a1*a3*mp*cos(q2+q3)+a1*a2*m3*cos(q2)+a1*ac2*m2*cos(q2)+2*a2*ac3*m3*cos(q3)+a1*a2*mp*cos(q2)+2*a2*a3*mp*cos(q3);
    m46=Izz3+ac3^2*m3+a3^2*mp+a1*ac3*m3*cos(q2+q3)+a1*a3*mp*cos(q2+q3)+a2*ac3*m3*cos(q3)+a2*a3*mp*cos(q3);
    m55=Izz2+Izz3+a2^2*m3+ac2^2*m2+ac3^2*m3+a2^2*mp+a3^2*mp+2*a2*ac3*m3*cos(q3)+2*a2*a3*mp*cos(q3);
    m56=mp*a3^2+a2*mp*cos(q3)*a3+m3*ac3^2+a2*m3*cos(q3)*ac3+Izz3;
    m66=mp*a3^2+m3*ac3^2+Izz3;
    c1=-a2*dpsi^2*m3*cos(psi+q1+q2)-a2*dq1^2*m3*cos(psi+q1+q2)-a2*dq2^2*m3*cos(psi+q1+q2)-ac2*dpsi^2*m2*cos(psi+q1+q2)-ac2*dq1^2*m2*cos(psi+q1+q2)-ac2*dq2^2*m2*cos(psi+q1+q2)-a2*dpsi^2*mp*cos(psi+q1+q2)-a2*dq1^2*mp*cos(psi+q1+q2)-a2*dq2^2*mp*cos(psi+q1+q2)-ac3*dpsi^2*m3*cos(psi+q1+q2+q3)-ac3*dq1^2*m3*cos(psi+q1+q2+q3)-ac3*dq2^2*m3*cos(psi+q1+q2+q3)-ac3*dq3^2*m3*cos(psi+q1+q2+q3)-a3*dpsi^2*mp*cos(psi+q1+q2+q3)-a3*dq1^2*mp*cos(psi+q1+q2+q3)-a3*dq2^2*mp*cos(psi+q1+q2+q3)-a3*dq3^2*mp*cos(psi+q1+q2+q3)-a1*dpsi^2*m2*cos(psi+q1)-a1*dpsi^2*m3*cos(psi+q1)-a1*dq1^2*m2*cos(psi+q1)-a1*dq1^2*m3*cos(psi+q1)-ac1*dpsi^2*m1*cos(psi+q1)-ac1*dq1^2*m1*cos(psi+q1)-a1*dpsi^2*mp*cos(psi+q1)-a1*dq1^2*mp*cos(psi+q1)-dpsi^2*ex*m1*cos(psi)-dpsi^2*ex*m2*cos(psi)-dpsi^2*ex*m3*cos(psi)-dpsi^2*ex*mp*cos(psi)-dpsi^2*ey*m1*sin(psi)-dpsi^2*ey*m2*sin(psi)-dpsi^2*ey*m3*sin(psi)-dpsi^2*ey*mp*sin(psi)-2*a2*dpsi*dq1*m3*cos(psi+q1+q2)-2*a2*dpsi*dq2*m3*cos(psi+q1+q2)-2*a2*dq1*dq2*m3*cos(psi+q1+q2)-2*ac2*dpsi*dq1*m2*cos(psi+q1+q2)-2*ac2*dpsi*dq2*m2*cos(psi+q1+q2)-2*ac2*dq1*dq2*m2*cos(psi+q1+q2)-2*a2*dpsi*dq1*mp*cos(psi+q1+q2)-2*a2*dpsi*dq2*mp*cos(psi+q1+q2)-2*a2*dq1*dq2*mp*cos(psi+q1+q2)-2*ac3*dpsi*dq1*m3*cos(psi+q1+q2+q3)-2*ac3*dpsi*dq2*m3*cos(psi+q1+q2+q3)-2*ac3*dpsi*dq3*m3*cos(psi+q1+q2+q3)-2*ac3*dq1*dq2*m3*cos(psi+q1+q2+q3)-2*ac3*dq1*dq3*m3*cos(psi+q1+q2+q3)-2*ac3*dq2*dq3*m3*cos(psi+q1+q2+q3)-2*a3*dpsi*dq1*mp*cos(psi+q1+q2+q3)-2*a3*dpsi*dq2*mp*cos(psi+q1+q2+q3)-2*a3*dpsi*dq3*mp*cos(psi+q1+q2+q3)-2*a3*dq1*dq2*mp*cos(psi+q1+q2+q3)-2*a3*dq1*dq3*mp*cos(psi+q1+q2+q3)-2*a3*dq2*dq3*mp*cos(psi+q1+q2+q3)-2*a1*dpsi*dq1*m2*cos(psi+q1)-2*a1*dpsi*dq1*m3*cos(psi+q1)-2*ac1*dpsi*dq1*m1*cos(psi+q1)-2*a1*dpsi*dq1*mp*cos(psi+q1);
    c2=dpsi^2*ey*m1*cos(psi)-a2*dq1^2*m3*sin(psi+q1+q2)-a2*dq2^2*m3*sin(psi+q1+q2)-ac2*dpsi^2*m2*sin(psi+q1+q2)-ac2*dq1^2*m2*sin(psi+q1+q2)-ac2*dq2^2*m2*sin(psi+q1+q2)-a2*dpsi^2*mp*sin(psi+q1+q2)-a2*dq1^2*mp*sin(psi+q1+q2)-a2*dq2^2*mp*sin(psi+q1+q2)-ac3*dpsi^2*m3*sin(psi+q1+q2+q3)-ac3*dq1^2*m3*sin(psi+q1+q2+q3)-ac3*dq2^2*m3*sin(psi+q1+q2+q3)-ac3*dq3^2*m3*sin(psi+q1+q2+q3)-a3*dpsi^2*mp*sin(psi+q1+q2+q3)-a3*dq1^2*mp*sin(psi+q1+q2+q3)-a3*dq2^2*mp*sin(psi+q1+q2+q3)-a3*dq3^2*mp*sin(psi+q1+q2+q3)-a1*dpsi^2*m2*sin(psi+q1)-a1*dpsi^2*m3*sin(psi+q1)-a1*dq1^2*m2*sin(psi+q1)-a1*dq1^2*m3*sin(psi+q1)-ac1*dpsi^2*m1*sin(psi+q1)-ac1*dq1^2*m1*sin(psi+q1)-a1*dpsi^2*mp*sin(psi+q1)-a1*dq1^2*mp*sin(psi+q1)-a2*dpsi^2*m3*sin(psi+q1+q2)+dpsi^2*ey*m2*cos(psi)+dpsi^2*ey*m3*cos(psi)+dpsi^2*ey*mp*cos(psi)-dpsi^2*ex*m1*sin(psi)-dpsi^2*ex*m2*sin(psi)-dpsi^2*ex*m3*sin(psi)-dpsi^2*ex*mp*sin(psi)-2*a2*dpsi*dq1*m3*sin(psi+q1+q2)-2*a2*dpsi*dq2*m3*sin(psi+q1+q2)-2*a2*dq1*dq2*m3*sin(psi+q1+q2)-2*ac2*dpsi*dq1*m2*sin(psi+q1+q2)-2*ac2*dpsi*dq2*m2*sin(psi+q1+q2)-2*ac2*dq1*dq2*m2*sin(psi+q1+q2)-2*a2*dpsi*dq1*mp*sin(psi+q1+q2)-2*a2*dpsi*dq2*mp*sin(psi+q1+q2)-2*a2*dq1*dq2*mp*sin(psi+q1+q2)-2*ac3*dpsi*dq1*m3*sin(psi+q1+q2+q3)-2*ac3*dpsi*dq2*m3*sin(psi+q1+q2+q3)-2*ac3*dpsi*dq3*m3*sin(psi+q1+q2+q3)-2*ac3*dq1*dq2*m3*sin(psi+q1+q2+q3)-2*ac3*dq1*dq3*m3*sin(psi+q1+q2+q3)-2*ac3*dq2*dq3*m3*sin(psi+q1+q2+q3)-2*a3*dpsi*dq1*mp*sin(psi+q1+q2+q3)-2*a3*dpsi*dq2*mp*sin(psi+q1+q2+q3)-2*a3*dpsi*dq3*mp*sin(psi+q1+q2+q3)-2*a3*dq1*dq2*mp*sin(psi+q1+q2+q3)-2*a3*dq1*dq3*mp*sin(psi+q1+q2+q3)-2*a3*dq2*dq3*mp*sin(psi+q1+q2+q3)-2*a1*dpsi*dq1*m2*sin(psi+q1)-2*a1*dpsi*dq1*m3*sin(psi+q1)-2*ac1*dpsi*dq1*m1*sin(psi+q1)-2*a1*dpsi*dq1*mp*sin(psi+q1);
    c3=-a1*ac3*dq2^2*m3*sin(q2+q3)-a1*ac3*dq3^2*m3*sin(q2+q3)-a1*a3*dq2^2*mp*sin(q2+q3)-a1*a3*dq3^2*mp*sin(q2+q3)-a2*dq1^2*ex*m3*sin(q1+q2)-a2*dq2^2*ex*m3*sin(q1+q2)-ac2*dq1^2*ex*m2*sin(q1+q2)-ac2*dq2^2*ex*m2*sin(q1+q2)-a2*dq1^2*ex*mp*sin(q1+q2)-a2*dq2^2*ex*mp*sin(q1+q2)-a1*dq1^2*ey*m2*cos(q1)-a1*dq1^2*ey*m3*cos(q1)-ac1*dq1^2*ey*m1*cos(q1)-a1*dq1^2*ey*mp*cos(q1)-a1*a2*dq2^2*m3*sin(q2)-a1*ac2*dq2^2*m2*sin(q2)-a2*ac3*dq3^2*m3*sin(q3)-a1*a2*dq2^2*mp*sin(q2)-a2*a3*dq3^2*mp*sin(q3)-a1*dq1^2*ex*m2*sin(q1)-a1*dq1^2*ex*m3*sin(q1)-ac1*dq1^2*ex*m1*sin(q1)-a1*dq1^2*ex*mp*sin(q1)-ac3*dq1^2*ey*m3*cos(q1+q2+q3)-ac3*dq2^2*ey*m3*cos(q1+q2+q3)-ac3*dq3^2*ey*m3*cos(q1+q2+q3)-a3*dq1^2*ey*mp*cos(q1+q2+q3)-a3*dq2^2*ey*mp*cos(q1+q2+q3)-a3*dq3^2*ey*mp*cos(q1+q2+q3)-ac3*dq1^2*ex*m3*sin(q1+q2+q3)-ac3*dq2^2*ex*m3*sin(q1+q2+q3)-ac3*dq3^2*ex*m3*sin(q1+q2+q3)-a3*dq1^2*ex*mp*sin(q1+q2+q3)-a3*dq2^2*ex*mp*sin(q1+q2+q3)-a3*dq3^2*ex*mp*sin(q1+q2+q3)-a2*dq1^2*ey*m3*cos(q1+q2)-a2*dq2^2*ey*m3*cos(q1+q2)-ac2*dq1^2*ey*m2*cos(q1+q2)-ac2*dq2^2*ey*m2*cos(q1+q2)-a2*dq1^2*ey*mp*cos(q1+q2)-a2*dq2^2*ey*mp*cos(q1+q2)-2*a2*dpsi*dq1*ey*m3*cos(q1+q2)-2*a2*dpsi*dq2*ey*m3*cos(q1+q2)-2*a2*dq1*dq2*ey*m3*cos(q1+q2)-2*ac2*dpsi*dq1*ey*m2*cos(q1+q2)-2*ac2*dpsi*dq2*ey*m2*cos(q1+q2)-2*ac2*dq1*dq2*ey*m2*cos(q1+q2)-2*a2*dpsi*dq1*ey*mp*cos(q1+q2)-2*a2*dpsi*dq2*ey*mp*cos(q1+q2)-2*a2*dq1*dq2*ey*mp*cos(q1+q2)-2*a1*ac3*dpsi*dq2*m3*sin(q2+q3)-2*a1*ac3*dpsi*dq3*m3*sin(q2+q3)-2*a1*ac3*dq1*dq2*m3*sin(q2+q3)-2*a1*ac3*dq1*dq3*m3*sin(q2+q3)-2*a1*ac3*dq2*dq3*m3*sin(q2+q3)-2*a1*a3*dpsi*dq2*mp*sin(q2+q3)-2*a1*a3*dpsi*dq3*mp*sin(q2+q3)-2*a1*a3*dq1*dq2*mp*sin(q2+q3)-2*a1*a3*dq1*dq3*mp*sin(q2+q3)-2*a1*a3*dq2*dq3*mp*sin(q2+q3)-2*a2*dpsi*dq1*ex*m3*sin(q1+q2)-2*a2*dpsi*dq2*ex*m3*sin(q1+q2)-2*a2*dq1*dq2*ex*m3*sin(q1+q2)-2*ac2*dpsi*dq1*ex*m2*sin(q1+q2)-2*ac2*dpsi*dq2*ex*m2*sin(q1+q2)-2*ac2*dq1*dq2*ex*m2*sin(q1+q2)-2*a2*dpsi*dq1*ex*mp*sin(q1+q2)-2*a2*dpsi*dq2*ex*mp*sin(q1+q2)-2*a2*dq1*dq2*ex*mp*sin(q1+q2)-2*a1*dpsi*dq1*ey*m2*cos(q1)-2*a1*dpsi*dq1*ey*m3*cos(q1)-2*ac1*dpsi*dq1*ey*m1*cos(q1)-2*a1*dpsi*dq1*ey*mp*cos(q1)-2*a1*a2*dpsi*dq2*m3*sin(q2)-2*a1*a2*dq1*dq2*m3*sin(q2)-2*a1*ac2*dpsi*dq2*m2*sin(q2)-2*a2*ac3*dpsi*dq3*m3*sin(q3)-2*a1*ac2*dq1*dq2*m2*sin(q2)-2*a2*ac3*dq1*dq3*m3*sin(q3)-2*a2*ac3*dq2*dq3*m3*sin(q3)-2*a1*a2*dpsi*dq2*mp*sin(q2)-2*a2*a3*dpsi*dq3*mp*sin(q3)-2*a1*a2*dq1*dq2*mp*sin(q2)-2*a2*a3*dq1*dq3*mp*sin(q3)-2*a2*a3*dq2*dq3*mp*sin(q3)-2*a1*dpsi*dq1*ex*m2*sin(q1)-2*a1*dpsi*dq1*ex*m3*sin(q1)-2*ac1*dpsi*dq1*ex*m1*sin(q1)-2*a1*dpsi*dq1*ex*mp*sin(q1)-2*ac3*dpsi*dq1*ey*m3*cos(q1+q2+q3)-2*ac3*dpsi*dq2*ey*m3*cos(q1+q2+q3)-2*ac3*dpsi*dq3*ey*m3*cos(q1+q2+q3)-2*ac3*dq1*dq2*ey*m3*cos(q1+q2+q3)-2*ac3*dq1*dq3*ey*m3*cos(q1+q2+q3)-2*ac3*dq2*dq3*ey*m3*cos(q1+q2+q3)-2*a3*dpsi*dq1*ey*mp*cos(q1+q2+q3)-2*a3*dpsi*dq2*ey*mp*cos(q1+q2+q3)-2*a3*dpsi*dq3*ey*mp*cos(q1+q2+q3)-2*a3*dq1*dq2*ey*mp*cos(q1+q2+q3)-2*a3*dq1*dq3*ey*mp*cos(q1+q2+q3)-2*a3*dq2*dq3*ey*mp*cos(q1+q2+q3)-2*ac3*dpsi*dq1*ex*m3*sin(q1+q2+q3)-2*ac3*dpsi*dq2*ex*m3*sin(q1+q2+q3)-2*ac3*dpsi*dq3*ex*m3*sin(q1+q2+q3)-2*ac3*dq1*dq2*ex*m3*sin(q1+q2+q3)-2*ac3*dq1*dq3*ex*m3*sin(q1+q2+q3)-2*ac3*dq2*dq3*ex*m3*sin(q1+q2+q3)-2*a3*dpsi*dq1*ex*mp*sin(q1+q2+q3)-2*a3*dpsi*dq2*ex*mp*sin(q1+q2+q3)-2*a3*dpsi*dq3*ex*mp*sin(q1+q2+q3)-2*a3*dq1*dq2*ex*mp*sin(q1+q2+q3)-2*a3*dq1*dq3*ex*mp*sin(q1+q2+q3)-2*a3*dq2*dq3*ex*mp*sin(q1+q2+q3);
    c4=a2*dpsi^2*ex*m3*sin(q1+q2)-a1*ac3*dq3^2*m3*sin(q2+q3)-a1*a3*dq2^2*mp*sin(q2+q3)-a1*a3*dq3^2*mp*sin(q2+q3)-a1*ac3*dq2^2*m3*sin(q2+q3)+ac2*dpsi^2*ex*m2*sin(q1+q2)+a2*dpsi^2*ex*mp*sin(q1+q2)+a1*dpsi^2*ey*m2*cos(q1)+a1*dpsi^2*ey*m3*cos(q1)+ac1*dpsi^2*ey*m1*cos(q1)+a1*dpsi^2*ey*mp*cos(q1)-a1*a2*dq2^2*m3*sin(q2)-a1*ac2*dq2^2*m2*sin(q2)-a2*ac3*dq3^2*m3*sin(q3)-a1*a2*dq2^2*mp*sin(q2)-a2*a3*dq3^2*mp*sin(q3)+a1*dpsi^2*ex*m2*sin(q1)+a1*dpsi^2*ex*m3*sin(q1)+ac1*dpsi^2*ex*m1*sin(q1)+a1*dpsi^2*ex*mp*sin(q1)+ac3*dpsi^2*ey*m3*cos(q1+q2+q3)+a3*dpsi^2*ey*mp*cos(q1+q2+q3)+ac3*dpsi^2*ex*m3*sin(q1+q2+q3)+a3*dpsi^2*ex*mp*sin(q1+q2+q3)+a2*dpsi^2*ey*m3*cos(q1+q2)+ac2*dpsi^2*ey*m2*cos(q1+q2)+a2*dpsi^2*ey*mp*cos(q1+q2)-2*a1*ac3*dpsi*dq2*m3*sin(q2+q3)-2*a1*ac3*dpsi*dq3*m3*sin(q2+q3)-2*a1*ac3*dq1*dq2*m3*sin(q2+q3)-2*a1*ac3*dq1*dq3*m3*sin(q2+q3)-2*a1*ac3*dq2*dq3*m3*sin(q2+q3)-2*a1*a3*dpsi*dq2*mp*sin(q2+q3)-2*a1*a3*dpsi*dq3*mp*sin(q2+q3)-2*a1*a3*dq1*dq2*mp*sin(q2+q3)-2*a1*a3*dq1*dq3*mp*sin(q2+q3)-2*a1*a3*dq2*dq3*mp*sin(q2+q3)-2*a1*a2*dpsi*dq2*m3*sin(q2)-2*a1*a2*dq1*dq2*m3*sin(q2)-2*a1*ac2*dpsi*dq2*m2*sin(q2)-2*a2*ac3*dpsi*dq3*m3*sin(q3)-2*a1*ac2*dq1*dq2*m2*sin(q2)-2*a2*ac3*dq1*dq3*m3*sin(q3)-2*a2*ac3*dq2*dq3*m3*sin(q3)-2*a1*a2*dpsi*dq2*mp*sin(q2)-2*a2*a3*dpsi*dq3*mp*sin(q3)-2*a1*a2*dq1*dq2*mp*sin(q2)-2*a2*a3*dq1*dq3*mp*sin(q3)-2*a2*a3*dq2*dq3*mp*sin(q3);
    c5=a1*ac3*dpsi^2*m3*sin(q2+q3)+a1*ac3*dq1^2*m3*sin(q2+q3)+a1*a3*dpsi^2*mp*sin(q2+q3)+a1*a3*dq1^2*mp*sin(q2+q3)+a2*dpsi^2*ex*m3*sin(q1+q2)+ac2*dpsi^2*ex*m2*sin(q1+q2)+a2*dpsi^2*ex*mp*sin(q1+q2)+a1*a2*dpsi^2*m3*sin(q2)+a1*a2*dq1^2*m3*sin(q2)+a1*ac2*dpsi^2*m2*sin(q2)+a1*ac2*dq1^2*m2*sin(q2)-a2*ac3*dq3^2*m3*sin(q3)+a1*a2*dpsi^2*mp*sin(q2)+a1*a2*dq1^2*mp*sin(q2)-a2*a3*dq3^2*mp*sin(q3)+ac3*dpsi^2*ey*m3*cos(q1+q2+q3)+a3*dpsi^2*ey*mp*cos(q1+q2+q3)+ac3*dpsi^2*ex*m3*sin(q1+q2+q3)+a3*dpsi^2*ex*mp*sin(q1+q2+q3)+a2*dpsi^2*ey*m3*cos(q1+q2)+ac2*dpsi^2*ey*m2*cos(q1+q2)+a2*dpsi^2*ey*mp*cos(q1+q2)+2*a1*ac3*dpsi*dq1*m3*sin(q2+q3)+2*a1*a3*dpsi*dq1*mp*sin(q2+q3)+2*a1*a2*dpsi*dq1*m3*sin(q2)+2*a1*ac2*dpsi*dq1*m2*sin(q2)-2*a2*ac3*dpsi*dq3*m3*sin(q3)-2*a2*ac3*dq1*dq3*m3*sin(q3)-2*a2*ac3*dq2*dq3*m3*sin(q3)+2*a1*a2*dpsi*dq1*mp*sin(q2)-2*a2*a3*dpsi*dq3*mp*sin(q3)-2*a2*a3*dq1*dq3*mp*sin(q3)-2*a2*a3*dq2*dq3*mp*sin(q3);
    c6=(ac3*m3+a3*mp)*(a1*dpsi^2*sin(q2+q3)+a1*dq1^2*sin(q2+q3)+a2*dpsi^2*sin(q3)+a2*dq1^2*sin(q3)+a2*dq2^2*sin(q3)+dpsi^2*ey*cos(q1+q2+q3)+dpsi^2*ex*sin(q1+q2+q3)+2*a1*dpsi*dq1*sin(q2+q3)+2*a2*dpsi*dq1*sin(q3)+2*a2*dpsi*dq2*sin(q3)+2*a2*dq1*dq2*sin(q3));
    % Matrices-------------------------------------------------------------
    Mb=[m11,m12,m13;...
        m12,m22,m23;...
        m13,m23,m33];

    Mba=[m14,m15,m16;...
        m24,m25,m26;...
        m34,m35,m36];

    Ma=[m44,m45,m46;...
        m45,m55,m56;...
        m46,m56,m66];
    
    cb=[c1;c2;c3];
    ca=[c4;c5;c6];
    
    M=[Mb,Mba;Mba',Ma];
    c=[cb;ca];
    
    %acc = inv(M)*(-c+tau);
    acc = tau;
    out = [x(7:12);acc]; 
    