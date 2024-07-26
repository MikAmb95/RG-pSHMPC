% This script impements a parametric SH-MPC for a 3DoF robotic arm on a
% spacecraft with a Reference Governot to handle input constraints



% last rev: 06/27/2024
clear all; close all; clc

% the system has 6 position (x,y,psi) for the chaser and (q1,q2,q3) robot
% joints. 

% Import CasADi v3.6.3
addpath('C:\Users\Administrator\Documents\ULB\Casadi')
import casadi.*


%Simualation Paramenters
T = 0.01; %samplig time [s] 
N = 250; % prediction horizon (number of steps)


% State constraints (position and velocity constraints)

%min
qmin = [-10 -10 -deg2rad(180) -deg2rad(180) -deg2rad(180) -deg2rad(180)];
qmin(7:12) = [-2 -2 -0.2 -0.2 -0.2 -0.2];

%max
qmax = [10 10 deg2rad(170) deg2rad(150) deg2rad(150) deg2rad(150)];
qmax(7:12) = -qmin(7:12);

%Input constraints (6 torques)
%tau_max =[20; 20;20;20;20;20]; 
%tau_min = -tau_max;

%Gain constraints (6 torques)
k_min =[1; 0.1;1;0.1;1;0.1;1; 0.1;1;0.1;1;0.1]; 
k_max =[200; 20;200;20;200;20;200; 20;200;20;200;20]; 



%Symbolic variable to build the Casadi solver 
xc = SX.sym('xc'); yc = SX.sym('yc'); psi = SX.sym('psi'); %base-position state 
q1 = SX.sym('q1'); q2 = SX.sym('q2'); q3 = SX.sym('q3'); %robot-position state
 
dxc = SX.sym('dxc'); dyc = SX.sym('dyc'); dpsi = SX.sym('dpsi'); %base-velocity state
dq1 = SX.sym('dq1'); dq2 = SX.sym('dq2'); dq3 = SX.sym('dq3'); %robot-velocity state

%pack states
states = [xc;yc;psi;q1;q2;q3;dxc;dyc;dpsi;dq1;dq2;dq3]; n_states = length(states);

%Simbolic variables for the system inputs
u1 = SX.sym('u1'); u2 = SX.sym('u2'); u3 = SX.sym('u3'); u4 = SX.sym('u4');
u5 = SX.sym('u5'); u6 = SX.sym('u6');

controls = [u1;u2;u3;u4;u5;u6]; n_controls = length(controls);


%In this example the decision variables are the gains of a PD controller
K = SX.sym('K',2,n_states/2); % Decision variables (controls)
P = SX.sym('P',n_states + n_states); % parameters (include the initial state and the reference state)

X = SX.sym('X',n_states,(N+1)); % A vector that represents the states over the optimization problem.

%Initialized variables to build the solver(s)
obj = 0; % Objective function
g = [];  % constraints vector

% Cost Function weigth matrices

Q = zeros(n_states,n_states); %state matix 
%Position:
Q(1,1) = 100; Q(2,2) = 100; Q(3,3) = 100; Q(4,4) = 100;
Q(5,5) = 100; Q(6,6) = 100;
%Velocity:
Q(7,7) = 1; Q(8,8) = 1; Q(9,9) = 1; Q(10,10) = 1;
Q(11,11) = 1; Q(12,12) = 1;

R = zeros(n_controls,n_controls); %input matrix
R(1,1) = 0.0005; R(2,2) = 0.0005; 
R(3,3) = 0.0005; R(4,4) = 0.0005; 
R(5,5) = 0.0005; R(6,6) = 0.0005;


st  = X(:,1); % initial state
g = [g;st-P(1:12)]; % initial condition constraints (the fist 12-states = X0)

obj_v =[];

% Here we do the mapping from a Gain Variables to the input Variables
U = SX.sym('U',n_controls,N); 

for k = 1:N
    
    st = X(:,k);
    
    %unpack the state in position and velocity
    x = st(1:n_states/2);
    dx = st(n_states/2+1:end);
    
    %unpack the reference in position and velocity
    ref_p  =  P(n_states+1:3*n_states/2);
    ref_v = P(3*n_states/2+1:end);
    
    %we build the 6 PD controllers
    con1 =  K(1,1)*(ref_p(1) -x(1))+K(2,1)*(ref_v(1)-dx(1));
    con2 = K(1,2)*(ref_p(2) -x(2))+K(2,2)*(ref_v(2)-dx(2));
    con3 = K(1,3)*(ref_p(3) -x(3))+K(2,3)*(ref_v(3)-dx(3));
    con4 = K(1,4)*(ref_p(4) -x(4))+K(2,4)*(ref_v(4)-dx(4));
    con5 = K(1,5)*(ref_p(5) -x(5))+K(2,5)*(ref_v(5)-dx(5));
    con6 = K(1,6)*(ref_p(6) -x(6))+K(2,6)*(ref_v(6)-dx(6));
    con = [con1;con2;con3;con4;con5;con6];
    
    %We map into U
    U(:,k) = con;

end


%here we create the objective functions for each iteration and the equality
%constraints due to the system dynamic
for k = 1:N
    k
    st = X(:,k);  con = U(:,k);   
    obj_v = [obj_v;(st-P(n_states+1:end))'*Q*(st-P(n_states+1:end)) + con'*R*con];
    st_next = X(:,k+1);
    f_value = FFS_dynamic_model(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end


% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(K,2*n_states/2,1)];

% create the solver
nlp_prob = struct('f', sum(obj_v), 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3 (0 = no iter details)
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver_struc = struct;
solver_struc.solver{1} = nlpsol('solver', 'ipopt', nlp_prob,opts);

%build vector of limits constraints
args = struct;
ns = n_states; nc = n_controls;
n_par = size(k_min,1);

%We call a function to build the inequality cnstrints
[lg,ug,lx,ux] = sh_constraints(N,n_states,n_par,qmin,qmax,k_min,k_max);
args.lbg = lg; args.ubg = ug;
args.lbx = lx; args.ubx = ux;


%This is for the SH, we pre-create all the solvers
Nf = N;
g3 = [];
for i = 2:Nf 
        i
        t1 = tic;
        N = N - 1;
        mpciter = i-1;
        OPT_variables2 = [reshape(X(:,1:end-mpciter),n_states*(N+1),1);reshape(K,2*n_states/2,1)];        
        obj2 = sum(obj_v(1:end-mpciter,1));
        g3 = g(1:end-n_states*mpciter,:);       
        nlp_prob2 = struct('f', obj2, 'x', OPT_variables2, 'g', g3, 'p', P);
        solver = nlpsol('solver', 'ipopt', nlp_prob2,opts);
        solver_struc.solver{i} = solver;
end
N = Nf;


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP
% clc,close all,clear all
% load('workspace_v2.mat')
% clc

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
tau_max =[50; 50;50;50;20;20]; 
tau_min = -tau_max;

x0 = [0; 0; 0; deg2rad(90);deg2rad(-90);deg2rad(-90);zeros(ns/2,1)];  % Initial Cofniguration
ref = [0.5; 0.5; deg2rad(10); deg2rad(80);deg2rad(-100);deg2rad(-80);zeros(12/2,1)];  % Desired Final Configuration
xs = x0 + 0.001*(ref-x0);

x_ = [];
%xs = [0.01;0.01;deg2rad(1);deg2rad(89);deg2rad(-91);deg2rad(-89);zeros(ns/2,1)];

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,nc);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables
k0 = [100 100 100 100 100 100;0.1 0.1 0.1 0.1 0.1 0.1];
sim_tim = 20; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
Kp_ = []; Kd_ = [];
tsim = N;


%check input constraints %
Su = zeros(n_controls*2,n_controls);
Su(1,1) = 1; Su(2,1) = -1;Su(3,2) = 1; Su(4,2) = -1;
Su(5,3) = 1; Su(6,3) = -1;Su(7,4) = 1; Su(8,4) = -1;
Su(9,5) = 1; Su(10,5) = -1;Su(11,6) = 1; Su(12,6) = -1;

su=[tau_max;-tau_min];

%%%first MPC %%%%%%
ss_error = norm((x0-xs),2)
args.p   = [x0;xs]; % set the values of the parameters vector
% initial value of the optimization variables
args.x0  = [reshape(X0',ns*(N+1),1);reshape(k0,12,1)];

solver = solver_struc.solver{mpciter+1}; %we call the new solver (pre-computed) at each iterataion
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);


k_opt = reshape(full(sol.x(ns*(N+1)+1:end))',2,6); % get controls only from the solution
%xx1(:,1:ns,mpciter+1)= reshape(full(sol.x(1:ns*(N+1)))',ns,N+1)'; % get solution TRAJECTORY
k0 = k_opt;

Kp = diag(k_opt(1,:)); Kd = diag(k_opt(2,:));
Kp_ = [Kp_;k_opt(1,:)]; Kd_ = [Kd_;k_opt(2,:)]; 


tau = Kp*(xs(1:6)-x0(1:6))+ Kd*(xs(7:12)-x0(7:12));


clear u
x0_ = x0; tau_ = tau;
for i = 1:N
% Apply the control and shift the solution
u(i,:) = eval_FB(x0,tau);
f_value = FFS_dynamic_model_nl(x0,tau);
x0 = x0+ (T*f_value);
x0 = full(x0);
tau = Kp*(xs(1:6)-x0(1:6))+ Kd*(xs(7:12)-x0(7:12));
end
x0 = x0_; tau = tau_;

u_cl= [u_cl ; tau']; %store input
t(mpciter+1) = t0;

% Apply the control and shift the solution
f_value = FFS_dynamic_model_nl(x0,tau);
x0 = x0+ (T*f_value);
x0 = full(x0);
t0 = t0 + T;

xx(:,mpciter+2) = x0;
X0 = reshape(full(sol.x(1:ns*(N+1)))',ns,N+1)'; % get solution TRAJECTORY
% Shift trajectory to initialize the next step
X0 = [X0(2:end,:);X0(end,:)];

mpciter = mpciter + 1;

tau_prev = tau;
xs_prev = xs;
k0_prev = k0;
%start loop %%%%

incr = 0;
XS_ = [];
XS_ = [XS_;xs'];
tic

k_incr = 1;
k_incr_prev = k_incr;

n_else = 0;
while(mpciter < tsim-2)
    
    
    ss_error = norm((x0-xs),2)
    
    N = N -1;
    [lg,ug,lx,ux] = sh_constraints(N,n_states,12,qmin,qmax,k_min,k_max);
    args.lbg = lg;
    args.ubg = ug;
    args.lbx = lx;
    args.ubx = ux;
    X0 = X0(1:N+1,:);     
    
    
    xs = xs_prev+k_incr*(ref-xs_prev);
    
    
    args.p   = [x0;xs]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',ns*(N+1),1);reshape(k0,12,1)];
    
    solver = solver_struc.solver{mpciter+1}; %we call the new solver (pre-computed) at each iterataion
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
    
    k_opt = reshape(full(sol.x(ns*(N+1)+1:end))',2,6); % get controls only from the solution
    %xx1(:,1:ns,mpciter+1)= reshape(full(sol.x(1:ns*(N+1)))',ns,N+1)'; % get solution TRAJECTORY
    %k0 = k_opt;
    
    Kp = diag(k_opt(1,:)); Kd = diag(k_opt(2,:));
    Kp_ = [Kp_;k_opt(1,:)]; Kd_ = [Kd_;k_opt(2,:)]; 
    
    
    tau = Kp*(xs(1:6)-x0(1:6))+ Kd*(xs(7:12)-x0(7:12));
   
    clear u
    x0_ = x0; tau_ = tau;
    for i = 1:N
    % Apply the control and shift the solution
    u(i,:) = eval_FB(x0,tau);
    f_value = FFS_dynamic_model_nl(x0,tau);
    x0 = x0+ (T*f_value);
    x0 = full(x0);
    tau = Kp*(xs(1:6)-x0(1:6))+ Kd*(xs(7:12)-x0(7:12));
    end
    x0 = x0_; tau = tau_;
    
    if (all(all(Su*u'<=su)))
        tau_next = tau;
        k_incr = 1;
        k0 = k_opt;
        n_else = 0;
    else
        tau_next = tau_prev;
        k_incr = k_incr_prev/2;
        xs = xs_prev;
        k0 = k0_prev;
        n_else = n_else +1;
    end
        
    
    u_cl= [u_cl ; tau_next']; %store input
    t(mpciter+1) = t0;
    
    % Apply the control and shift the solution
    f_value = FFS_dynamic_model_nl(x0,tau_next);
    x0 = x0+ (T*f_value);
    x0 = full(x0);
    t0 = t0 + T;
    
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:ns*(N+1)))',ns,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    
    mpciter
    mpciter = mpciter + 1;
    
    tau_prev = tau_next;
    xs_prev = xs;
    k0_prev = k0;
    k_incr_prev = k_incr;
    
    XS_ = [XS_;xs'];
    
end
toc

ss_error = norm((x0-xs),2)

figure(1)
grid on
subplot(3,2,1)
plot(t,xx(1,1:end-1),t,XS_(:,1)); ylabel('x[m]');xlabel('time[s]');grid on
subplot(3,2,2)
plot(t,xx(2,1:end-1),t,XS_(:,2)); ylabel('y[m]');xlabel('time[s]');grid on
subplot(3,2,3)
plot(t,rad2deg(xx(3,1:end-1)),t,rad2deg(XS_(:,3))); ylabel('theta[deg]');xlabel('time[s]');grid on
subplot(3,2,4)
plot(t,rad2deg(xx(4,1:end-1)),t,rad2deg(XS_(:,4))); ylabel('qr1[deg]');xlabel('time[s]');grid on
subplot(3,2,5)
plot(t,rad2deg(xx(5,1:end-1)),t,rad2deg(XS_(:,5))); ylabel('qr2[deg]');xlabel('time[s]');grid on
subplot(3,2,6)
plot(t,rad2deg(xx(6,1:end-1)),t,rad2deg(XS_(:,6))); ylabel('qr3[deg]');xlabel('time[s]');grid on

figure(2)
subplot(3,2,1)
plot(t,xx(7,1:end-1),t,0*ones(1,Nf-2)); ylabel('x')
subplot(3,2,2)
plot(t,xx(8,1:end-1),t,0*ones(1,Nf-2)); ylabel('y')
subplot(3,2,3)
plot(t,xx(9,1:end-1),t,0*ones(1,Nf-2)); ylabel('theta')
subplot(3,2,4)
plot(t,xx(10,1:end-1),t,0*ones(1,Nf-2)); ylabel('qr1')
subplot(3,2,5)
plot(t,xx(11,1:end-1),t,0*ones(1,Nf-2)); ylabel('qr2')
subplot(3,2,6)
plot(t,xx(12,1:end-1),t,0*ones(1,Nf-2)); ylabel('qr3')

U_nl = [];
for i = 1:tsim-2
    u = eval_FB(xx(:,i),u_cl(i,:)');
    U_nl = [U_nl;u'];
end

figure(3)
subplot(3,2,1)
plot(t,U_nl(:,1),'r',t,[tau_max(1);tau_min(1)]*ones(1,Nf-2),'k--'); ylabel('u1');xlabel('time[s]');grid on
subplot(3,2,2)
plot(t,U_nl(:,2),'r',t,[tau_max(2);tau_min(2)]*ones(1,Nf-2),'k--'); ylabel('u2');xlabel('time[s]');grid on
subplot(3,2,3)
plot(t,U_nl(:,3),'r',t,[tau_max(3);tau_min(3)]*ones(1,Nf-2),'k--'); ylabel('u3');xlabel('time[s]');grid on
subplot(3,2,4)
plot(t,U_nl(:,4),'r',t,[tau_max(4);tau_min(4)]*ones(1,Nf-2),'k--'); ylabel('u4');xlabel('time[s]');grid on
subplot(3,2,5)
plot(t,U_nl(:,5),'r',t,[tau_max(5);tau_min(5)]*ones(1,Nf-2),'k--'); ylabel('u5');xlabel('time[s]');grid on
subplot(3,2,6)
plot(t,U_nl(:,6),'r',t,[tau_max(6);tau_min(6)]*ones(1,Nf-2),'k--'); ylabel('u6');xlabel('time[s]');grid on


%print_system_config


function [lbg,ubg,lbx,ubx] = sh_constraints(N,ns,np,qmin,qmax,k_min,k_max)

%INPUT: N = Prediction Horizon; ns =  dim states; np = dimension par
%       qmin = state low limit; qmax = state up limit;
%       k_min = gain low limit; k_max = gain up limit;

lbg(1:ns*(N+1)) = 0;  % -1e-20  % Equality constraints
ubg(1:ns*(N+1)) = 0;  % 1e-20   % Equality constraints

for i =1:ns
    lbx(i:ns:ns*(N+1),1) = qmin(i); %state q1 lower bound
    ubx(i:ns:ns*(N+1),1) = qmax(i); %state q1 upper bound
end



for i =1:np
    lbx(ns*(N+1)+i:np:ns*(N+1)+np,1) = k_min(i); %u1 lower bound
    ubx(ns*(N+1)+i:np:ns*(N+1)+np,1) = k_max(i); %v upper bound
end

end




