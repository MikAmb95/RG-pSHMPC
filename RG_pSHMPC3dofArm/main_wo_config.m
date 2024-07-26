clc,close all,clear all
load('workspace_v0.mat')
clc

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;


x0 = [0; 0; 0; deg2rad(90);deg2rad(-90);deg2rad(-90);zeros(ns/2,1)];  % Initial Cofniguration
%xs = [2.2; 0.3; deg2rad(40); deg2rad(90);deg2rad(-90);deg2rad(-90);zeros(ns/2,1)];  % Desired Final Configuration
xs = x0;
xs(1:6) = xs(1:6)+1e-3; 
x_ = [];
%xs = [0.01;0.01;deg2rad(1);deg2rad(89);deg2rad(-91);deg2rad(-89);zeros(ns/2,1)];

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,nc);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables
k0 = [100 100 100 100 100 100;1 1 1 1 1 1];
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
while(mpciter < tsim-2)
    
    
    ss_error = norm((x0-xs),2)
    
    N = N -1;
    [lg,ug,lx,ux] = sh_constraints(N,n_states,12,qmin,qmax,k_min,k_max);
    args.lbg = lg;
    args.ubg = ug;
    args.lbx = lx;
    args.ubx = ux;
    X0 = X0(1:N+1,:);     
    
    if incr <=10
    xs(1:6) = xs(1:6)+1e-3;
    end


    args.p   = [x0;xs]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',ns*(N+1),1);reshape(k0,12,1)];
    
    solver = solver_struc.solver{mpciter+1}; %we call the new solver (pre-computed) at each iterataion
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
    
    k_opt = reshape(full(sol.x(ns*(N+1)+1:end))',2,6); % get controls only from the solution
    %xx1(:,1:ns,mpciter+1)= reshape(full(sol.x(1:ns*(N+1)))',ns,N+1)'; % get solution TRAJECTORY
    
    
    Kp = diag(k_opt(1,:)); Kd = diag(k_opt(2,:));
    Kp_ = [Kp_;k_opt(1,:)]; Kd_ = [Kd_;k_opt(2,:)]; 
    
    
    tau = Kp*(xs(1:6)-x0(1:6))+ Kd*(xs(7:12)-x0(7:12));
   
    
%     u = eval_FB(x0,tau);
%     if (all(all(Su*u<=su)))
%         tau_next = tau;
%         incr = incr + 1;
%         k0 = k_opt;
%     else
%         tau_next = tau_prev;
%         xs(1:6) = xs_prev(1:6);
%         k0 = k0_prev;
%     end
    incr = incr + 1;  
    tau_next = tau;
    k0 = k_opt;
    
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
    
    XS_ = [XS_;xs'];
    
end
toc

ss_error = norm((x0-xs),2)

figure(1)
subplot(3,2,1)
plot(t,xx(1,1:end-1),t,XS_(:,1)); ylabel('x')
subplot(3,2,2)
plot(t,xx(2,1:end-1),t,XS_(:,2)); ylabel('y')
subplot(3,2,3)
plot(t,xx(3,1:end-1),t,XS_(:,3)); ylabel('theta')
subplot(3,2,4)
plot(t,xx(4,1:end-1),t,XS_(:,4)); ylabel('qr1')
subplot(3,2,5)
plot(t,xx(5,1:end-1),t,XS_(:,5)); ylabel('qr2')
subplot(3,2,6)
plot(t,xx(6,1:end-1),t,XS_(:,6)); ylabel('qr3')

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
plot(t,U_nl(:,1),'r',t,[tau_max(1);tau_min(1)]*ones(1,Nf-2),'k--'); ylabel('u1')
subplot(3,2,2)
plot(t,U_nl(:,2),'r',t,[tau_max(2);tau_min(2)]*ones(1,Nf-2),'k--'); ylabel('u2')
subplot(3,2,3)
plot(t,U_nl(:,3),'r',t,[tau_max(3);tau_min(3)]*ones(1,Nf-2),'k--'); ylabel('u3')
subplot(3,2,4)
plot(t,U_nl(:,4),'r',t,[tau_max(4);tau_min(4)]*ones(1,Nf-2),'k--'); ylabel('u4')
subplot(3,2,5)
plot(t,U_nl(:,5),'r',t,[tau_max(5);tau_min(5)]*ones(1,Nf-2),'k--'); ylabel('u5')
subplot(3,2,6)
plot(t,U_nl(:,6),'r',t,[tau_max(6);tau_min(6)]*ones(1,Nf-2),'k--'); ylabel('u6')


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




