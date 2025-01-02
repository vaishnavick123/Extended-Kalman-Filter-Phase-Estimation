%% SIMULATION OF PHASE ESTIMATION OF ECG SIGNAL USING EXTENDED KALMAN FILTER (EKF)
T = 2*pi; % period
ytild = @(t) sin(4*pi*t)+sin(2*pi*t + pi/2);
dytild = @(t) 4*pi*cos(4*pi*t) + 2*pi*cos(2*pi*t+pi/2);
Tx = 1;
d = 3;
J = 10;
Ktemp = -d:J;
K = (1.0001)* Ktemp / Ktemp(end-d);
x = [-1 0 1 1.5 2 2.5 3 3.5 3.5 4.5];
spl = spmak(K,x);
htrue = @(t) fnval(spl,t);

atrue = [0 1]; % true amplitude parameters
anoise = .1; % noise amplitude
L = 5000;

rng('default')
rng(999)
tdat = linspace(0,Tx,L)';
xdat = atrue(1)+atrue(2)*ytild(htrue(tdat))+anoise*randn(L,1);
fs = 1./mean(diff(tdat));
dt = mean(diff(tdat));

% EKF: Raw Phase Estimator
x0 = zeros(4,1); x0(3) = 1;
P0 = diag([1e2 1e2 1e-12 1e-12]);
Q = diag([1e-5 1e0 1e-12 1e-12]);
R = 1;
x1_n_rawekf = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,0);


% Plots
figure(1)
clf
subplot(311)
hold on
plot(tdat,htrue(tdat),'black')
plot(tdat,x1_n_rawekf(1,:),'b')
legend('true','ekf')
title('Phase Estimates')

subplot(312)
hold on
plot(tdat,x1_n_rawekf(1,:)'-htrue(tdat),'b')
yline(0,'--')
legend('ekf')
title('Phase Errors')

subplot(313)
hold on
plot(tdat,xdat,'r.','MarkerSize',.1)
plot(tdat,ytild(x1_n_rawekf(1,:)),'b','LineWidth',1)
legend('data','ekf')
title('Phase Reconstructed Signals')

mse_ekf = sum ( (x1_n_rawekf(1,:)'-htrue(tdat)).^2) / length(tdat);
disp(mse_ekf)





%% Functions
function [x1_n,P1_n] = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,varargin)

phase_ord = length(x0)-3;
if phase_ord > 2
    error('assume phase linear or quadratic (no higher orders for now)')
end

x1_n = zeros(4,length(tdat));
P1_n = [];
% raw_phase = zeros(1,length(tdat));
% freq = zeros(1,length(tdat));
dt = mean(diff(tdat));
if nargin > 9 && ~isempty(varargin{2})
    dt = varargin{2};
end
xprev = x0;
Pprev = P0;

ndat = length(xdat);
if nargin > 10
    ndat = varargin{3};
end

for k = 1:ndat
    zk = xdat(k);
    Ak = diag(ones(1,length(x0))); 
    if phase_ord == 1
        Ak(1,2) = dt;
    elseif phase_ord == 2
        Ak(1,2) = dt; Ak(1,3) = .5*dt^2;
        Ak(2,3) = dt;
    end
    xk_pred = Ak * xprev;
    
    Ck = [];
    innovk = [];
    if phase_ord == 1
        Ck = [xk_pred(3)*dytild(xk_pred(1)) 0 ytild(xk_pred(1)) 1];
        innovk = zk - (xk_pred(3)*ytild(xk_pred(1))+xk_pred(4));
    elseif phase_ord == 2
        Ck = [xk_pred(4)*dytild(xk_pred(1)) 0 0 ytild(xk_pred(1)) 1];
        innovk = zk - (xk_pred(4)*ytild(xk_pred(1))+xk_pred(5));
    end
    
    [xk,Pk] = kalman(xprev,Pprev,0,zk,Ak,0,Ck,Q,R,innovk);
    
    if nargin > 8
        if xk(2) < varargin{1}
            xk(2) = varargin{1};
           % xk(1) = xprev(1) + dt*xk(2);
%         elseif xk(1) < xprev(1)
%             xk(1) = xk_pred(1);
        end
    end
    
    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
    
    xprev = xk;
    Pprev = Pk;
    
end

end

function [x_post,P_post,x_pred,P_pred,K] = kalman(xprev,Pprev,u,y,A,B,C,Q,R,varargin)

% INPUTS: 
%   xprev -- Previous state estimate
%   Pprev -- Previous state estimate covariance
%   u -- Input/control vector
%   y -- Current observed value
%   A -- State transition matrix (or Jacobian of state transition function)
%   B -- Input/control matrix
%   C -- Measurement matrix (or Jacobian of measurement function)
%   Q -- State noise covariance matrix
%   R -- Measurement noise covariance matrix
%   [optional...] 
%   innov -- For EKF, input innovation using nonlinear measurement function
%   xpred -- For EKF, input prediction using nonlinear state transition
%   function
%
% OUTPUTS:
%   x_post -- Posterior state estimate 
%   P_post -- Posterior state estimate covariance
%   x_pred -- Predicted state estimate
%   P_pred -- Predicted state estimate covariance
%   K -- Optimal Kalman gain matrix


% Time update
x_pred = A*xprev+B*u;
% For EKF, Jacobian (A) of state transition function and nonlinear
% prediction is inputted
if nargin > 10
    if ~isempty(varargin{2})
        x_pred = varargin{2};
    end
end

P_pred = A*Pprev*A'+Q;

% Measurement update
S = C*P_pred*C'+R;
K = P_pred*C'*inv(S);
dim = size(K*C,1);
ident = diag(ones(1,dim));

innov = (y-C*x_pred);
% For EKF, Jacobian (C) of measurement function and nonlinear innovation is inputted
if nargin > 9
    if ~isempty(varargin{1})
        innov = varargin{1};
    end
end

x_post = x_pred+K*innov;
P_post = (ident-K*C)*P_pred*(ident-K*C)'+K*R*K';


end







