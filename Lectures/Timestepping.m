clear all;
clc;
w = 4;
A = [0, 1; -w^2,0];
T = 10;
dt = 0.1;
nsteps = T/dt;
t = linspace(0,T,nsteps+1);

% Exact solution
y_exact = cos(w*t);

% Euler explicit
z_e = zeros(2,nsteps+1);
z_e(1,1) = 1;
z_e(2,1) = 0;

for k=1:nsteps
    z_e(:,k+1) = z_e(:,k) + A*dt*z_e(:,k);
end
error_e = max(abs(y_exact-z_e(1,:)))

% Euler implicit
z_i = zeros(2,nsteps+1);
z_i(1,1) = 1;
z_i(2,1) = 0;

B = eye(2) - dt*A;
for k=1:nsteps
    z_i(:,k+1) = B\z_i(:,k);
end
error_i = max(abs(y_exact-z_i(1,:)))

% Trapezidal
z_t = zeros(2,nsteps+1);
z_t(1,1) = 1;
z_t(2,1) = 0;

C = eye(2)-dt*A/2;
D = eye(2)+dt*A/2;
for k=1:nsteps
    z_t(:,k+1) = C\(D*z_t(:,k));
end
error_t = max(abs(y_exact-z_t(1,:)))

%RK4
z_r = zeros(2,nsteps+1);
z_r(1,1) = 1;
z_r(2,1) = 0;

for k=1:nsteps
    k1 = dt*A*z_r(:,k);
    k2 = dt*A*(z_r(:,k)+k1/2);
    k3 = dt*A*(z_r(:,k)+k2/2);
    k4 = dt*A*(z_r(:,k)+k3);
    z_r(:,k+1)  = z_r(:,k) + (k1+2*k2+2*k3+k4)/6;
end
error_r = max(abs(y_exact-z_r(1,:)))

figure;
plot(t,z_e(1,:),'r');
hold on
plot(t,z_i(1,:),'b');
hold on
plot(t,z_t(1,:),'g');
hold on
plot(t,z_r(1,:),'m+');
hold on
plot(t,y_exact,'k');
axis([0,T,-1.5,1.5]);
legend({'Euler explicit','Euler implicit','Trapezoidal','RK4'},'FontSize', 14)