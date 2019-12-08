clc; clear;

l1 = 1;
l2 = 0.5;

m1 = 8;
m2 = 3;
m3 = 1;

% initialize joint variables
syms t real
q1 = 2*sin(t);
q2 = cos(t);
q3 = 3*sin(t);

I1 = [(1/12)*m1*l1^2 0 0;
       0 (1/12)*m1*l1^2 0;
       0 0 0];

I2 = [0 0 0;
       0 (1/12)*m2*l2^2 0;
       0 0 (1/12)*m2*l2^2];
  
I3 = [(1/12)*m3*q3^2 0 0;
       0 (1/12)*m3*q3^2 0;
       0 0 0];

% T = trotz(theta)*transl(0,0,d)*transl(a,0,0)*trotx(alpha);
T01 = trotz(q1) * transl(0, 0, l1)*trotx(pi/2);
T12 = trotz(q2) * transl(l2, 0, 0)*troty(pi/2);
T23 = transl(0, 0, q3);
T02 = T01*T12;
T03 = T02*T23;

R01 = T01(1:3, 1:3);
R12 = T12(1:3, 1:3);
R23 = T23(1:3, 1:3);
R02 = T02(1:3, 1:3);
R03 = T03(1:3, 1:3);

z01 = R01(1:3, 3);
z02 = R02(1:3, 3);
z03 = R03(1:3, 3);

r01 = T01(1:3, 4);
r12 = T12(1:3, 4);
r23 = T23(1:3, 4);

rc01 = r01/2;
rc12 = r12/2;
rc23 = r23/2;

q = [q1 q2 q3];
dq = [diff(q1,t), diff(q2,t), diff(q3,t)];
ddq = [diff(q1,t,2), diff(q2,t,2), diff(q3,t,2)];
z0 = [0 0 1]';
omega0 = [0 0 0]';
domega0 = [0 0 0]';

dp0 = [0 0 0]';
ddp0 = [0 0 0]';

g0 = -[0 0 1]' * 9.82;

omega1 = R01'*(omega0 + dq(1)*z0);
omega2 = R12'*(omega1 + dq(2)*z0);
omega3 = R23'*omega2; % because it is the prismatic joint

domega1 = R01'*(domega0 + ddq(1)*z0 + cross(dq(1)*omega0, z0));
domega2 = R12'*(domega1 + ddq(2)*z0 + cross(dq(2)*omega1, z0));
domega3 = R23'*domega2;

dp1 = R01'*dp0 + cross(omega1, r01);
dp2 = R12'*dp1 + cross(omega2, r12);
dp3 = R23'*(dp2 + dq(3)*z0) + cross(omega3, r23);

dpc1 = dp1 + cross(omega1, -rc01);
dpc2 = dp2 + cross(omega2, -rc12);
dpc3 = dp3 + cross(omega3, -rc23);

ddp1 = R01'*ddp0 + cross(domega1, r01) + cross(omega1, cross(omega1, r01));
ddp2 = R12'*ddp1 + cross(domega2, r12) + cross(omega2, cross(omega2, r12));
ddp3 = R23'*(ddp2 + ddq(3)*z0) + cross(2*dq(3)*omega3, R23'*z0) + cross(domega3, r23) + cross(omega3, cross(omega3, r23));

ddpc1 = ddp1 + cross(domega1, -rc01) + cross(omega1, cross(omega1, -rc01));
ddpc2 = ddp2 + cross(domega2, -rc12) + cross(omega2, cross(omega2, -rc12));
ddpc3 = ddp3 + cross(domega3, -rc23) + cross(omega3, cross(omega3, -rc23));

centrifugal1 = cross(omega1, cross(omega1, r01)) + cross(omega1, cross(omega1, -rc01));
centrifugal2 = cross(omega2, cross(omega2, r12)) + cross(omega2, cross(omega2, -rc12));
centrifugal3 = cross(omega3, cross(omega3, r23)) + cross(omega3, cross(omega3, -rc23));

coriolis3 = cross(2*dq(3)*omega3, R23'*z0);

g1 = R01'*g0;
g2 = R12'*g1;
g3 = R23'*g2;

time = 0:0.01:10;

myplot(1, dp1, time, 'Time', 'End-effector velocity')
myplot(2, dp2, time, 'Time', 'End-effector velocity')
myplot(3, dp3, time, 'Time', 'End-effector velocity')

myplot(4, dpc1, time, 'Time', 'Center of mass velocity')
myplot(5, dpc2, time, 'Time', 'Center of mass velocity')
myplot(6, dpc3, time, 'Time', 'Center of mass velocity')

myplot(7, ddp1, time, 'Time')
myplot(8, ddp2, time, 'Time')
myplot(9, ddp3, time, 'Time')

myplot(10, ddpc1, time, 'Time')
myplot(11, ddpc2, time, 'Time')
myplot(12, ddpc3, time, 'Time', 'Center of mass acceleration', '$\ddot{p}_{Cx}$', '$\ddot{p}_{Cy}$', '$\ddot{p}_{Cz}$', 'ddpc3')

myplot(13, omega1, time, 'Time')
myplot(14, omega2, time, 'Time')
myplot(15, omega3, time, 'Time')

myplot(16, domega1, time, 'Time')
myplot(17, domega2, time, 'Time')
myplot(18, domega3, time, 'Time')

myplot(19, centrifugal1, time, 'Time', 'Centrifugal term')
myplot(20, centrifugal2, time, 'Time', 'Centrifugal term')
myplot(21, centrifugal3, time, 'Time', 'Centrifugal term')

myplot(22, coriolis3, time, 'Time')

myplot(23, g1, time, 'Time', 'Gravity term')
myplot(24, g2, time, 'Time', 'Gravity term')
myplot(25, g3, time, 'Time', 'Gravity term')

f3 = m3*ddpc3 - m3*g3;
f2 = R23*f3 + m2*ddpc2 - m2*g2;
f1 = R12*f2 + m1*ddpc1 - m1*g1;
  
mu3 = -cross(f3, rc23) + I3*domega3 + cross(omega3, I3*omega3);
mu2 = -cross(f2, rc12) + R23*mu3 + cross(R23*f3, -rc12) + I2*domega2 + cross(omega2, I2*omega2);
mu1 = -cross(f1, rc01) + R12*mu2 + cross(R12*f2, -rc01) + I1*domega1 + cross(omega1, I1*omega1);

myplot(26, f1, time, 'Time', 'Force')
myplot(27, f2, time, 'Time', 'Force')
myplot(28, f3, time, 'Time', 'Force')

myplot(29, mu1, time, 'Time', 'Force')
myplot(30, mu2, time, 'Time', 'Force')
myplot(31, mu3, time, 'Time', 'Force')

tau1 = dot(mu1, R01'*z0);
tau2 = dot(mu2, R12'*z0);
tau3 = dot(f3, R23'*z0);
tau = [tau1 tau2 tau3]';

myplot(32, tau, time, 'Time', 'Generic forces')
myplot(33, q', time, 'Time', 'Joint variables')

f_centrifugal3 = m3*centrifugal3;
f_centrifugal2 = m2*centrifugal2;
f_centrifugal1 = m1*centrifugal1;

mu_centrifugal3 = -cross(f_centrifugal3, rc23);
mu_centrifugal2 = -cross(f_centrifugal2, rc12) + R23*mu_centrifugal3 + cross(R23*f_centrifugal3, -rc12);
mu_centrifugal1 = -cross(f_centrifugal1, rc01) + R12*mu_centrifugal2 + cross(R12*f_centrifugal2, -rc01);

tau1 = dot(mu_centrifugal1, R01'*z0);
tau2 = dot(mu_centrifugal2, R12'*z0);
tau3 = dot(f_centrifugal3, R23'*z0);
tau = [tau1 tau2 tau3]';
myplot(34, tau, time, 'Time')

f_coriolis3 = m3*coriolis3;
f_coriolis2 = R23*f3;
f_coriolis1 = R12*f2;

mu_coriolis3 = -cross(f_coriolis3, rc23);
mu_coriolis2 = -cross(f_coriolis2, rc12) + R23*mu_coriolis3 + cross(R23*f_coriolis3, -rc12);
mu_coriolis1 = -cross(f_coriolis1, rc01) + R12*mu_coriolis2 + cross(R12*f_coriolis2, -rc01);

tau1 = dot(mu_coriolis1, R01'*z0);
tau2 = dot(mu_coriolis2, R12'*z0);
tau3 = dot(f_coriolis3, R23'*z0);
tau = [tau1 tau2 tau3]';
myplot(35, tau, time, 'Time', 'Сoriolis term', '$\tau_1$', '$\tau_2$', '$\tau_3$', 'tau_coriolis')

f_gravity3 = -m3*g3;
f_gravity2 = R23*f_gravity3 - m2*g2;
f_gravity1 = R12*f_gravity2 - m1*g1;
  
mu_gravity3 = -cross(f_gravity3, rc23) + I3*domega3 + cross(omega3, I3*omega3);
mu_gravity2 = -cross(f_gravity2, rc12) + R23*mu_gravity3 + cross(R23*f_gravity3, -rc12) + I2*domega2 + cross(omega2, I2*omega2);
mu_gravity1 = -cross(f_gravity1, rc01) + R12*mu_gravity2 + cross(R12*f_gravity2, -rc01) + I1*domega1 + cross(omega1, I1*omega1);

tau1 = dot(mu_gravity1, R01'*z0);
tau2 = dot(mu_gravity2, R12'*z0);
tau3 = dot(f_gravity3, R23'*z0);
tau = [tau1 tau2 tau3]';
myplot(36, tau, time, 'Time', 'Gravity term')
