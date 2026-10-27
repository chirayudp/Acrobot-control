clear
clc


m1 = 1;
m2 = 1;
l1 = 1;
l2 = 2;
I1 = 0.083;
I2 = 0.33;
g = 9.81;
taumax = 20;


syms q1 q2 q1dot q2dot real


M = [ I1 + I2 + m1*(l1/2)^2 + m2*(l1^2 + (l2/2)^2 + l1*l2*cos(q2)), I2 + m2*((l2/2)^2 + (l1*l2/2)*cos(q2));
    I2 + m2*((l2/2)^2 + (l1*l2/2)*cos(q2)),                        I2 + m2*(l2/2)^2 ];


C = [ -m2*l1*l2*sin(q2)*q2dot,        -m2*(l1*l2/2)*sin(q2)*q2dot;
    m2*(l1*l2/2)*sin(q2)*q1dot,     0 ];


tau_g1 = -(m1*l1/2 + m2*l1)*g*sin(q1) - (m2*l2/2)*g*sin(q1+q2);
tau_g2 = -(m2*l2/2)*g*sin(q1+q2);
tau_g = [tau_g1; tau_g2];


B = [0; 1];


disp('Mass matrix M(q):')
disp(M)

disp('Coriolis matrix C(q,qdot):')
disp(C)

disp('Gravity torque vector tau_g(q):')
disp(tau_g)


q1_vals = [0, 0, pi, pi];
q2_vals = [0, pi, 0, pi];
names = ["down-down", "down-foldedup", "upright", "up-foldeddown"];


Gjac = jacobian(tau_g, [q1, q2]);


figure


for k = 1:4

    q1e = q1_vals(k);
    q2e = q2_vals(k);

    Me = double(subs(M, {q1,q2}, {q1e,q2e}));
    Ge = double(subs(Gjac, {q1,q2}, {q1e,q2e}));

    A21 = Me \ Ge;
    Bvec = Me \ B;

    Amat = [zeros(2), eye(2); A21, zeros(2)];
    Bmat = [0; 0; Bvec];
    Cmat = eye(2,4);
    Dmat = [0;0];

    fprintf('\nEquilibrium %s : q1 = %.3f , q2 = %.3f\n', names(k), q1e, q2e)

    disp('A matrix:')
    disp(Amat)

    disp('B matrix:')
    disp(Bmat)

    poles = eig(Amat);
    disp('Poles of linearized system:')
    disp(poles.')

    sysLin = ss(Amat, Bmat, Cmat, Dmat);
    tfLin = tf(sysLin);

    disp('Transfer functions theta1/tau and theta2/tau:')
    tfLin

    subplot(2,2,k)
    pzmap(sysLin)
    title(names(k))
    grid on

end

sgtitle('Pole-zero maps at each equilibrium')