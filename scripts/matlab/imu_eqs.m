%% Setup

syms R1 R2 R3 dV12 dt12 dt23 v1 v2 dP12 dP23 p1 p2 p3
syms Rwi gI G gI_x JVa12 JPa12 JPa23
syms R1_c R2_c R3_c p1_c p2_c p3_c Rcb pcb
syms s g ba dxy

%% Equations

p1 = s*p1_c + R1_c*pcb;
p2 = s*p2_c + R2_c*pcb;
p3 = s*p3_c + R3_c*pcb;

eq_p3 = p3 == p2 + v2*dt23 + 0.5*g*dt23^2 + R2*dP23;
eq_p2 = p2 == p1 + v1*dt12 + 0.5*g*dt12^2 + R1*dP12;
eq_v2 = v2 == v1 + g*dt12 + R1*dV12;

%% Initial Gravity

eqn = rhs(isolate(eq_p3, v2)) == subs(rhs(eq_v2), v1, rhs(isolate(eq_p2, v1)));
eqn = collect(lhs(eqn) - rhs(eqn), [s, g]);

lambda = simplify(diff(eqn, s));
beta = simplify(diff(eqn, g));

c = eqn;
c = subs(c, s, 0);
c = subs(c, g, 0);
gamma = simplify(-c);

disp('lambda:')
disp(lambda)

disp('beta:')
disp(beta)

disp('gamma:')
disp(gamma)

disp('Linear equation?')
disp(isAlways((lambda)*s + (beta)*g - gamma == eqn))

eqn = subs(eqn, g, Rwi*gI*G - Rwi*gI_x*G*dxy);
eqn = subs(eqn, dV12, dV12 + JVa12*ba);
eqn = subs(eqn, dP12, dP12 + JPa12*ba);
eqn = subs(eqn, dP23, dP23 + JPa23*ba);

eqn = collect(eqn, [s, dxy, ba]);

lambda = simplify(diff(eqn, s));
phi = simplify(diff(eqn, dxy));
xi = simplify(diff(eqn, ba));

c = eqn;
c = subs(c, s, 0);
c = subs(c, dxy, 0);
c = subs(c, ba, 0);
psi = simplify(-c);

disp('lambda:')
disp(lambda)

disp('phi:')
disp(phi)

disp('xi:')
disp(xi)

disp('psi:')
disp(psi)

disp('Linear equation?')
disp(isAlways((lambda)*s + (phi)*dxy + (xi)*ba - psi == eqn))

%% Refinement

eqn_p2 = collect(lhs(eq_p2) - rhs(eq_p2), [v1, v2, g, s]);
eqn_v2 = collect(lhs(eq_v2) - rhs(eq_v2), [v1, v2, g, s]);

H = [
 simplify(diff(eqn_p2, v1)), simplify(diff(eqn_p2, v2)), simplify(diff(eqn_p2, g)), simplify(diff(eqn_p2, s));
 simplify(diff(eqn_v2, v1)), simplify(diff(eqn_v2, v2)), simplify(diff(eqn_v2, g)), simplify(diff(eqn_v2, s));
];

a = eqn_p2;
a = subs(a, v1, 0);
a = subs(a, v2, 0);
a = subs(a, g, 0);
a = subs(a, s, 0);

b = eqn_v2;
b = subs(b, v1, 0);
b = subs(b, v2, 0);
b = subs(b, g, 0);
b = subs(b, s, 0);

z = [
 simplify(-a);
 simplify(-b);
];

disp('H:')
disp(H)

disp('z:')
disp(z)

disp('Linear equation?')
disp(isAlways(H(1, :)*[v1; v2; g; s] - z(1) == eqn_p2))
disp(isAlways(H(2, :)*[v1; v2; g; s] - z(2) == eqn_v2))

%% Velocities

eqn = collect(lhs(eq_p2) - rhs(eq_p2), v1);
A = simplify(diff(eqn, v1));
b = simplify(-subs(eqn, v1, 0));

disp('A:')
disp(A)

disp('b:')
disp(b)

disp('Linear equation?')
disp(isAlways(A*v1 - b == eqn))
