
m=0.1;
g=9.81;
l=1;

Jxx=3;
Jyy=3;
Jzz=5;

k1 = 1;
k2 = 1;
k3 = 1;
k4 = 1;

b1 = 0.1;
b2 = 0.1;
b3 = 0.1;
b4 = 0.1;
params = [m;g;l;Jxx;Jyy;Jzz;k1;k2;k3;k4;b1;b2;b3;b4];

initial_states = [0;0;0;0;0;0;0.1;0;0.2;0;0.3;0];
w1 = 0.4;
w2 = 0.2;
w3 = 0.4;
w4 = 0.9;
w = [w1; w2; w3; w4];

sim('quad_dynamics_sim_e.slx');

states_e = [x_out, x_dot_out, y_out, y_dot_out, z_out, z_dot_out, phi_out, phi_dot_out, theta_out, theta_dot_out, psi_out, psi_dot_out];
ddot_e = [phi_ddot_out, theta_ddot_out, psi_ddot_out];

draw3Dmotion(t, x_out, y_out, z_out, phi_out, theta_out, psi_out);















