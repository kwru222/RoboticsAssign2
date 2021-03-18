
function [t1,t2,t3]=inverse_dynamics_one(a,b,c,d,e,f,g,h,i)
syms p1(t) p2(t) p3(t) Dp1(t) Dp2(t) Dp3(t) LL1 LL2 LL3 LLm_one LLm_two LLm_three;
%%%%%%%%%%%%%%%%%%%%........ Rotational and Translational
%%%%%%%%%%%%%%%%%%%%frames.................................................
Rz1 = [cos(p1(t)) -sin(p1(t)) 0 0; sin(p1(t)) cos(p1(t)) 0 0; 0 0 1 0; 0 0 0 1];
Tz1 = [1 0 0 0; 0 1 0 0; 0 0 1 LL1; 0 0 0 1];
Rx1 = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
T01 = Rz1*Tz1*Rx1*[1 0 0 LL1/2;0 1 0 0;0 0 1 0;0 0 0 1];
R01 = T01(1:3,1:3);

Rz2 = [cos(p2(t)) -sin(p2(t)) 0 0; sin(p2(t)) cos(p2(t)) 0 0; 0 0 1 0; 0 0 0 1];
Tx2 = [1 0 0 LL2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
T02 = T01*Rz2*Tx2*[1 0 0 LL1/2;0 1 0 0;0 0 1 0;0 0 0 1];
R02 = T02(1:3,1:3);

Rz3 = [cos(p3(t)) -sin(p3(t)) 0 0; sin(p3(t)) cos(p3(t)) 0 0; 0 0 1 0; 0 0 0 1];
Tx3 = [1 0 0 LL3; 0 1 0 0; 0 0 1 0; 0 0 0 1];
T03 = T02*Rz3*Tx3*[1 0 0 LL1/2;0 1 0 0;0 0 1 0;0 0 0 1];
R03 = T03(1:3,1:3);
d1=T01(1:3,4);
d2=T02(1:3,4);
d3 = T03(1:3,4);
d1_d0=d1;
d2_d0=d2;
d2_d1=d2-d1;
d3_d0=d3;
d3_d1=d3-d1;
d3_d2=d3-d2;


%%%%%%..........Moments of Inertia as obtained from the urdf file.......... 
I_one = [1 0 0; 0 .083 0; 0 0 1];
I_two = [1 0 0; 0 .083 0; 0 0 1];
I_three = [1 0 0; 0 .33 0; 0 0 1];

DQ = [Dp1(t); Dp2(t); Dp3(t)];
m_one = 1;
m_two = 1;
m_three = 1;
g=9.81;



%%%%%%%%%.............. Calculating the
%%%%%%%%%Jacobians...........................
J11 = [cross([0; 0; 1], d1_d0) [0 0;0 0;0 0]];
J12 = [cross(R01(:,3), d2_d0) cross(R02(:,3), d2_d1) [0;0;0]];
J13 = [cross([0;0;1], d3_d0) cross(R01(:,3), d3_d1) cross(R02(:,3), d3_d2)];

J21 = [[0; 0; 1] [0;0;0] [0;0;0]];
J22 = [[0; 0; 1] R01(:,3) [0;0;0]];
J23 = [[0;0;1] R01(:,3) R02(:,3)];

J = [J11 J12 J13; J21 J22 J23];
J=simplify(J);


%%%%%%%...............Kinetic Energy.........
KE = 0.5*transpose(DQ)*(transpose(J21)*R01*I_one*transpose(R01)*J21+m_one*transpose(J11)*J11+transpose(J22)*R02*I_two*transpose(R02)*J22+m_two*transpose(J12)*J12+transpose(J23)*R03*I_three*transpose(R03)*J23+m_three*transpose(J13)*J13)*DQ;

%%%%%%%...............Potential Energy calculation......................
PE= (m_one*g*LL1/2)+m_two*g*(LL1+LL2*sin(p2(t))/2)+m_three*g*(LL1+LL2*sin(p2(t))+(LL3/2)*sin(p3(t)));

%%%%%%%...............Calculating the Lagrangian.......................
Lagrange=KE-PE;
Lagrange=simplify(Lagrange);
D=(transpose(J21)*R01*I_one*transpose(R01)*J21+m_one*transpose(J11)*J11+transpose(J22)*R02*I_two*transpose(R02)*J22+m_two*transpose(J12)*J12+transpose(J23)*R03*I_three*transpose(R03)*J23+m_three*transpose(J13)*J13);
dLL1_dp1=diff(Lagrange,Dp1(t));
dLL1_dp2=diff(Lagrange,Dp2(t));
dLL1_dp3=diff(Lagrange,Dp3(t));
dLL1_p1=diff(Lagrange,p1(t));
dLL1_p2=diff(Lagrange,p2(t));
dLL1_p3=diff(Lagrange,p3(t));

%%%%%%%%%.................Calculating the Taos.............................
tao_one=diff(dLL1_dp1,t)-dLL1_p1;
tao_two=diff(dLL1_dp2,t)-dLL1_p2;
tao_three=diff(dLL1_dp2,t)-dLL1_p3;
t1=subs(tao_one,[p1(t),p2(t),p3(t),Dp1(t),Dp2(t),Dp3(t),diff(p1(t), t),diff(p2(t), t),diff(p2(t), t),diff(Dp1(t), t),diff(Dp2(t), t),diff(Dp3(t), t),LL1,LL2,LL3],[a,b,c,d,e,f,d,e,f,g,h,i,1,1,1])
t2=subs(tao_two,[p1(t),p2(t),p3(t),Dp1(t),Dp2(t),Dp3(t),diff(p1(t), t),diff(p2(t), t),diff(p2(t), t),diff(Dp1(t), t),diff(Dp2(t), t),diff(Dp3(t), t),LL1,LL2,LL3],[a,b,c,d,e,f,d,e,f,g,h,i,1,1,1])
t3=subs(tao_three,[p1(t),p2(t),p3(t),Dp1(t),Dp2(t),Dp3(t),diff(p1(t), t),diff(p2(t), t),diff(p2(t), t),diff(Dp1(t), t),diff(Dp2(t), t),diff(Dp3(t), t),LL1,LL2,LL3],[a,b,c,d,e,f,d,e,f,g,h,i,1,1,1]) 
end