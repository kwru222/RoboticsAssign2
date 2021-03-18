syms p1(t) p2(t) p3(t) Dp1(t) Dp2(t) Dp3(t) LL1 LL2 LL3 LLm_one LLm_two LLm_three;

%%%%%%%.............Taking the user inputs...............
prompt1 = 'Input the value of p1(radian): ';
p1 = input(prompt1)

prompt2 = 'Input the value of p2(radian): ';
p2 = input(prompt2)

prompt3 = 'Input the value of p3(radian): ';
p3 = input(prompt3)


prompt4 = 'Input the value of dp1: ';
dp1 = input(prompt4)


prompt5 = 'Input the value of dp2: ';
dp2 = input(prompt5)

prompt6 = 'Input the value of dp3: ';
dp3 = input(prompt6)

prompt7 = 'Input the value of ddp1: ';
ddp1 = input(prompt7)

prompt8 = 'Input the value of ddp2: ';
ddp2 = input(prompt8)

prompt9 = 'Input the value of ddp3: ';
ddp3 = input(prompt9)


[t1 t2 t3]=inverse_dynamics_one(p1,p2,p3,dp1,dp2,dp3,ddp1,ddp2,ddp3)
