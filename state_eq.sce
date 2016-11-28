
g= 9.80665

Ixx = 0.078
Iyy = 0.078
Izz = 0.064
Ixz = 0.001
mass = 0.806
l = 0.2

den = (Ixx*Izz - Ixz^2)
ApL = Izz/den
ApN = Ixz/den
AqM = 1/Iyy
ArL = Ixz/den
ArN = Ixx/den

K1 = 4.321*g
K2 = 3.426*g
K3 = 3.422*g
K4 = 3.854*g

KT1 = 0.06329
KT2 = 0.04349
KT3 = 0.02104
KT4 = 0.05735


B31 = -K1/mass
B32 = -K2/mass
B33 = -K3/mass
B34 = -K4/mass

B41 = -ApL*l*K1 + ApN*KT1
B42 =  ApL*l*K2 + ApN*KT2
B43 = -ApN*KT3 
B44 = -ApN*KT4

B53 =  AqM*l*K3
B54 = -AqM*l*K4

B61 = -ArL*l*K1 + ArN*KT1
B62 =  ArL*l*K2 + ArN*KT2
B63 = -ArN*KT3
B64 = -ArN*KT4

A=[ 0 0 0 0 0 0 0 -g 0;
    0 0 0 0 0 0 g  0 0;
    0 0 0 0 0 0 0  0 0;
    0 0 0 0 0 0 0  0 0;
    0 0 0 0 0 0 0  0 0;
    0 0 0 0 0 0 0  0 0;
    0 0 0 1 0 0 0  0 0;
    0 0 0 0 1 0 0  0 0;
    0 0 0 0 0 1 0  0 0 ]
    
B=[ 0 0 0 0;
    0 0 0 0;
    B31 B32 B33 B34;
    B41 B42 B43 B44;
    0   0   B53 B54;
    B61 B62 B63 B64;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0 ]
    
Csim = eye(9,9)
C=[ 1 1 1 1 1 1 1 1 1;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0]

D = zeros(9,4)

Q=diag([1,1,1,1,1,1,100,100,1]);R=diag([1,1,1,1]);     //Usual notations x'Qx + u'Ru
Big=sysdiag(Q,R);    //Now we calculate C1 and D12
[w,wp]=fullrf(Big);C1=wp(:,1:9);D12=wp(:,10:$);   //[C1,D12]'*[C1,D12]=Big
P=syslin('c',A,B,C1,D12);    //The plant (continuous-time)
[K,X]=lqr(P)
disp(spec(A+B*K))    //check stability
disp(norm(A'*X+X*A-X*B*inv(R)*B'*X+Q,1))  //Riccati check


L=ppol(A',C',[-1 -1 -1 -1 -1 -1 -1 -1 -1])
