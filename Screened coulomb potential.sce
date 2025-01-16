//solve Schrodinger equation for screened coulomb potential
clc;clear;clf;
h=1973;m=0.511e6;e=3.795;a1=3;a2=5;a3=7//define constants
rmin=1e-10;rmax=10;n=1000;//range for r is defined
r=linspace(rmin,rmax,n);//grid of r is created
d=r(2)-r(1);//cal spacing between grid points
//create potential energy matrix
V1=zeros(n,n);
for i=1:n
    V1(i,i)=(-((e^2)*exp((-r(i))/a1))/r(i));
end
V2=zeros(n,n);
for i=1:n
    V2(i,i)=(-((e^2)*exp((-r(i))/a2))/r(i));
end
V3=zeros(n,n);
for i=1:n
    V3(i,i)=(-((e^2)*exp((-r(i))/a3))/r(i));
end
//create kinetic energy matrix
K=eye(n,n)^(-2);
for i=1:(n-1)
    K(i,i+1)=1;
    K(i+1,i)=1;
end
//evaluate HAMILTONIAN matrix H
H1=(-(h^2)/(2*m*d^2))*K+V1;
H2=(-(h^2)/(2*m*d^2))*K+V2;
H3=(-(h^2)/(2*m*d^2))*K+V3;
//cal eigen values and eigen vectors of H matrix
[Y1,EV1]=spec(H1);
[Y2,EV2]=spec(H2);
[Y3,EV3]=spec(H3);
//extracting the eigen values(at diagonal)
E1=diag(EV1);
E2=diag(EV2);
E3=diag(EV3);
//display energy eigen values(at diagonal)
format(6);
//first value is not valid
disp("ground state energy :"+string(E1(2))+"eV");
disp("1st excited state energy :"+string(E1(3))+"eV");
disp("ground state energy :"+string(E2(2))+"eV");
disp("1st excited state energy :"+string(E2(3))+"eV");
disp("ground state energy :"+string(E3(2))+"eV");
disp("1st excited state energy :"+string(E3(3))+"eV");
subplot(3,1,1);
//r is a row vector so we transpose it to match the dimension of Y
plot(r',[abs(Y1(:,2)),abs(Y1(:,3))],'linewidth',3);
legend('groundstate','1st excited state',4);
xlabel('r','fontsize',5);
ylabel('wavefunction','fontsize',5);
title("a=5A",'position',[5 0.08]);
subplot(3,1,2);
plot(r',[abs(Y2(:,2)),abs(Y2(:,3))],'linewidth',3);
legend('groundstate','1st excited state',4);
xlabel('r','fontsize',5);
ylabel('wavefunction','fontsize',5);
title("a-5A",'position',[5 0.08]);
subplot(3,1,3);
plot(r',[abs(Y1(:,2)),abs(Y1(:,3))],'linewidth',3);
legend('groundstate','1st excited state',4);
xlabel('r','fontsize',5);
ylabel('wavefunction','fontsize',5);
title("a-7A",'position',[5 0.08]);
