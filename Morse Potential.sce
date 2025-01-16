clc;clear;clf;
h=197.3;m=940*10^6;D=0.755501;a=1.44;
b2=10;r0=0.131349;
rmin=0.01;rmax=10;n=100;
r=linspace(rmin,rmax,n);
d=r(2)-r(1);
V=zeros(n,n);
for i=1:n;
    rp=(r(i)-r0)/r(i);
    V(i,i)=D*(exp(-2*a*rp)-exp(-a*rp));
end
K=eye(n,n)*(-2);
for i=1:(n-1)
    K(i,i+1)=1; K(i+1,i)=1;
end
H=(-(h^2)/(2*m*d^2))*K+V;
[U,EV]=spec(H);
E=diag(EV);
disp("ground state energy(eV):",E(1));
disp("1st excited state energy(eV):",E(2));
plot(r',[abs(U(:,1)),abs(U(:,2))],'linewidth',5);
legend('Ground state','1st excited state',1);
xlabel('r');
ylabel('|wavefunction|');
