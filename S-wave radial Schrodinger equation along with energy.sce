// Solve the s-wave radial Schrodinger equation for an atom:
clc; clear; clf;
mcsq=0.511*10^6; hbarc=1973;
e=3.795; a=700;
r_min=0; r_max=40; n=100;
h=(r_max-r_min)/n
alpha=-(hbarc^2)/(2*mcsq*h^2)
T=zeros(n-1,n-1)
for i=1:n-1
    T(i,i)=-2
    if i<n-1 then
        T(i,i+1)=1
        T(i+1,i)=1
    end
end
T_matrix=alpha*T
V_matrix=zeros(n-1,n-1)
for i=1:n-1
    r(i)=r_min+i*h
    V_matrix(i,i)=-(e^(2)/r(i))*(exp(-r(i)/a))
    V(i)=V_matrix(i,i)
end
H_matrix=T_matrix+V_matrix
[u,eigen]=spec(H_matrix)
disp("Ground state enrgy (in eV):",eigen(1,1))
disp("First excited state energy (in eV):",eigen(2,2))
disp("Second sxcited state energy (in eV):",eigen(3,3))
disp("Third excited state energy (in eV):",eigen(4,4))
u=u/sqrt(h)
subplot(2,1,2)
plot(r,u(:,1))
plot(r,u(:,2),'*')
e=legend('Ground','1st')
subplot(2,1,1)
for i=1:n-1
    plot(r(i),V(i),'*')
    xlabel("r")
    ylabel("Potential")
end
