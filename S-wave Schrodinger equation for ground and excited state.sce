// S wave equation for ground and 1st excited state of H atom
clc; clear; clf;
mcsq=0.511*10^6
hbarc=1973; n=900;
e=3.795; r_min=1e-10; r_max=18;
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
    V_matrix(i,i)=-(e^2)/r(i)
end
H_matrix=T_matrix+V_matrix
[u,eigen]=spec(H_matrix)
disp("Ground state energy (in eV):", eigen(1,1))
disp("First excited state energy (in eV):", eigen(2,2))
disp("Second excited state energy (in eV):", eigen(3,3))
disp("Third excited state energy (in eV):", eigen(4,4))
u=u/sqrt(h)
plot(r,[abs(u(:,1)),abs(u(:,2)),abs(u(:,3))])
legend('Ground state','First Excited State','Second Excited State','Second Excited State')
xlabel('r');ylabel('wavefunction');xgrid();
