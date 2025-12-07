clear all; clc
    
mpc=pglib_opf_case30_as; type=mpc.bus(:,2);
P=-mpc.bus(:,3); P(mpc.gen(:,1))=P(mpc.gen(:,1))+mpc.gen(:,2);
Q=mpc.bus(:,4)*0; Q(mpc.bus(:,2)==1,1)=-mpc.bus(mpc.bus(:,2)==1,4);
for i=1:1:size(mpc.gen,1)
if mpc.bus(mpc.gen(i,1),2)==1
    Q(mpc.gen(i,1))=Q(mpc.gen(i,1))+mpc.gen(i,3);
    end
end
Ym=makeYbus(mpc); Ym=full(Ym); V=mpc.bus(:,8); Vdegs=deg2rad(mpc.bus(:,9));

tol = 1e-5; max_iter = 10;
[V,Vdegs, P, Q] = NTRaph(V, Vdegs, type, P/mpc.baseMVA, Q/mpc.baseMVA, full(Ym), tol, max_iter);
