function [u] = fct_Steady_DeterHead_3D(Wcoef,idx,ns,tp,nrp,nlp,nnfp,nip,ntp,RBC,LBC,Bd_V)
C = sparse((nrp+nlp+nnfp+ntp),(ntp));
for ith = 1:(nrp+nlp)
    C(ith,ith) = 1;
end
for ith = (nrp+nlp+1):(nrp+nlp+nnfp) % No-Flow Boundaries
    C(ith,(idx(ith,1:ns))) = (Bd_V(ith,4)).*Wcoef(1,1:ns,ith) + (Bd_V(ith,5)).*Wcoef(2,1:ns,ith) + (Bd_V(ith,6)).*Wcoef(3,1:ns,ith);
end
for ith = (1):(ntp) % Governing Equation
    dYdx = (Wcoef(1,1:ns,ith)*tp((idx(ith,1:ns)),4));
    dYdy = (Wcoef(2,1:ns,ith)*tp((idx(ith,1:ns)),4));
    dYdz = (Wcoef(3,1:ns,ith)*tp((idx(ith,1:ns)),4));
    Cia = Wcoef(4,1:ns,ith) + Wcoef(7,1:ns,ith) + Wcoef(9,1:ns,ith);
    Cib = dYdx*Wcoef(1,1:ns,ith);
    Cic = dYdy*Wcoef(2,1:ns,ith);
    Cid = dYdz*Wcoef(3,1:ns,ith);
    C(ith+nrp+nlp+nnfp,(idx(ith,1:ns))) = Cia + Cib + Cic + Cid;
end
f(1:nrp,1) = RBC; % Dirichlet Boundaries
f((nrp+1):(nrp+nlp),1) = LBC; % Dirichlet Boundaries
f((nrp+nlp+1):(nrp+nlp+nnfp),1) = 0; % No-Flow Boundaries
f((nrp+nlp+nnfp+1):(nrp+nlp+nnfp+ntp),1) = 0; % No Sink/Source

u=C\f;