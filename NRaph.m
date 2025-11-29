clc; clear all;

% Bus types (using a 3 bus example)
% 1 = PQ, 2 = PV, 3 = Slack. Assume unknowns to be 0 except voltages which are assumed to be 1.
% PQ : V and Degrees are unknown. (Load Bus)
% PV : Q and Degrees are unknown. (Generation Bus)
% Slack : Set initial guess (V,deg) = (1,0)

type = input('Input type vector:\n')

% ---------------- Ybus ---------------------
Ym = input('Input Admittance Matrix:\n')
% -------------- Specified powers (injections) --------------
nbus = length(type);
% Note: Using (:) ensures column vectors
V = input('Input Voltage vector:\n');
P = input('Input Power vector:\n');
Q = input('Input Q vector:\n');
Vdegs = input('Input Bus Phase vector:\n');

%% Identify PV+PQ, PQ, PV, Slack
PvPq = find(type~=3);
PQb = find(type==1);
PVb = find(type==2);
slack = find(type==3);

g = real(Ym); b = imag(Ym); 
G = g; B = b;

max_iter = 20;
tol = 1e-6;

% Pre-allocate Jacobian sub-matrices for speed/safety
n_pvpq = length(PvPq);
n_pq = length(PQb);

tic
for iter = 1:max_iter
    
    % Reset calculation vectors every iteration
    Pcalc = zeros(nbus,1);
    Qcalc = zeros(nbus,1);

    for k = 1:nbus
        for m = 1:nbus
            theta_km = Vdegs(k) - Vdegs(m);
            Pcalc(k) = Pcalc(k) + V(k)*V(m)*( g(k,m)*cos(theta_km) + b(k,m)*sin(theta_km) );
            Qcalc(k) = Qcalc(k) + V(k)*V(m)*( g(k,m)*sin(theta_km) - b(k,m)*cos(theta_km) );
        end
    end

    dP = P(PvPq) - Pcalc(PvPq); % operate on both bus types
    dQ = Q(PQb) - Qcalc(PQb); % operate only on PQ bus

    F = [dP; dQ];

    if max(abs(F)) < tol
        fprintf('Converged in %d iterations.\n', iter-1);
        break;
    end

    %--- Building Jacobian ---
    % Saves time to do it this way. Build the full square matrix and then build the specific matrices using the proper dimensions. 
    for i = 1:nbus
        for k = 1:nbus
            del = Vdegs(i) - Vdegs(k);

            if i == k
                Hfull(i,i) =  Qcalc(i) + B(i,i)*V(i)^2;
                Nfull(i,i) = -Pcalc(i)/V(i) - G(i,i)*V(i);
                Jfull(i,i) = -Pcalc(i) + G(i,i)*V(i)^2;
                Lfull(i,i) = -Qcalc(i)/V(i) + B(i,i)*V(i)^2;
            else
                Hfull(i,k) = -V(i)*V(k)*( G(i,k)*sin(del) - B(i,k)*cos(del));
                Nfull(i,k) = -V(i)*( G(i,k)*cos(del) + B(i,k)*sin(del) );
                Jfull(i,k) =  V(i)*V(k)*( G(i,k)*cos(del) + B(i,k)*sin(del) );
                Lfull(i,k) = -V(i)*( G(i,k)*sin(del) - B(i,k)*cos(del) );
            end
        end
    end
    % (n-1)x(n-1):H (n-1)x(n-m-1):N (n-m-1)x(n-1):J (n-m-1)x(n-m-1):L
    JH = Hfull(PvPq, PvPq);
    JN = Nfull(PvPq, PQb);
    JJ = Jfull(PQb, PvPq);
    JL = Lfull(PQb, PQb);
    totalJacobian = [JH, JN; 
        JJ, JL];
    %---

    Vpq = V(PQb);
    Dv  = diag(Vpq);
    %H and J correspond to angle, so do not multiply by Dv.
    Hs = JH;
    Ns = JN * Dv;
    Js = JJ;
    Ls = JL * Dv;

    Jmat = [Hs Ns; Js Ls];

    % Unknown vector u = [Δθ; ΔV/V]
    u = -Jmat \ F;          

    d_delta = u(1:length(PvPq));
    rel_dV  = u(length(PvPq)+1:end);     % ΔV / V
    dV      = Vpq .* rel_dV;             % convert to absolute ΔV

    % Update states
    Vdegs(PvPq) = Vdegs(PvPq) + d_delta;
    V(PQb)      = V(PQb)      + dV;
end
toc

% Final Calculation
Pcalc = zeros(nbus,1);
Qcalc = zeros(nbus,1);

for i = 1:nbus
    for k = 1:nbus
        thetaik = Vdegs(i) - Vdegs(k);
        Pcalc(i) = Pcalc(i) + V(i)*V(k)*( G(i,k)*cos(thetaik) + B(i,k)*sin(thetaik) );
        Qcalc(i) = Qcalc(i) + V(i)*V(k)*( G(i,k)*sin(thetaik) - B(i,k)*cos(thetaik) );
    end
end
Pslack = Pcalc(slack);
Qslack = Qcalc(slack);

fprintf('\nBus Voltages:\n');
for k = 1:nbus
    fprintf('Bus %d: |V| = %.4f pu, angle = %.4f deg\n', k, V(k), rad2deg(Vdegs(k)));
end

fprintf('\nSlack bus power injection (bus %d): P = %.4f pu, Q = %.4f pu\n', slack, Pslack, Qslack);

