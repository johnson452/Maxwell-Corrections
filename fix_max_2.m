function M_Eq = fix_max_2(M_Eq,n,u,T,app)
% fix Maxwellian distribtion

% Compute the moments
kb = app.kb;
m = app.m0;
[n_c,u_c,T_c] = moments(M_Eq,app);
Q_prime = Q_to_Q_prime(n,u,T,kb,m);
Q = [n,u,T];

% Number of fixing iterations:
niter = 0;
tol = 1e-15;

% Inital guess
alphak = Q_to_alpha(n_c,u_c,T_c,kb,m);

%Initial difference
Q_comp = [n_c,u_c,T_c];
dQ1_diff = rel_diff(Q_comp(1),Q(1));
dQ2_diff = rel_diff(Q_comp(2),Q(2));
dQ3_diff = rel_diff(Q_comp(3),Q(3));
max_diff = max([dQ1_diff,dQ2_diff,dQ3_diff]);
fprintf("(Fix Max) Iter: %d, error dQ1: %1.16e, error dQ2: %1.16e, error dQ3: %1.16e\n",...
    niter,dQ1_diff,dQ2_diff,dQ3_diff);

% Iterate to correct the points
while (niter < 15 && max_diff > tol)

    % Run the moment correction routine:
    M_Eq = alpha_maxwellian(alphak,app.grid.v);
    [n_c,nu_c,nu_sq_c] = simple_moments(M_Eq,app);
    Q_prime_discreet = [n_c;nu_c;nu_sq_c];
    Jij = Jacobian(alphak,M_Eq,app);
    corr = inv(Jij)*(Q_prime_discreet - Q_prime);
    alpha_kp = alphak - corr;

    % Advance niter
    niter = niter + 1;

    %Step the interation alpha
    alphak = alpha_kp;

    % Compute the error
    m0 = app.m0;
    kb = app.kb;
    u_c = nu_c./n_c;
    T_c = (m0./(kb*n_c)).*(nu_sq_c - n_c.*u_c.*u_c);
    Q_comp = [n_c,u_c,T_c];

    dQ1_diff = rel_diff(Q_comp(1),Q(1));
    dQ2_diff = rel_diff(Q_comp(2),Q(2));
    dQ3_diff = rel_diff(Q_comp(3),Q(3));
    max_diff = max([dQ1_diff,dQ2_diff,dQ3_diff]);
    fprintf("(Fix Max) Iter: %d, error dQ1: %1.16e, error dQ2: %1.16e, error dQ3: %1.16e\n",...
        niter,dQ1_diff,dQ2_diff,dQ3_diff);


end

end


% Compute Jij
function [Jij] = Jacobian(alpha,g_alpha,app)

% Initialize the jacobian:
Jij = zeros(3,3);
moment_kernel = zeros(3,app.grid.Nv);

% Compute the equilibrium
moment_kernel(1,:) = (1/alpha(1)).*g_alpha;
moment_kernel(2,:) = (-(app.grid.v - alpha(3)).^2).*g_alpha;
moment_kernel(3,:) = (2*alpha(2)*(app.grid.v - alpha(3))).*g_alpha;

% For all components of the prime:
for j = 1:3
    [n_dg,nu_dg,nu_sq_dg] = simple_moments(moment_kernel(j,:),app);
    Jij(:,j) = [n_dg;nu_dg;nu_sq_dg];
end

end


% Maxwellian in terms of alpha
function [g_alpha] = alpha_maxwellian(alpha,v)

% Compute the distribution
g_alpha = alpha(1)*exp(-alpha(2)*(v - alpha(3)).^2);

end


% Compute Q to Q'
function [alpha] = Q_to_Q_prime(n,u,T,kb,m)

% Build the alpha vector
alpha = zeros(3,1);
alpha(1) = n;
alpha(2) = n*u;
alpha(3) = n*u*u + n*kb*T/m;

end


% Compute Q to alpha
function [alpha] = Q_to_alpha(n,u,T,kb,m)

% Build the alpha vector
alpha = zeros(3,1);
alpha(1) = n*sqrt(m/(2*pi*kb*T));
alpha(2) = m/(2*kb*T);
alpha(3) = u;

end


% Simple moments where we do not polynomial reconstruct in velocity space
function [n,nu,nu_sq] = simple_moments(f,app)

%sz f
sz_f = size(f);

% Grab quantities
v = app.grid.v;
dv = app.grid.dv;

% Build moment arrays
n = zeros(sz_f(1),1);
nu = zeros(sz_f(1),1);
nu_sq = zeros(sz_f(1),1);

% Iterate over the whole grid
for i = 1:sz_f(1)
    n(i) = dv*sum(f(i,:));
    nu(i) = dv*sum(v.*f(i,:));
    nu_sq(i) = dv*sum(v.*v.*f(i,:));
end

end