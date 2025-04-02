function [M_Eq,error] = fix_max(M_Eq,n,u,T,app)


%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT WORKING %
fprintf("NOT WORKING\n");
%%%%%%%%%%%%%%%%%%%%%%%%%%


% fix Maxwellian distribtion

% Compute the moments
[n_c,u_c,T_c] = moments(M_Eq,app);

% Start the error array
error = [];

% Compute differences
n_diff = rel_diff(n_c,n);
u_diff = rel_diff(u_c,u);
T_diff = rel_diff(T_c,T);
max_diff = max([n_diff,u_diff,T_diff]);

% Number of fixing iterations:
niter = 0;
tol = 1e-15;

% Intial vector
M0 = [n,u,T];
Mk = [n_c,u_c,T_c];
dMk = [0,0,0];

% Print intial correction:
fprintf("(Fix Max) Iter: %d, error n: %1.3e, error u: %1.3e, error T: %1.3e\n",...
    niter,n_diff,u_diff,T_diff);
error(:,1) = [n_diff;u_diff;T_diff];

% Iterate to correct the points
while (niter < 15 && max_diff > tol)

    % Run the moment correction routine:
    max_diff = max([n_diff,u_diff,T_diff]);
    ddMk = M0 - Mk;
    dMkp = dMk + ddMk;
    Mk = M0 + dMkp;

    % Compute the new distribution
    M_Eq = maxwellian(Mk(1),Mk(2),Mk(3),app.grid_struct.v,app);

    % Recompute moments
    [n_c,u_c,T_c] = moments(M_Eq,app);
    Mk = [n_c,u_c,T_c];
    dMk = dMkp;

    % Compute the error
    n_diff = rel_diff(Mk(1),n);
    u_diff = rel_diff(Mk(2),u);
    T_diff = rel_diff(Mk(3),T);
    fprintf("(Fix Max) Iter: %d, error n: %1.3e, error u: %1.3e, error T: %1.3e\n",...
        niter,n_diff,u_diff,T_diff);

    % Build the error array
    error(:,niter+2) = [n_diff;u_diff;T_diff];

    % Advance niter
    niter = niter + 1;
end

end

