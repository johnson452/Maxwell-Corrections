function [M_Eq,error] = fix_max_3(M_Eq,n,u,T,app)

%%% ANDERSON %%%
m = 4; %Stages
N_moment_size = 3;
k_max = 10;
total_res = 1e-14;

% Start the error array
error = [];

% Comput the inital f(x0), x0
[n_c,u_c,T_c] = moments(M_Eq,app);
M_exact = [n;u;T];
M_computed = [n_c;u_c;T_c];

% Compute the error
n_diff = rel_diff(n_c,n);
u_diff = rel_diff(u_c,u);
T_diff = rel_diff(T_c,T);
%fprintf("(Fix Max) Iter: %d, n_c: %1.8e, u_c: %1.8e, T_c: %1.8e\n",...
%    1,n_c,u_c,T_c);
fprintf("(Fix Max) Iter: %d, error n: %1.8e, error u: %1.8e, error T: %1.8e\n",...
    1,n_diff,u_diff,T_diff);
error(:,1) = [n_diff;u_diff;T_diff];

% Essentially take the first step:
% Compute vector of iterates x
x0 = M_exact;
dMkp = M_exact - M_computed;
M_computed = M_exact + dMkp;
x = [x0, M_computed];

% Compute the residuals (g = 0 if solved)
% M_exact - M_computed
g = zeros(N_moment_size,2);
for i = 1:2

    % Compute the new distribution
    % f(x) exal
    M_Eq = maxwellian(x(1,i),x(2,i),x(3,i),app.grid_struct.v,app);
    [n_c,u_c,T_c] = moments(M_Eq,app);
    comp_moms = [n_c,u_c,T_c]';

    % Compute the residuals
    g(:,i) = comp_moms - M_exact;
end

% Matrix increments in residuals (Differences)
G_k = g(:,2) - g(:,1); % Matrix of increments in residuals
X_k = x(:,2) - x(:,1); % Matrix of increments in x

% Compute the error
n_diff = rel_diff(n_c,n);
u_diff = rel_diff(u_c,u);
T_diff = rel_diff(T_c,T);
%fprintf("(Fix Max) Iter: %d, n_c: %1.8e, u_c: %1.8e, T_c: %1.8e\n",...
%    2,n_c,u_c,T_c);
fprintf("(Fix Max) Iter: %d, error n: %1.8e, error u: %1.8e, error T: %1.8e\n",...
    2,n_diff,u_diff,T_diff);
error(:,2) = [n_diff;u_diff;T_diff];

% Iterate from the second step onwards while we
% havent maxed out the k steps or the total res
k = 2;
while k < k_max && max(abs(g(:,k))) > total_res

    % Stage number to effectively use based on
    % available history
    m_k = min(k, m);

    % Solve the optimization problem by QR decomposition
    % R = qr(A) returns the upper-triangular R factor of the
    % QR decomposition A = Q*R.
    [Q,R] = qr(G_k);

    % Compute inv(R) *( transpose (Q) * g_k)
    gamma_k = R \ ( Q' * g(:,k) );

    % Compute the new iterate and new residual
    x( :, k + 1 ) = x(:, k) + g(:, k) - (X_k + G_k) * gamma_k;
    % Compute the new distribution
    % f(x_{k+1}) exal
    M_Eq = maxwellian(x(1,k+1),x(2,k+1),x(3,k+1),app.grid_struct.v,app);
    [n_c,u_c,T_c] = moments(M_Eq,app);
    comp_moms = [n_c,u_c,T_c]';
    g(:,  k + 1 ) = comp_moms - M_exact;

    % Update the increment matrices with new elements
    X_k = [X_k, x(:,  k + 1 ) - x(:,  k )];
    G_k = [G_k, g(:,  k + 1 ) - g(:,  k )];

    % resize the matrices based on history
    n_indx = size(X_k, 2);
    if n_indx > m_k
        X_k = X_k(:, n_indx - m_k + 1:end);
        G_k = G_k(:, n_indx - m_k + 1:end);
    end

    % Update the index:
    k = k + 1;

    % Compute the error
    n_diff = rel_diff(n_c,n);
    u_diff = rel_diff(u_c,u);
    T_diff = rel_diff(T_c,T);
    %fprintf("(Fix Max) Iter: %d, n_c: %1.8e, u_c: %1.8e, T_c: %1.8e\n",...
    %    k,n_c,u_c,T_c);
    fprintf("(Fix Max) Iter: %d, error n: %1.8e, error u: %1.8e, error T: %1.8e\n",...
        k,n_diff,u_diff,T_diff);
    error(:,k) = [n_diff;u_diff;T_diff];
end

end