function [Z_vecs_mat,eta_vecs_array,T_mat_array,I_mat_array,theta_vec_mat] = NicheNeutralGibbsFast(X_mat,G,N_its)

% This function implements the Gibbs sampling algorithm for the niche neutral model.

% Variable names:

% X_mat = N by S counts matrix;
% G = the number of niches;
% N_its = the number of MCMC samples.

% Storage for the outputs of interest:

% Z_vecs_mat = J_1 + ... + J_N (= J) by N_its matrix of niche indicator
% variables;
% eta_vecs_array = G by S+1 by N_its array of niche metacommunities;
% T_mat_array = N by S by G by N_its array of table numbers;
% I_mat_array = N by G ny N_its array of immigration rates;
% theta_vec_mat = G by N_its matrix of niche biodiversity parameters.

% Identifying N and S:

[N,S] = size(X_mat);

% Create a species label vector of the whole data:

J = sum(sum(X_mat));
X_vec = zeros(J,1);
count = 0;
for i = 1:N
    for j = 1:S
        X_vec(count + 1 : count + X_mat(i,j)) = j;
        count = count + X_mat(i,j);
    end
end

% Create a sample label vector of the whole data:

X_sample = zeros(J,1);
count2 = 0;
for i = 1:N
    for j = 1:S
        X_sample(count2 + 1 : count2 + X_mat(i,j)) = i;
        count2 = count2 + X_mat(i,j);
    end
end

% Specifying the hyperparameters:

alpha = 0.01;
zeta = 0.01;
tau = 0.01;
nu = 0.01;
gamma = 0.01;

% Initialise the outputs of interest and create an array with the counts of reads allocated to niches in samples:

Z_vecs_mat = zeros(J,N_its);
Z_vecs_mat(:,1) = unidrnd(G,J,1);
eta_vecs_array = zeros(G,S+1,N_its);
gamma_rnd_mat2 = gamrnd(1,1,G,S+1);
eta_vecs_array(:,:,1) = gamma_rnd_mat2 ./ repmat(sum(gamma_rnd_mat2,2),1,S+1);
T_mat_array = zeros(N,S,G,N_its);
n_array = zeros(N,S,G);
for i = 1:N
    for j = 1:S
        for k = 1:G
            n_array(i,j,k) = sum(X_sample == i & X_vec == j & Z_vecs_mat(:,1) == k);
        end
    end
end
for i = 1:N
    for j = 1:S
        for k = 1:G
            if n_array(i,j,k) == 0
                T_mat_array(i,j,k,1) = 0;
            else
                T_mat_array(i,j,k,1) = 1;
            end
        end
    end
end
log_Stirling_array = zeros(N,S,G);
I_mat_array = zeros(N,G,N_its);
I_mat_array(:,:,1) = gamrnd(10^2/100,100/10,N,G);
theta_vec_mat = zeros(G,N_its);
theta_vec_mat(:,1) = gamrnd(10^2/100,100/10,G,1);

% Precalculate Stirling numbers required:

Stirling_mat = zeros(max(max(X_mat)),max(max(X_mat)));
Log_Normalisation_vec = zeros(max(max(X_mat)),1);

for i = 1:max(max(X_mat))
[Stirling_mat(i,1:i),Log_Normalisation_vec(i)] = stirling(i);
end

% The Gibbs sampling algorithm itself:

for n_it = 2:N_its
    for l = 1:J
        i = X_sample(l);
        j = X_vec(l);
        old_k = Z_vecs_mat(l,n_it-1);
        n_array(i,j,old_k) = n_array(i,j,old_k) - 1;
        log_prob_vec = zeros(G,1);
        for k = 1:G
            log_prob_vec(k) = log(n_array(i,j,k) + (I_mat_array(i,k,n_it-1) * eta_vecs_array(k,j,n_it-1))) + log(sum(n_array(i,:,k)) + gamma) - log(sum(n_array(i,:,k)) + I_mat_array(i,k,n_it-1)) - log(sum(sum(n_array(i,:,:))) + (G*gamma));
        end
        D = max(log_prob_vec);
        log_prob_vec_adj = log_prob_vec - D;
        log_norm_const = - log(sum(exp(log_prob_vec_adj))) - D;
        norm_log_prob_vec = log_prob_vec + log_norm_const;
        prob_vec = safeexp(norm_log_prob_vec);
        prob_vec = prob_vec / sum(prob_vec);
        cdf = cumsum(prob_vec);
        u = rand;
        Z_vecs_mat(l,n_it) = sum(cdf<=u) + 1;
        n_array(i,j,Z_vecs_mat(l,n_it)) = n_array(i,j,Z_vecs_mat(l,n_it)) + 1;
    end
    for k = 1:G
        Dirichlet_params3 = [sum(T_mat_array(:,:,k,n_it-1)) theta_vec_mat(k,n_it-1)];
        gamma_rnd_vec3 = gamrnd(Dirichlet_params3,1);
        eta_vecs_array(k,:,n_it) = gamma_rnd_vec3 / sum(gamma_rnd_vec3);
    end
    for i = 1:N
        for j = 1:S
            for k = 1:G
                if n_array(i,j,k) == 0
                    T_mat_array(i,j,k,n_it) = 0;
                    log_Stirling_array(i,j,k) = 0;
                elseif n_array(i,j,k) == 1
                    T_mat_array(i,j,k,n_it) = 1;
                    log_Stirling_array(i,j,k) = 0;
                else
                    Stirling_vec = Stirling_mat(n_array(i,j,k),1:n_array(i,j,k));
                    Log_Normalisation = Log_Normalisation_vec(n_array(i,j,k));
                    log_prob_vec2 = gammaln(I_mat_array(i,k,n_it-1) * eta_vecs_array(k,j,n_it)) - gammaln(n_array(i,j,k) + (I_mat_array(i,k,n_it-1) * eta_vecs_array(k,j,n_it))) + Log_Normalisation + log(Stirling_vec') + ((1:n_array(i,j,k))' * log(I_mat_array(i,k,n_it-1) * eta_vecs_array(k,j,n_it)));
                    D2 = max(log_prob_vec2);
                    log_prob_vec_adj2 = log_prob_vec2 - D2;
                    log_norm_const2 = - log(sum(exp(log_prob_vec_adj2))) - D2;
                    norm_log_prob_vec2 = log_prob_vec2 + log_norm_const2;
                    prob_vec2 = safeexp(norm_log_prob_vec2);
                    prob_vec2 = prob_vec2 / sum(prob_vec2);
                    cdf2 = cumsum(prob_vec2);
                    u2 = rand;
                    T_mat_array(i,j,k,n_it) = sum(cdf2 <= u2) + 1;
                    log_Stirling_array(i,j,k) = Log_Normalisation + log(Stirling_vec(T_mat_array(i,j,k,n_it)));
                end
            end
        end
    end
    n_i_dot_k = squeeze(sum(n_array,2));
    w_mat = betarnd(I_mat_array(:,:,n_it-1) + 1,n_i_dot_k);
    s_mat = binornd(1,n_i_dot_k ./ (n_i_dot_k + I_mat_array(:,:,n_it-1)));
    I_mat_array(:,:,n_it) = gamrnd(tau + squeeze(sum(T_mat_array(:,:,:,n_it),2)) - s_mat,(nu - log(w_mat)).^-1);
    T_dot_dot_k = squeeze(sum(sum(T_mat_array(:,:,:,n_it))));
    xi_vec = betarnd(theta_vec_mat(:,n_it-1) + 1,T_dot_dot_k);
    rho_vec = binornd(1,T_dot_dot_k ./ (T_dot_dot_k + theta_vec_mat(:,n_it-1)));
    theta_vec_mat(:,n_it) = gamrnd(alpha + S - rho_vec,(zeta - log(xi_vec)).^-1);
end