function [eta_vec_mat,phi_vec_array,Z_mat,X_vecs,X_mat_original,X_mat] = NicheNeutralSim(N,S,G,J,theta_vec,I_mat,beta_vec_mat)

eta_vec_mat = zeros(G,S+1);
for k = 1:G
    eta_prime_vec = betarnd(1,theta_vec(k),1,S);
    eta_vec_mat(k,1) = eta_prime_vec(1);
    for j = 2:S
        eta_vec_mat(k,j) = eta_prime_vec(j) * prod(1 - eta_prime_vec(1:j-1));
    end
    eta_vec_mat(k,S+1) = 1 - sum(eta_vec_mat(k,1:S));
end

phi_vec_array = zeros(N,G,S+1);
phi_vec_array_prime = zeros(N,G,S);
for i = 1:N
    for k = 1:G
        for j = 1:S
            phi_vec_array_prime(i,k,j) = betarnd(I_mat(i,k) * eta_vec_mat(k,j), I_mat(i,k) * (1 - sum(eta_vec_mat(k,1:j))));
        end
    end
end
phi_vec_array(:,:,1) = phi_vec_array_prime(:,:,1);
for j = 2:S
    for i = 1:N
        for k = 1:G
            phi_vec_array(i,k,j) = phi_vec_array_prime(i,k,j) * prod(1 - phi_vec_array_prime(i,k,1:j-1));
        end
    end
end
for i = 1:N
    for k = 1:G
        phi_vec_array(i,k,S+1) = 1 - sum(phi_vec_array(i,k,1:S));
    end
end

Z_mat = zeros(N,J);
for i = 1:N
    for l = 1:J
        Z_ind_vec = mnrnd(1,beta_vec_mat(i,:));
        Z_mat(i,l) = find(Z_ind_vec);
    end
end

X_vecs = zeros(N,J);
for i = 1:N
    for l = 1:J
        X_ind_vec = mnrnd(1,squeeze(phi_vec_array(i,Z_mat(i,l),1:S)/sum(phi_vec_array(i,Z_mat(i,l),1:S)))');
        X_vecs(i,l) = find(X_ind_vec);
    end
end

X_mat_original = zeros(N,S);
for i = 1:N
    for j = 1:S
        X_mat_original(i,j) = sum(X_vecs(i,:) == j);
    end
end

X_mat = X_mat_original;
X_mat (:,~any(X_mat,1)) = [];