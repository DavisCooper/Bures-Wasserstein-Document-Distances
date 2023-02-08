function [M,X,b] = bwdml(M_0, X_0, u, S, C, gamma, maxItrs, tol)

m = size(C,1);
l = zeros(m,1);
l_ = zeros(m,1);
b = C(:,4);

valid = ones(m,1);
for i=1:m
   i1 = C(i,1);
   i2 = C(i,2);
   u_ij = u(:,i1) - u(:,i2); 
   S_ij = S(:,:,i1) - S(:,:,i2);
   if norm(u_ij) + norm(S_ij)  < 10e-10
      valid(i) = 0;
   end
end
C = C(valid>0,:);
m = size(C,1);

try
    sqrt_X_0 = sqrtm(X_0);
    X = sqrt_X_0;

    sqrt_M_0 = sqrtm(M_0);
    M = sqrt_M_0;
    
    for i=1:maxItrs
        for k=1:m 
            i1 = C(k,1);
            i2 = C(k,2);
            [U,~,V] = svd(S(:,:,i2)*X^2*S(:,:,i1));
            U = U*V';
            V_ij = S(:,:,i1) - S(:,:,i2)*U; 
            u_ij = u(:,i1) - u(:,i2); 
            sgn = C(k,3);
            b_0 = C(k,4);
            
            if sgn*(u_ij'*M*M*u_ij + trace(V_ij*V_ij'*X*X)) >= sgn*b(k)
                [M,X,b(k),a] = bwd_proj(M,X,b(k),u_ij,V_ij,b_0,gamma,sgn,l(k));
                l(k) = l(k) - a;
            end
        end
        
        norm_sum = norm(l) + norm(l_);
        if norm_sum == 0
            break
        else
            conv = norm(l_ - l,1)/norm_sum;
            if conv < tol
                break
            end
            l_ = l;
            
        end
        
    end
    
    M = sqrt_M_0\M*M/sqrt_M_0;
    X = sqrt_X_0\X*X/sqrt_X_0;
catch ex
    msgText = getReport(ex)
end