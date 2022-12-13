 function [M,X,b,l] = bwd_proj(M,X,b,u_ij,V_ij,b_0,gamma,sgn,l)

p = trace(V_ij*V_ij'*X*X) + u_ij'*M*M*u_ij;
q = trace(V_ij*V_ij'*X) + u_ij'*M*u_ij;

c_2 = p/(b_0*gamma^2) - q^2;
c_1 = sgn*2*(p/(gamma*sqrt(b_0)*sqrt(b))  + q);
c_0 = p/b - 1;

desc = sqrt(c_1^2 - 4*c_2*c_0);
l_(1) = (-c_1+real(desc))/(2*c_2);
l_(2) = (-c_1-real(desc))/(2*c_2);

if sgn == 1
    k = l_ < 1/q & l_ > -gamma*sqrt(b_0)/sqrt(b);
else
    k = l_ > -1/q & l_ < gamma*sqrt(b_0)/sqrt(b);
end

l = min([l_(k),l]);
X = X + sgn*l*X*(V_ij*V_ij')*X/(1-sgn*l*q);
M = M + sgn*l*M*u_ij*u_ij'*M/(1-sgn*l*q);
b = b_0*b./(sqrt(b_0) + sgn*sqrt(b).*l/gamma).^2;
    
 end