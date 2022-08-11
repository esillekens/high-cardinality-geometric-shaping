function s = tr_step(g,Bk,D)
% Cauchy point
%s0 = -D/sqrt(sum(abs(g).^2))*g;
%if g'*Bk*g <= 0
%    tau = 1;
%else
%    tau = min(sum(abs(g).^3).^(1/3)./(D*g'*Bk*g),1);
%end
%s0 = tau*s0;

% Dogleg


% % pB = -inv(Bk)*g; %newton point
% pB = -Bk\g;
% if sqrt(sum(abs(pB).^2))<=D
%     % use the Newton point if it falls within the trust radius
%     s0 = pB;
% else
%     @(tau)deal(sum((z+tau*d).^2)-D.^2,0)
% %     pU = - (g.'*g)/(g.'*Bk*g)*g;% cauchy point
% %     optfun = @(tau) abs( sum(dogleg(tau,pU,pB).^2) - D.^2 ).^2;
% %     tau2 =golden(optfun,1,2,1e-9);
% %     if sum(dogleg(tau2,pU,pB).^2)-D.^2 > 1e-9
% %         tau = golden(optfun,0,1,1e-9);
% %         while sum(dogleg(tau,pU,pB).^2)>=D.^2
% %             tau = (1-1e-3)*tau;
% %         end
% %     else
% %         tau = tau2;
% %     end
% %     s0 = dogleg(tau,pU,pB);
% 
%     s0 = -D/sqrt(sum(abs(g).^2))*g;
%     if g'*Bk*g <= 0
%         tau = 1;
%     else
%         tau = min(sum(abs(g).^3).^(1/3)./(D*g'*Bk*g),1);
%     end
%     s0 = tau*s0;
% end

% CG  Steilhaug

e = D*1e-9;

z= 0*g;
r = g;
d = -g;

if sqrt(sum(r.^2))<e
    s0 = e*g;
% keyboard
else
    for j = 1:500
        if (d'*Bk*d)<0
            mk = @(s) g'*s +1/2*s'*Bk*s;
            tau0 = sqrt((sum(D.^2)-sum(z.^2))./sum(d.^2));
%try
            tau = fmincon(@(tau) mk(z+tau*d),tau0,[],[],[],[],[],[],@(tau)deal(sum((z+tau*d).^2)-D.^2,0),optimoptions('fmincon','Display','none'));
%catch
%tau =0;
%end
            s0 = z+tau*d;
            if (sum(s0.^2)-D^2)>1e-6
                %keyboard
            end
            break
        end
        alpha = r'*r /(d'*Bk*d);
        zp1 =z+alpha*d;
        if sqrt(sum(zp1.^2))>D
            tau = golden(@(tau) (sum((z+tau*d).^2)-D.^2)^2,0,alpha,1e-6);
            s0 = z+tau*d;
            if (sum(s0.^2)-D^2)>1e-6
                %keyboard
            end
            break
        end
        rp1 = r+ alpha*Bk*d;
        if sqrt(sum( rp1.^2)) < e
            s0 = zp1;
            if sum(s0.^2)>D^2
                %keyboard
            end
            break;
        end
        beta = rp1'*rp1 / (r'*r);
        d = -rp1 + beta*d;
        
        r = rp1;
        z = zp1;
    end
    if (j==500) && all(isfinite(d))
        mk = @(s) g'*s +1/2*s'*Bk*s;
            tau0 = sqrt((sum(D.^2)-sum(z.^2))./sum(d.^2));
            if ~isfinite(tau0)
                keyboard
            end
% try
        tau = fmincon(@(tau) mk(zp1+tau*d),tau0,[],[],[],[],[],[],@(tau)deal(sum((z+tau*d).^2)-D.^2,0),optimoptions('fmincon','Display','none'));
% catch
% tau =0;
% end
        s0 = zp1+tau*d;
        if (sum(s0.^2)-D^2)>1e-6
            %keyboard
        end
    else
        % Cauchy point
        s0 = -D/sqrt(sum(abs(g).^2))*g;
    end
end

s = s0;
end

function s = dogleg(tau,pU,pB)
s = zeros(size(pU,1),size(tau,2));   
if any(tau<=1)
s(:,tau<=1) = tau(tau<=1).*pU;
end
if any(tau>1)
s(:,tau>1) = pU + (tau(tau>1)-1).*(pB-pU);
end
end

function x = golden(fun,a,b,tol)
r = (sqrt(5)-1)/2;

c = b - r*(b-a);
d = a + r*(b-a);

for iter = 1:100
    if fun(c) < fun(d)
        b = d;
    else
        a = c;
    end
        c = b - r*(b-a);
    d = a + r*(b-a);
    if abs(b-a)<tol
        break
    end
end

x = (b+a)/2;
end

function t = inters(z, d, D)
% Solve the scalar quadratic equation ||z + t d|| == D.
a = d'*d;
b = 2*z'*d;
c = z'*z - D.^2;

sqrt_discr = sqrt(b.^2 - 4*a*c);

% t = sort([(-b + sqrt_discr)/2*a,(-b - sqrt_discr)/2*a]);
aux = b +sign(b)*sqrt_discr;
t = sort([-aux/(2*a), -2*c/aux]);
end

