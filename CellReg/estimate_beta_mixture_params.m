function [p,q] = estimate_beta_mixture_params(assignments,data,maximal_distance)
% Estimates the parameters of an unstandard beta distribution
% (one that ranges from [0,maximal_distance] using maximum
% likelihood estimation via Newton-Raphson.
%Note this is also used to fit the beta distribution in the spatial
%Correlation model. This is done by assigning maximal_distance to 1.


%Values =maximal_distance will cause log(0)..bad.
data(data==maximal_distance)=0.9999*maximal_distance;
g1=sum(assignments.*log(data/maximal_distance))/sum(assignments);
g2=sum(assignments.*log((maximal_distance-data)/maximal_distance))/sum(assignments);
%Initial guess via moment-matching

sample_mean = sum(assignments.*data)/sum(assignments);
sample_var = sum(assignments.*(data-sample_mean).^2)/sum(assignments);

xbar = sample_mean/maximal_distance;
ssq= sample_var/(maximal_distance)^2;

p = xbar*((xbar*(1-xbar)/ssq)-1);
q = (1-xbar)*((xbar*(1-xbar)/ssq)-1);

x = [p;q];
%Newton-Raphson
for i=1:100
    grad = [psi(p)-psi(p+q)-g1;psi(q)-psi(p+q)-g2];
    
    hess = [psi(1,p)-psi(1,p+q), -psi(1,p+q);-psi(1,p+q),psi(1,q)-psi(1,p+q)];
    x= x-hess\grad;
    p = x(1);
    q=x(2);
end

end

