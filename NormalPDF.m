function p = NormalPDF(x,mu,sigma)
%p = NormalPDF(x,lambda,mu,sigma)
%returns the relative likelihood p of observation x
%according to the Normal probability density function
%given parameters mu and sigma
%
%Richard S.J. Tol, 2 June 2020

p = exp(-0.5*((x-mu)/sigma)^2)/sigma/sqrt(2*pi);

end