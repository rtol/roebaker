function p = RoeBakerPDF(x,lambda,mu,sigma)
%p = RoeBakerPDF(x,lambda,mu,sigma)
%returns the relative likelihood p of observation x
%according to the Roe-Baker probability density function
%given parameters lambda, mu and sigma
%
%Richard S.J. Tol, 2 June 2020

if x > 0,
    p = exp(-0.5*((1-lambda/x-mu)/sigma)^2)*lambda/sigma/sqrt(2*pi)/x/x;
else
    p = 0;

end

