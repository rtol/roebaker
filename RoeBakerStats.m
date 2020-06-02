function [grid, PDF, CDF, mean, mode, median, p025, p050, p167, p833, p950, p975] = RoeBakerStats(mn,mx,steps,lambda,mu,sigma)
%[grid, PDF, CDF, mean, mode, median, p025, p050, p167, p833, p950, p975] = RoeBakerStats(mn,mx,steps,lambda,mu,sigma)
%returns statistics of the Roe-Baker probability density function
%for parameters lambda, mu and sigma, for the range mn-mx in steps steps.
%grid has the values of x, pN are the percentiles
%
%Richard S.J. Tol, 2 June 2020

for i = 1:steps+1,
    grid(i) = mn + (i-1)*(mx-mn)/steps;
    PDF(i) = RoeBakerPDF(grid(i),lambda,mu,sigma);
end

PDF = PDF/sum(PDF);

CDF = PDF;
mean = grid(1)*PDF(1);
for i = 2:steps+1,
    CDF(i) = CDF(i-1) + PDF(i);
    mean = mean + grid(i)*PDF(i);
end 

mode = grid(PDF==max(PDF));
aux = (CDF-0.5).^2;
median = grid(aux==min(aux));
if size(median,2) > 1,
    median = (mx-mn)/2;
end
aux = (CDF-0.025).^2;
p025 = grid(aux==min(aux));
if size(p025,2) > 1,
    p025 = mn;
end
aux = (CDF-0.050).^2;
p050 = grid(aux==min(aux));
if size(p050,2) > 1,
    p050 = mn;
end
aux = (CDF-0.167).^2;
p167 = grid(aux==min(aux));
if size(p167,2) > 1,
    p167 = mn;
end
aux = (CDF-0.833).^2;
p833 = grid(aux==min(aux));
if size(p833,2) > 1,
    p833 = mx;
end
aux = (CDF-0.950).^2;
p950 = grid(aux==min(aux));
if size(p950,2) > 1,
    p950 = mx;
end
aux = (CDF-0.975).^2;
p975 = grid(aux==min(aux));
if size(p975,2) > 1,
    p975 = mx;
end

end