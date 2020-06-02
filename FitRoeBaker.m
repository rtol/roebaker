%<a href="https://www.nature.com/articles/ngeo3017">Knutti et al. (2017)</a>
%survey estimates of the equilibrium and transient climate sensitivity,
%the change in the global annual mean surface air temperature due to a
%doubling of the atmospheric concentration of carbon dioxide.
%This script takes their data and fits either a Normal distribution (if the
%reported confidence interval is symmetric) or <a
%href="https://science.sciencemag.org/content/318/5850/629.abstract">Roe-Baker
%(2007)</a> distribution (if not).
%For most studies, Knutti reports one of three central tendencies and one
%of three confidence intervals, so that there are 18 different cases.
%For the Normal distribution, mu is set equal to the mean, mode or median
%and sigma follows from the width of the confidence interval.
%For the Roe-Baker distribution, lambda, mu and sigma are found by
%minimizing the squared distance between the modelled and reported central
%tendency, upper bound, and lower bound.
%Each of the fitted distributions is plotted, and the composite ones that
%is unweighted average.
%
%Richard S.J. Tol, 2 June 2020

clear all
load Knutti
lambda0 = 0.30;
global NSteps MinCS MaxCS
NSteps = 3000;
MinCS = 0;
MaxCS = 15;

%%
for i=1:182,
    if ~isnan(Knutti.mean(i+1)) & ~isnan(Knutti.p50(i+1)) & ~isnan(Knutti.p95(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mean(i+1) Knutti.p50(i+1) Knutti.p95(i+1)];
            fun = @(param)RBMean90(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mean(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.645;
        end
    elseif ~isnan(Knutti.mode(i+1)) & ~isnan(Knutti.p50(i+1)) & ~isnan(Knutti.p95(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mode(i+1) Knutti.p50(i+1) Knutti.p95(i+1)];
            fun = @(param)RBMode90(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mode(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.645;
        end
    elseif ~isnan(Knutti.median(i+1)) & ~isnan(Knutti.p50(i+1)) & ~isnan(Knutti.p95(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.median(i+1) Knutti.p50(i+1) Knutti.p95(i+1)];
            fun = @(param)RBMedian90(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.median(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.645;
        end
    elseif ~isnan(Knutti.mean(i+1)) & ~isnan(Knutti.p25(i+1)) & ~isnan(Knutti.p975(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mean(i+1) Knutti.p25(i+1) Knutti.p975(i+1)];
            fun = @(param)RBMean95(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mean(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.96;
        end
    elseif ~isnan(Knutti.mode(i+1)) & ~isnan(Knutti.p25(i+1)) & ~isnan(Knutti.p975(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mode(i+1) Knutti.p25(i+1) Knutti.p975(i+1)];
            fun = @(param)RBMode95(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mode(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.96;
        end
    elseif ~isnan(Knutti.median(i+1)) & ~isnan(Knutti.p25(i+1)) & ~isnan(Knutti.p975(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.median(i+1) Knutti.p25(i+1) Knutti.p975(i+1)];
            fun = @(param)RBMedian95(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.median(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/1.96;
        end
    elseif ~isnan(Knutti.mean(i+1)) & ~isnan(Knutti.p17(i+1)) & ~isnan(Knutti.p73(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mean(i+1) Knutti.p17(i+1) Knutti.p73(i+1)];
            fun = @(param)RBMean67(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mean(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/0.97;
        end
    elseif ~isnan(Knutti.mode(i+1)) & ~isnan(Knutti.p17(i+1)) & ~isnan(Knutti.p73(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.mode(i+1) Knutti.p17(i+1) Knutti.p73(i+1)];
            fun = @(param)RBMode67(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.mode(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/0.97;
        end
    elseif ~isnan(Knutti.median(i+1)) & ~isnan(Knutti.p17(i+1)) & ~isnan(Knutti.p73(i+1)),
        if Knutti.symmetric(i+1)=="false"
            target = [Knutti.median(i+1) Knutti.p17(i+1) Knutti.p73(i+1)];
            fun = @(param)RBMedian67(param,target);
            param0 = [lambda0 1-lambda0/target(1) 0.1];

            paropt = fminsearch(fun,param0);
        
            KnuttiPar(i,1) = paropt(1);
            KnuttiPar(i,2) = paropt(2);
            KnuttiPar(i,3) = paropt(3);
        else
            KnuttiPar(i,1) = 0;
            KnuttiPar(i,2) = Knutti.median(i+1);
            KnuttiPar(i,3) = Knutti.distu(i+1)/0.97;
        end
    end
end

%%
composite = zeros(1,NSteps+1);
hold on
for i=1:182,
    if KnuttiPar(i,1)>0,
        [grid,RB(i,:),RBc(i,:),mean(i),mode(i),median(i),p025(i),p050(i),p167(i),p833(i),p950(i),p975(i)] = RoeBakerStats(MinCS,MaxCS,NSteps,KnuttiPar(i,1),KnuttiPar(i,2),KnuttiPar(i,3));
        composite = composite + RB(i,:);
        plot(grid,RB(i,:))
    elseif KnuttiPar(i,2)>0,
        [grid,RB(i,:),RBc(i,:),mean(i),mode(i),median(i),p025(i),p050(i),p167(i),p833(i),p950(i),p975(i)] = NormalStats(MinCS,MaxCS,NSteps,KnuttiPar(i,2),KnuttiPar(i,3));
        composite = composite + RB(i,:);
        plot(grid,RB(i,:))
    end
end
composite = composite/sum(composite);
plot(grid,composite,'k','LineWidth',2)
ylabel('Probability density')
xlabel('Equilibrium climate sensitivity')
hold off

%%
plot(grid,composite,'k','LineWidth',2)
ylabel('Probability density')
xlabel('Equilibrium climate sensitivity')