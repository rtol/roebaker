function R2 = RBMedian95(par,target)

global NSteps MinCS MaxCS

[grid,RB,RBc,mean,mode,median,p025,p050,p167,p833,p950,p975] = RoeBakerStats(MinCS,MaxCS,NSteps,par(1),par(2),par(3));

R2 = (median-target(1))^2 + (p025-target(2))^2 + (p975-target(3))^2;

end


