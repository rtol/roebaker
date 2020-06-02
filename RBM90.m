function R2 = RBM90(par,target)

[grid,RB,RBc,mean,mode,median,p025,p050,p167,p833,p950,p975] = RoeBakerStats(0,10,par(1),par(2),par(3));

R2 = (mean-target(1))^2 + (p050-target(2))^2 + (p950-target(3))^2;

end

