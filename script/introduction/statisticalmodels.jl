
using DataFrames, Distributions, GLM, NLreg, RCall
using Statistics, StatsBase, StatsModels, Tables


Formaldehyde = rcopy(R"Formaldehyde")


typeof(Formaldehyde)


Formaldehyde.carb


propertynames(Formaldehyde)


getproperty(Formaldehyde, :optden)


DataFrames.describe(Formaldehyde)


R"""
suppressPackageStartupMessages(library(ggplot2))
(p <- ggplot(Formaldehyde, aes(carb, optden)) + geom_point() +
      xlab("Carbohydrate concentration") + ylab("Optical Density"))
"""


f1 = @formula(optden ~ 1 + carb);
m1 = fit(LinearModel, f1, Formaldehyde)


y, X1 = modelcols(apply_schema(f1, schema(Formaldehyde)), Formaldehyde);
X1


r²(m1)   # to type the superscript 2 type \^2<tab> or use the name r2


R"""
qplot(x=$(predict(m1)), y=$(residuals(m1)), geom="point",
      ylab="Residuals", xlab="Fitted values")
"""


R"p + geom_abline(intercept=$(coef(m1)[1]), slope=$(coef(m1)[2]))"


f2 = @formula(optden ~ 1 + carb + abs2(carb))
m2 = fit(LinearModel, f2, Formaldehyde)


y, X2 = modelcols(apply_schema(f2, schema(Formaldehyde)), Formaldehyde)
X2


rss1 = sum(abs2, residuals(m1))


rss2 = sum(abs2, residuals(m2))


ftest(m2.model, m1.model)   # add a delegate method to StatsModels


InsectSprays = rcopy(R"datasets::InsectSprays")


R"""
(p2 <- ggplot(InsectSprays, aes(spray, count)) + geom_violin() +
       geom_point() + xlab("Spray") +
       ylab("Count of insects in a standard area") + coord_flip())
"""


R"p2 + scale_y_sqrt()"


f3 = @formula sqrt(count) ~ 1 + spray;
m3 = fit(LinearModel, f3, InsectSprays)


Int.(m3.mm.m)   # show Int values to emphasize differences


f3a = @formula sqrt(count) ~ 0 + spray
m3a = fit(LinearModel, f3a, InsectSprays)


Int.(m3a.mm.m)


m4 = fit(LinearModel, @formula(sqrt(count) ~ 1), InsectSprays)
ftest(m3.model, m4.model)


ftest(m3a.model, m4.model)


by(InsectSprays, :spray, mn = :count => x -> mean(sqrt.(x)))


f5 = @formula count ~ 1 + spray
m5 = fit(GeneralizedLinearModel, f5, InsectSprays, Poisson())


deviance(m5)


f6 = @formula count ~ 1
m6 = fit(GeneralizedLinearModel, f6, InsectSprays, Poisson())


deviance(m6)


dof(m5)


ccdf(Chisq(dof(m5) - dof(m6)), deviance(m6) - deviance(m5))


Theoph = rcopy(R"Theoph")


R"""
(tpl <- ggplot(Theoph, aes(Time, conc)) + geom_point() +
    geom_line() + facet_wrap(~ Subject) +
    xlab("Time since drug adminstration (hr)") +
    ylab("Serum concentration (mg/L)"))
"""


function sdOral1C(p, d)
    k  = exp(p.lk)   # transform from logarithm to rate constant
    ka = exp(p.lka)
    V  = exp(p.lV)
    t  = d.Time
    (d.Dose/V) * ((ka - k)/ka) * (exp(-t*k) - exp(-t*ka))
end


subj6col = first(groupby(Theoph, :Subject))


subj6row = Tables.rowtable(subj6col)


pars = (lk = -2.5, lka = 0.5, lV = -1.0);


pred0 = [sdOral1C(pars, d) for d in subj6row]


m7 = fit(NLregModel, subj6row, :conc, sdOral1C, pars)


using MixedModelsTutorials
MixedModelsTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

