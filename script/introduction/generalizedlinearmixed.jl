
using DataFrames, RCall, MixedModels


RCall.ijulia_setdevice(MIME("image/svg+xml"), width=7, height=5)


R"""
library(lattice)
suppressPackageStartupMessages(library(mlmRev))
""";


R"""
xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception, groups = livch,
    type = c("g", "smooth"), ylab = "Proportion", xlab = "Centered age",
    auto.key = list(space = "top", points = FALSE, lines = TRUE, columns = 4))
"""


const form1 = @formula use ~ 1 + age + abs2(age) + livch + urban + (1|district);
const contra = rcopy(R"Contraception");
m1 = fit!(GeneralizedLinearMixedModel(form1, contra,
    Bernoulli()), fast=true)


contra.ch = categorical(ifelse.(contra.livch .== "0", "N", "Y"))


f2 = @formula use ~ age + abs2(age) + urban + ch + (1|district);
m2 = fit(MixedModel, f2, contra, Bernoulli(), fast=true)


R"""
Contraception$ch <- factor(Contraception$livch != 0, labels = c("N","Y"))
xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception, groups = ch, type = c("g", "smooth"),
       auto.key = list(space = "top", points = FALSE, lines = TRUE, columns = 2),
       ylab = "Proportion", xlab = "Centered age")
"""


f3 = @formula use ~ 1 + age * ch + abs2(age) + urban + (1|district);
m3 = fit(MixedModel, f3, contra, Bernoulli(), fast = true)


using LinearAlgebra, Gadfly
sym3 = SymTridiagonal(zeros(3), sqrt.(1:2))
ev = eigen(sym3);
show(ev.values)
show(abs2.(ev.vectors[1,:]))


function gausshermitenorm(k)
    ev = eigen(SymTridiagonal(zeros(k), sqrt.(1:k-1)))
    ev.values, abs2.(ev.vectors[1,:])
end


gausshermitenorm(3)


using Gadfly
gh9=gausshermitenorm(9)
plot(x=gh9[1], y=gh9[2], Geom.hair, Geom.point, Guide.ylabel("Weight"), Guide.xlabel(""))


plot(x=gh9[1], y=gh9[2], Geom.hair, Geom.point, Scale.y_log2, Guide.ylabel("Weight (log scale)"), Guide.xlabel(""))


using MixedModels
GHnorm(3)


μ = 2; σ = 3; ghn3 = GHnorm(3);
sum(@. ghn3.w * abs2(μ + σ * ghn3.z))  # should be μ² + σ² = 13


const devc0 = map!(abs2, m1.devc0, m1.u[1]);  # start with uᵢ²
const devresid = m1.resp.devresid;   # n-dimensional vector of deviance residuals
const refs = first(m1.LMM.reterms).refs;  # n-dimensional vector of indices in 1:q
for (dr, i) in zip(devresid, refs)
    devc0[i] += dr
end
show(devc0)


using FreqTables
freqtable(contra, :d)'


const devc = m1.devc;
const xvals = -5.0:2.0^(-4):5.0;
const uv = vec(m1.u[1]);
const u₀ = vec(m1.u₀[1]);
results = zeros(length(devc0), length(xvals))
for (j, u) in enumerate(xvals)
    fill!(devc, abs2(u))
    fill!(uv, u)
    MixedModels.updateη!(m1)
    for (dr, i) in zip(devresid, refs)
        devc[i] += dr
    end
    copyto!(view(results, :, j), devc)
end


plot(x=xvals, y=view(results, 1, :), Geom.line, Guide.xlabel("u₁"), Guide.ylabel("Deviance contribution"))


plot(x=xvals, y=view(results, 3, :), Geom.line, Guide.xlabel("u₃"), Guide.ylabel("Deviance contribution"))


m1.u₀[1]


const s = inv.(m1.LMM.L[Block(1,1)].diag);
s'


for (j, z) in enumerate(xvals)
    @. uv = u₀ + z * s
    MixedModels.updateη!(m1)
    @. devc = abs2(uv) - devc0
    for (dr, i) in zip(devresid, refs)
        devc[i] += dr
    end
    copyto!(view(results, :, j), devc)
end


plot(x=xvals,y=view(results, 1, :),Geom.line,Guide.xlabel("Scaled and shifted u₁"),Guide.ylabel("Shifted deviance contribution"))


plot(x=xvals,y=view(results, 3, :),Geom.line,Guide.xlabel("Scaled and shifted u₃"),Guide.ylabel("Shifted deviance contribution"))


for (j, z) in enumerate(xvals)
    @. uv = u₀ + z * s
    MixedModels.updateη!(m1)
    @. devc = abs2(uv) - devc0
    for (dr, i) in zip(devresid, refs)
        devc[i] += dr
    end
    copyto!(view(results, :, j), @. exp(-devc/2))
end


plot(x=xvals,y=view(results, 1, :),Geom.line,Guide.xlabel("Scaled and shifted u₁"),Guide.ylabel("Conditional density"))


plot(x=xvals,y=view(results, 3, :),Geom.line,Guide.xlabel("Scaled and shifted u₃"),Guide.ylabel("Conditional density"))


for (j, z) in enumerate(xvals)
    @. uv = u₀ + z * s
    MixedModels.updateη!(m1)
    @. devc = abs2(uv) - devc0
    for (dr, i) in zip(devresid, refs)
        devc[i] += dr
    end
    copyto!(view(results, :, j), @. exp((abs2(z) - devc)/2))
end


plot(x=xvals,y=view(results, 1, :),Geom.line,Guide.xlabel("Scaled and shifted u₁"),Guide.ylabel("Kernel ratio"))


plot(x=xvals,y=view(results, 3, :),Geom.line,Guide.xlabel("Scaled and shifted u₃"),Guide.ylabel("Kernel ratio"))


using MixedModelsTutorials
MixedModelsTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

