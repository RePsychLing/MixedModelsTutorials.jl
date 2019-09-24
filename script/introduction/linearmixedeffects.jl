
using LinearAlgebra, MixedModels, RCall
RCall.ijulia_setdevice(MIME("image/svg+xml"), width=7, height=5)


R"""
require(lme4, quietly=TRUE)
require(lattice, quietly=TRUE)
xyplot(Reaction ~ Days | Subject, sleepstudy,
       type = c("g","p","r"), layout = c(9,2),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)",
       aspect = "xy")
"""


sleepstudy = rcopy(R"lme4::sleepstudy")


f1 = @formula Reaction ~ 1 + Days + (1+Days|Subject)
m1 = fit(MixedModel, f1, sleepstudy)


first(ranef(m1))


fixef(m1)


λ = first(m1.λ)


σ² = varest(m1)


Σ = σ² * λ * λ'


Σ[2,1] / sqrt(Σ[1,1] * Σ[2,2])


b = vec(first(ranef(m1)))


Int.(m1.X)  # display as Int to reduce clutter


Int.(first(m1.reterms))


Λ = kron(I(18), first(m1.λ))


show(m1.θ)


m1.optsum


UpperTriangular(m1.L')


logdet(m1)


pwrss(m1)


abs2(first(m1.L[Block(3,3)]))


varest(m1)


logdet(m1) + dof_residual(m1)*(1 + log(2π * varest(m1)))


Symmetric(m1.A, :L)


sum(abs2, sleepstudy.Reaction)


R"""
df <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fclow <- subset(df, `(Intercept)` < 251)
fchigh <- subset(df, `(Intercept)` > 251)
cc1 <- as.data.frame(coef(fm2)$Subject)
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
ff <- fixef(fm2)
with(df,
     print(xyplot(`(Intercept)` ~ Days, aspect = 1,
                  x1 = B, y1 = A,
                  panel = function(x, y, x1, y1, subscripts, ...) {
                      panel.grid(h = -1, v = -1)
                      x1 <- x1[subscripts]
                      y1 <- y1[subscripts]
                      larrows(x, y, x1, y1, type = "closed", length = 0.1,
                              angle = 15, ...)
                      lpoints(x, y,
                              pch = trellis.par.get("superpose.symbol")$pch[2],
                              col = trellis.par.get("superpose.symbol")$col[2])
                      lpoints(x1, y1,
                              pch = trellis.par.get("superpose.symbol")$pch[1],
                              col = trellis.par.get("superpose.symbol")$col[1])
                      lpoints(ff[2], ff[1],
                              pch = trellis.par.get("superpose.symbol")$pch[3],
                              col = trellis.par.get("superpose.symbol")$col[3])
                      ltext(fclow[,2], fclow[,1], row.names(fclow),
                            adj = c(0.5, 1.7))
                      ltext(fchigh[,2], fchigh[,1], row.names(fchigh),
                            adj = c(0.5, -0.6))
                  },
                  key = list(space = "top", columns = 3,
                  text = list(c("Mixed model", "Within-group", "Population")),
                  points = list(col = trellis.par.get("superpose.symbol")$col[1:3],
                  pch = trellis.par.get("superpose.symbol")$pch[1:3]))
                  )))
""";


R"""
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(9,2), type = c("g", "p", "r"),
             coef.list = df[,3:4],
             panel = function(..., coef.list) {
                 panel.xyplot(...)
                 panel.abline(as.numeric(coef.list[packet.number(),]),
                              col.line = trellis.par.get("superpose.line")$col[2],
                              lty = trellis.par.get("superpose.line")$lty[2]
                              )
                 panel.abline(fixef(fm2),
                              col.line = trellis.par.get("superpose.line")$col[4],
                              lty = trellis.par.get("superpose.line")$lty[4]
                              )
             },
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)",
             key = list(space = "top", columns = 3,
             text = list(c("Within-subject", "Mixed model", "Population")),
             lines = list(col = trellis.par.get("superpose.line")$col[c(2:1,4)],
             lty = trellis.par.get("superpose.line")$lty[c(2:1,4)]))))
""";


using MixedModelsTutorials
MixedModelsTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

