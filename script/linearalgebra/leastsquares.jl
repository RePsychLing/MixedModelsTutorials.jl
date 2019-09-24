
A = rand(4,3)   # simulate a matrix with elements selected at random from (0,1)


B = rand(3, 4)


A'   # "lazy" transpose


A*B


using LinearAlgebra, RCall, StatsModels, Tables


Formaldehyde = rcopy(R"Formaldehyde")


R"""
library(ggplot2)
qplot(x=carb, y=optden, data=Formaldehyde, geom="point")
"""


X = hcat(ones(size(Formaldehyde, 1)), Formaldehyde.carb)


y = Formaldehyde.optden


β = X\y     # least squares estimate


r = y - X*β   #residual


X'r   # not exactly zero but very small entries


using GLM


m1 = fit(LinearModel, @formula(optden ~ 1 + carb), Formaldehyde)


f1 = @formula(optden ~ 1 + carb)
y, X = modelcols(apply_schema(f1, schema(Formaldehyde)), Formaldehyde)


X


X


xpx = X'X


ch = cholesky(xpx)


ch.U'ch.U


ch.U'ch.U ≈ xpx


v = ldiv!(ch.U', X'y)


βc = ldiv!(ch.U, copy(v))     # solution from the Cholesky factorization


βc ≈ β


ldiv!(ch, X'y)


Xy = hcat(X, y)              # augmented model matrix


cha = cholesky(Xy'Xy)        # augmented Cholesky factor


RXX = UpperTriangular(view(cha.U, 1:2, 1:2))


RXX ≈ ch.U


rXy = cha.U[1:2, end]     # creates a copy ("view" doesn't copy)


rXy ≈ v


βac = ldiv!(RXX, copy(rXy))     # least squares solution from the augmented Cholesky


βac ≈ β


abs2(cha.U[end,end]) ≈ sum(abs2, y - X*β)   # check on residual sum of squares


Formrows = Tables.rowtable(Formaldehyde)


chr = cholesky(zeros(3, 3) + I)  # initialize


fill!(chr.factors, 0);   # zero out the contents
chr


fill!(chr.factors, 0)
for r in Formrows
    lowrankupdate!(chr, [1.0, r.carb, r.optden])
end
chr


qrX = qr(X)


qrX.R'qrX.R


qrX.R'qrX.R ≈ xpx

