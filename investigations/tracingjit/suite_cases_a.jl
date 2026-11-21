# Cases, part A: the standard library, statistics, and small utility packages.
# add!(name, package, f, R -> args; wrt, extra = [R -> other args], expect_failure)

dd(R, n) = randn(R, n, n) + n * I           # a well conditioned matrix
spd(R, n) = (A = randn(R, n, n); A' * A + n * I)
pos(R, n) = rand(R, n) .+ 0.5                 # positive vector

# ------------------------------------------------------------------------------------------------ LinearAlgebra (generic code)
const LA = "LinearAlgebra"
add!("A\\b 8x8 generic LU", LA, (A, b) -> A \ b, R -> (dd(R, 8), randn(R, 8)))
add!("det 5x5", LA, A -> det(A), R -> (dd(R, 5),))
add!("logdet SPD 6x6", LA, A -> logdet(A), R -> (spd(R, 6),))
add!("inv 5x5", LA, A -> vec(inv(A)), R -> (dd(R, 5),))
add!("A^5 4x4", LA, A -> vec(A^5), R -> (randn(R, 4, 4) / 2,))
add!("exp(A) 4x4", LA, A -> vec(exp(A)), R -> (randn(R, 4, 4) / 2,))
add!("norm 2 of 8-vector", LA, x -> norm(x), R -> (randn(R, 8),))
add!("norm 1 plus Inf of 8-vector", LA, x -> norm(x, 1) + norm(x, Inf), R -> (randn(R, 8),))
add!("opnorm 1 of 5x5", LA, A -> opnorm(A, 1), R -> (randn(R, 5, 5),))
add!("opnorm Inf of 5x5", LA, A -> opnorm(A, Inf), R -> (randn(R, 5, 5),))
add!("dot and cross", LA, (a, b) -> vcat(dot(a, b), cross(a, b)), R -> (randn(R, 3), randn(R, 3)))
add!("tr(A*B) and kron", LA, (A, B) -> vcat(tr(A * B), vec(kron(A, B))[1:6]), R -> (randn(R, 3, 3), randn(R, 3, 3)))
add!("cholesky solve 6x6", LA, (A, b) -> cholesky(Symmetric(A)) \ b, R -> (spd(R, 6), randn(R, 6)))
add!("least squares 8x4", LA, (A, b) -> A \ b, R -> (randn(R, 8, 4), randn(R, 8)))
add!("lower triangular solve 6x6", LA, (L, b) -> LowerTriangular(L) \ b, R -> (dd(R, 6), randn(R, 6)))
add!("SymTridiagonal solve 8", LA, (d, e, b) -> SymTridiagonal(d, e) \ b, R -> (4 .+ rand(R, 8), randn(R, 7), randn(R, 8)))
add!("qr: abs diag of R 6x6", LA, A -> abs.(diag(qr(A).R)), R -> (dd(R, 6),))
add!("normalize and angle", LA, (a, b) -> (u = normalize(a); v = normalize(b); acos(dot(u, v))), R -> (randn(R, 4), randn(R, 4)))
add!("matrix chain A*B*A", LA, (A, B) -> vec(A * B * A), R -> (randn(R, 4, 4), randn(R, 4, 4)))
add!("eigvals sym 3x3 (LAPACK only)", LA, A -> eigvals(Symmetric(A)), R -> (randn(R, 3, 3),); expect_failure=true)
add!("svdvals 3x3 (LAPACK only)", LA, A -> svdvals(A), R -> (randn(R, 3, 3),); expect_failure=true)
add!("opnorm 2 (needs svdvals)", LA, A -> opnorm(A, 2), R -> (randn(R, 3, 3),); expect_failure=true)

# ------------------------------------------------------------------------------------------------ Statistics (standard library)
const ST = "Statistics"
add!("mean, var, std of 20", ST, x -> [mean(x), var(x), std(x)], R -> (randn(R, 20),))
add!("median of 21 (odd)", ST, x -> median(x), R -> (randn(R, 21),))
add!("median of 20 (even)", ST, x -> median(x), R -> (randn(R, 20),))
add!("quantile 0.3 of 20", ST, x -> quantile(x, 0.3), R -> (randn(R, 20),))
add!("quantiles [.1,.5,.9] of 30", ST, x -> quantile(x, [0.1, 0.5, 0.9]), R -> (randn(R, 30),))
add!("cor and cov of two 15-vectors", ST, (x, y) -> [cor(x, y), cov(x, y)], R -> (randn(R, 15), randn(R, 15)))
add!("5th order statistic by sort", ST, x -> sort(x)[5], R -> (randn(R, 12),))
add!("max minus min", ST, x -> maximum(x) - minimum(x), R -> (randn(R, 12),))
add!("sum of abs and cumsum", ST, x -> vcat(sum(abs, x), cumsum(x)[end]), R -> (randn(R, 12),))
add!("sort descending", ST, x -> sort(x; rev=true), R -> (randn(R, 10),))
add!("partialsort k=3", ST, x -> partialsort(x, 3), R -> (randn(R, 10),))
add!("mean of abs deviations", ST, x -> mean(abs.(x .- median(x))), R -> (randn(R, 15),))
add!("clamp and sum", ST, (x, lo, hi) -> sum(clamp.(x, lo, hi)), R -> (randn(R, 12), -0.5, 0.7))

# ------------------------------------------------------------------------------------------------ NaNMath
const NM = "NaNMath"
add!("NaNMath.sum", NM, x -> NaNMath.sum(x), R -> (randn(R, 10),))
add!("NaNMath.max and min", NM, (a, b) -> [NaNMath.max(a, b), NaNMath.min(a, b)], R -> (randn(R), randn(R)))
add!("NaNMath.pow and log, sqrt", NM, (a, b) -> [NaNMath.pow(a, b), NaNMath.log(a), NaNMath.sqrt(a)], R -> (rand(R) + 0.5, randn(R)))
add!("NaNMath.maximum", NM, x -> NaNMath.maximum(x), R -> (randn(R, 10),))
add!("NaNMath.mean", NM, x -> NaNMath.mean(x), R -> (randn(R, 10),))

# ------------------------------------------------------------------------------------------------ Distances
const DS = "Distances"
add!("euclidean, sqeuclidean", DS, (x, y) -> [Distances.euclidean(x, y), Distances.sqeuclidean(x, y)], R -> (randn(R, 6), randn(R, 6)))
add!("cityblock, chebyshev", DS, (x, y) -> [Distances.cityblock(x, y), Distances.chebyshev(x, y)], R -> (randn(R, 6), randn(R, 6)))
add!("cosine and correlation distance", DS, (x, y) -> [Distances.cosine_dist(x, y), Distances.corr_dist(x, y)], R -> (randn(R, 6), randn(R, 6)))
add!("braycurtis", DS, (x, y) -> Distances.braycurtis(x, y), R -> (pos(R, 6), pos(R, 6)))
add!("KL divergence", DS, (p, q) -> Distances.kl_divergence(p ./ sum(p), q ./ sum(q)), R -> (pos(R, 5), pos(R, 5)))
add!("Jensen-Shannon", DS, (p, q) -> Distances.js_divergence(p ./ sum(p), q ./ sum(q)), R -> (pos(R, 5), pos(R, 5)))
add!("weighted euclidean", DS, (x, y, w) -> Distances.weuclidean(x, y, w), R -> (randn(R, 5), randn(R, 5), pos(R, 5)))
add!("Mahalanobis", DS, (x, y, Q) -> Distances.evaluate(Distances.Mahalanobis(Q * Q' + I), x, y), R -> (randn(R, 3), randn(R, 3), randn(R, 3, 3)); wrt=[1, 2])
add!("rmsd and meanad", DS, (x, y) -> [Distances.rmsd(x, y), Distances.meanad(x, y)], R -> (randn(R, 6), randn(R, 6)))
add!("haversine", DS, (a, b) -> Distances.Haversine(6371.0)(a, b), R -> ([10.0 * rand(R), 40.0 + rand(R)], [11.0 + rand(R), 41.0 + rand(R)]))
add!("pairwise euclidean 4 points", DS, X -> vec(Distances.pairwise(Distances.Euclidean(), X; dims=2)), R -> (randn(R, 2, 4),))

# ------------------------------------------------------------------------------------------------ StatsBase
const SB = "StatsBase"
add!("weighted mean and var", SB, (x, w) -> (ww = StatsBase.weights(w); [mean(x, ww), var(x, ww; corrected=false)]), R -> (randn(R, 10), pos(R, 10)))
add!("skewness and kurtosis", SB, x -> [StatsBase.skewness(x), StatsBase.kurtosis(x)], R -> (randn(R, 20),))
add!("mad and iqr", SB, x -> [StatsBase.mad(x), StatsBase.iqr(x)], R -> (randn(R, 20),))
add!("zscore", SB, x -> StatsBase.zscore(x), R -> (randn(R, 8),))
add!("geomean and harmmean", SB, x -> [StatsBase.geomean(x), StatsBase.harmmean(x)], R -> (pos(R, 8),))
add!("percentile 30", SB, x -> StatsBase.percentile(x, 30), R -> (randn(R, 20),))
add!("entropy and crossentropy", SB, (p, q) -> [StatsBase.entropy(p ./ sum(p)), StatsBase.crossentropy(p ./ sum(p), q ./ sum(q))], R -> (pos(R, 5), pos(R, 5)))
add!("mean_and_var", SB, x -> collect(StatsBase.mean_and_var(x)), R -> (randn(R, 10),))
add!("third moment", SB, x -> StatsBase.moment(x, 3), R -> (randn(R, 12),))

# ------------------------------------------------------------------------------------------------ Polynomials
const PL = "Polynomials"
add!("polynomial value, degree 5", PL, (c, x) -> Polynomials.Polynomial(c)(x), R -> (randn(R, 6), 0.7))
add!("polynomial derivative at x", PL, (c, x) -> Polynomials.derivative(Polynomials.Polynomial(c))(x), R -> (randn(R, 6), 0.7))
add!("polynomial integral at 1", PL, c -> Polynomials.integrate(Polynomials.Polynomial(c))(1.0), R -> (randn(R, 6),))
add!("product of two polynomials", PL, (a, b, x) -> (Polynomials.Polynomial(a) * Polynomials.Polynomial(b))(x), R -> (randn(R, 4), randn(R, 3), 0.4))
add!("least squares polynomial fit", PL, (ys, x0) -> Polynomials.fit(collect(0.0:7.0), ys, 2)(x0), R -> (randn(R, 8), 2.5))
add!("Chebyshev series value", PL, (c, x) -> Polynomials.ChebyshevT(c)(x), R -> (randn(R, 6), 0.3))
add!("polynomial roots (companion eigvals)", PL, c -> real.(Polynomials.roots(Polynomials.Polynomial(vcat(c, 1.0)))), R -> (randn(R, 3),); expect_failure=true)
