using Revise
using CMPSKit
using KrylovKit
using OptimKit
using LinearAlgebra
using JLD2
using TensorOperations

χ = 4
k = 1.
μ1 = 2.0
μ2 = 2.0
c1 = 1.0
c2 = 1.0
c12 = -0.5
c21 = -0.5

alg1 = LBFGS(; verbosity = 4, maxiter = 40000, gradtol = 1e-5);
linalg = GMRES(krylovdim = 80)

Id = 1*Matrix(I,χ,χ)
KL = Constant(randn(χ^2,χ^2))
KL = 0.5*(KL-KL')
R1 = Constant(kron(Id,randn(χ,χ)))
R2 = Constant(kron(randn(χ,χ),Id))
RLs = (R1,R2)
QL = KL
for R in RLs
    mul!(QL, R', R, -1/2, 1)
end
#Put them in cMPS form
Ψ = InfiniteCMPS(QL, (R1,R2); gauge = :left)

h = k * (∂ψ[1]'*∂ψ[1] + ∂ψ[2]'*∂ψ[2]) - μ1 * ψ[1]'*ψ[1] - μ2 * ψ[2]'*ψ[2] + c1 * (ψ[1]')^2*ψ[1]^2 + c2 * (ψ[2]')^2*ψ[2]^2 + c12 * ψ[1]'*ψ[2]'*ψ[2]*ψ[1] + c21 * ψ[2]'*ψ[1]'*ψ[1]*ψ[2] 
H = ∫(h, (-Inf,+Inf))

Ψ, ρR, E, e, normgrad, numfg, history = groundstate_tensprod(H, Ψ; optalg = alg1, linalg = linalg);

@save "data_tens_prod_$(μ1)_c=$(c1)_c12=$(c12)_D=$(χ)" Ψ E history