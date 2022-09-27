using ApproxOperator, YAML, LinearAlgebra, Printf

config = YAML.load_file("./yml/bar.yml")
elements, nodes = importmsh("./msh/bar.msh",config)

nₚ = getnₚ(elements["Ω"])
nₑ = nₚ - 1

set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

coefficient = (:EA=>1.0,:η=>1e-6,:k=>1.0,:l=>0.01)
ops = [Operator(:∫v²uₓuₓdx,coefficient...),
       Operator(:∫vₓvₓvvdx,coefficient...),
       Operator(:UPDATE_PFM_1D,coefficient...),
       Operator(:∫vgdΓ,coefficient...,:α=>1e11)]

prescribe!(elements["Γᵍ"],:g=>(x,y,z)->4.0*x)
prescribe!(elements["Ω"],:ℋ =>(x,y,z)->0.0)

u = zeros(nₚ)
v = ones(nₚ)
push!(nodes,:u=>u)
push!(nodes,:v=>v)
Δu = zeros(nₚ)
global normΔ = 1.0
tol = 1e-13
maxiter = 1000
global iter = 0
while normΔ > tol && iter ≤ maxiter
    global iter += 1
    # elasticity
    k = zeros(nₚ,nₚ)
    f = zeros(nₚ)
    ops[1](elements["Ω"],k,f)
    ops[4](elements["Γᵍ"],k,f)
    d = k\f
    normΔu = norm(u .- d)
    # println(normΔu )
    u .= d

    # phase field
    k = zeros(nₚ,nₚ)
    f = zeros(nₚ)
    ops[2](elements["Ω"],k,f)
    d = k\f
    normΔv = norm(v .- d)
    v .= k\f

    # update variables
    global normΔ = normΔu + normΔv 

    @printf("iter = %3i, normΔ = %10.2e\n", iter, normΔ)
end 
# ops[3](elements["Ω"])


# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# ops[1](elements["Ω"],k,f)
# ops[4](elements["Γᵍ"],k,f)
# u .= k\f

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# ops[2](elements["Ω"],k,f)
# d = k\f
# err = k*ones(nₚ) .- f 