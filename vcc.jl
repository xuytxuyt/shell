using ApproxOperator, Tensors
import BenchmarkExample: BenchmarkExample

include("import_vcc.jl")

elements, nodes = import_vcc("msh/scordelislo_1.msh")

# cs = BenchmarkExample.cylindricalCoordinate(25.0)
cs = BenchmarkExample.sphericalCoordinate(25.0)


𝒂₁ = cs.𝒂₁
𝒂₂ = cs.𝒂₂
𝒂₃ = cs.𝒂₃
𝒂¹ = cs.𝒂¹
𝒂² = cs.𝒂²
𝒂³ = cs.𝒂³
Γ¹₁₁ = cs.Γ¹₁₁
Γ¹₁₂ = cs.Γ¹₁₂
Γ²₂₁ = cs.Γ²₁₂
Γ²₂₂ = cs.Γ²₂₂
𝐽 = cs.𝐽

# 𝒂₁(x) = Vec{3}((1.0,0.0,0.0))
# 𝒂₂(x) = Vec{3}((0.0,1.0,0.0))
# 𝒂₃(x) = Vec{3}((0.0,0.0,1.0))
# 𝒂¹(x) = Vec{3}((1.0,0.0,0.0))
# 𝒂²(x) = Vec{3}((0.0,1.0,0.0))
# 𝒂³(x) = Vec{3}((0.0,0.0,1.0))
# Γ¹₁₁(x) = 0.0
# Γ¹₁₂(x) = 0.0
# Γ²₂₁(x) = 0.0
# Γ²₂₂(x) = 0.0

# 𝑓(x) = 1.0+2.0*x[1]+3.0*x[2] + 4.0*x[1]^2 + 5.0*x[1]*x[2] + 6.0*x[2]^6
# 𝑓(x) = 1.0
𝑓(x) = x[1]
∂ₐ𝑓(x) = gradient(𝑓,x)

function ∇𝑓(x)
    grad = ∂ₐ𝑓(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return grad[1]*𝒂¹(x_) + grad[2]*𝒂²(x_)
end

# 𝑣ᵃ(x) = Vec{2}((4.0+5.0*x[1]+6.0*x[2],7.0+8.0*x[1]+9.0*x[2]))
𝑣ᵃ(x) = Vec{2}((1.0,0.0))
function 𝒗(x)
    v = 𝑣ᵃ(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return v[1]*𝒂₁(x_) + v[2]*𝒂₂(x_)
end
𝑣ᵃₐ(x) = divergence(𝑣ᵃ,x)
function div𝒗(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    ∇𝒗 = gradient(𝒗,x_)
    return 𝒂¹(x_)⋅∇𝒗[:,1] + 𝒂²(x_)⋅∇𝒗[:,2]
end

nₐ(x) = Vec{2}((1/2^0.5,1/2^0.5))
function 𝒏(x)
    n = nₐ(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return n[1]*𝒂¹(x_) + n[2]*𝒂²(x_)
end

function Γᵞᵧₐ(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return Vec{2}((Γ¹₁₁(x_)+Γ²₂₁(x_),Γ¹₁₂(x_)+Γ²₂₂(x_)))
end

x = Vec{2}((rand(),rand()))
(check1 = ∇𝑓(x)⋅𝒗(x) - ∂ₐ𝑓(x)⋅𝑣ᵃ(x);check1 ≈ 0.0) ? println("Check 1 right!: $check1") : println("Check 1 wrong: $check1")
(check2 = 𝒗(x)⋅𝒏(x) - 𝑣ᵃ(x)⋅nₐ(x);check2 ≈ 0.0) ? println("Check 2 right!: $check2") : println("Check 2 wrong: $check2")
(check3 = div𝒗(x) - 𝑣ᵃₐ(x) - Γᵞᵧₐ(x)⋅𝑣ᵃ(x);check3 ≈ 0.0) ? println("Check 3 right!: $check3") : println("Check 3 wrong: $check3")

err = 0.0
temp1 = 0.0
temp2 = 0.0
for a in elements["Ω"]
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        x_ = Vec{2}((ξ.x,ξ.y))
        global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝑤
        global temp1 += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝑤
        # println(Γᵞᵧₐ(x_))
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*𝑣ᵃₐ(x_))*𝑤
        # global temp += 𝑓(x_)*Γᵞᵧₐ(x_)⋅𝑣ᵃ(x_)*𝑤
        # println(temp)
        
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝐽(x_)*𝑤
        # println(𝑓(x_)*div𝒗(x_))
    end
    # println(temp)
end

for a in elements["Γ"]
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        x_ = Vec{2}((ξ.x,ξ.y))
        nₐ_ = Vec{2}((ξ.n₁,ξ.n₂))
        global err -= 𝑓(x_)*𝑣ᵃ(x_)⋅nₐ_*𝑤
        global temp2 -= 𝑓(x_)*𝑣ᵃ(x_)⋅nₐ_*𝑤
    end
end

println(err)