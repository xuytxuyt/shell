using ApproxOperator, Tensors
import BenchmarkExample: BenchmarkExample

include("import_vcc.jl")


elements, nodes = import_vcc("msh/scordelislo_16.msh")
cs = BenchmarkExample.cylindricalCoordinate(25.0)
# elements, nodes = import_vcc("msh/sphericalshell_1.msh")
# cs = BenchmarkExample.sphericalCoordinate(10.0)


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
Γ¹₁₁₁ = cs.∂₁Γ¹₁₁
Γ¹₁₁₂ = cs.∂₂Γ¹₁₁
Γ¹₁₂₁ = cs.∂₁Γ¹₁₂
Γ¹₁₂₂ = cs.∂₂Γ¹₁₂
Γ²₂₁₁ = cs.∂₁Γ²₁₂
Γ²₂₁₂ = cs.∂₂Γ²₁₂
Γ²₂₂₁ = cs.∂₁Γ²₂₂
Γ²₂₂₂ = cs.∂₂Γ²₂₂
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

# 𝑓(x) = 1.0+2.0*x[1]+3.0*x[2] + 4.0*x[1]^2 + 5.0*x[1]*x[2] + 6.0*x[2]^2
# 𝑓(x) = 1.0
𝑓(x) = (x[1]+x[2])^2
function ∂ₐ𝑓(x)
    x_ = Vec{2}((x[1],x[2]))
    gradient(𝑓,x_)
end
function ∂ₐᵦ𝑓(x)
    x_ = Vec{2}((x[1],x[2]))
    gradient(∂ₐ𝑓,x_)
end
𝑔(x) = (x[1]+x[2])^3
∂ₐ𝑔(x) = gradient(𝑔,x)

function ∇𝑓(x)
    grad = ∂ₐ𝑓(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return grad[1]*𝒂¹(x_) + grad[2]*𝒂²(x_)
end

𝑣ᵃ(x) = Vec{2}((4.0+5.0*x[1]+6.0*x[2],7.0+8.0*x[1]+9.0*x[2]))
# 𝑣ᵃ(x) = Vec{2}((1.0,0.0))
# Aᵃᵇ(x) = SymmetricTensor{2,2}((1.0+2.0*x[1]+3.0*x[2],2.0+3.0*x[1]+4.0*x[2],3.0+4.0*x[1]+5.0*x[2]))
# Aᵃᵇ(x) = Tensor{2,2}((1.0*x[1]^1,2.0*x[2]^1,3.0,4.0))
Aᵃᵇ(x) = Tensor{2,2}((1.0*x[1]^2,2.0*x[2]^2,2.0*x[2]^2,3.0*x[1]*x[2]))
function ∂Aᵃᵇ(x)
    x_ = Vec{2}((x[1],x[2]))
    gradient(Aᵃᵇ,x_)
end
function divAᵃᵇ(x)
    x_ = Vec{2}((x[1],x[2]))
    gradA = ∂Aᵃᵇ(x_)
    return Vec{2}(gradA[:,1,1] + gradA[:,2,2])
end
function ∂ₐᵦAᵃᵇ(x)
    x_ = Vec{2}((x[1],x[2]))
    divergence(divAᵃᵇ,x_)
end

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

function Γᵞᵧᵦₐ(x)
    x_ = Vec{3}((x[1],x[2],0.0))
    return Tensor{2,2}((Γ¹₁₁₁(x_)+Γ²₂₁₁(x_),Γ¹₁₂₁(x_)+Γ²₂₂₁(x_),Γ¹₁₁₂(x_)+Γ²₂₁₂(x_),Γ¹₁₂₂(x_)+Γ²₂₂₂(x_)))
end

x = Vec{2}((rand(),rand()))
x_ = Vec{3}((rand(),rand(),rand()))
(check1 = ∇𝑓(x)⋅𝒗(x) - ∂ₐ𝑓(x)⋅𝑣ᵃ(x);check1 ≈ 0.0) ? println("Check 1 right!: $check1") : println("Check 1 wrong: $check1")
(check2 = 𝒗(x)⋅𝒏(x) - 𝑣ᵃ(x)⋅nₐ(x);check2 ≈ 0.0) ? println("Check 2 right!: $check2") : println("Check 2 wrong: $check2")
(check3 = div𝒗(x) - 𝑣ᵃₐ(x) - Γᵞᵧₐ(x)⋅𝑣ᵃ(x);check3 ≈ 0.0) ? println("Check 3 right!: $check3") : println("Check 3 wrong: $check3")
(check4 = 𝒂¹(x_) - 𝒂₂(x_)×𝒂₃(x_)/𝐽(x_);check4 ≈ Vec{3}((0.0,0.0,0.0))) ? println("Check 4 right!: $check4") : println("Check 4 wrong: $check4")
(check5 = 𝒂²(x_) - 𝒂₃(x_)×𝒂₁(x_)/𝐽(x_);check5 ≈ Vec{3}((0.0,0.0,0.0))) ? println("Check 5 right!: $check5") : println("Check 5 wrong: $check5")
(check6 = 𝒂³(x_) - 𝒂₁(x_)×𝒂₂(x_)/𝐽(x_);check6 ≈ Vec{3}((0.0,0.0,0.0))) ? println("Check 6 right!: $check6") : println("Check 6 wrong: $check6")

err = 0.0
temp1 = 0.0
temp2 = 0.0
for a in elements["Ω"]
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        x_ = Vec{2}((ξ.x,ξ.y))
        global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝑤
        # global temp1 += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝑤
        # println(Γᵞᵧₐ(x_))
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*𝑣ᵃₐ(x_))*𝑤
        # global temp += 𝑓(x_)*Γᵞᵧₐ(x_)⋅𝑣ᵃ(x_)*𝑤
        # println(temp)
        
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝐽(x_)*𝑤
        # println(𝑓(x_)*div𝒗(x_))

        # global err += (∂ₐ𝑓(x_)[1]*𝑔(x_) + 𝑓(x_)*∂ₐ𝑔(x_)[1])*𝑤
        # global err += ∂ₐᵦ𝑓(x_) ⊡ Aᵃᵇ(x_)*𝑤 - 𝑓(x_)*(∂ₐᵦAᵃᵇ(x_) + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_) + 2*Γᵞᵧₐ(x_) ⋅ divAᵃᵇ(x_) + (Γᵞᵧₐ(x_)⊗Γᵞᵧₐ(x_)) ⊡ Aᵃᵇ(x_))*𝑤
        # global temp1 += ∂ₐᵦ𝑓(x_) ⊡ Aᵃᵇ(x_)*𝑤 - 𝑓(x_)*(∂ₐᵦAᵃᵇ(x_) + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_) + 2*Γᵞᵧₐ(x_) ⋅ divAᵃᵇ(x_) + (Γᵞᵧₐ(x_)⊗Γᵞᵧₐ(x_)) ⊡ Aᵃᵇ(x_))*𝑤
        # global temp1 -= 𝑓(x_)*(∂ₐᵦAᵃᵇ(x_) + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_) + 2*Γᵞᵧₐ(x_) ⋅ divAᵃᵇ(x_) + (Γᵞᵧₐ(x_)⊗Γᵞᵧₐ(x_)) ⊡ Aᵃᵇ(x_))*𝑤
    end
    # println(temp)
end

for a in elements["Γ"]
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        x_ = Vec{2}((ξ.x,ξ.y))
        nₐ_ = Vec{2}((ξ.n₁,ξ.n₂))
        nᵃ_ = Vec{2}((ξ.n¹,ξ.n²))
        sₐ_ = Vec{2}((ξ.s₁,ξ.s₂))
        sᵃ_ = Vec{2}((ξ.s¹,ξ.s²))
        # println((ξ.∂₁n₁,ξ.∂₁n₂,ξ.∂₂n₁,ξ.∂₂n₂))
        ∂ᵦnₐ_ = Tensor{2,2}((ξ.∂₁n₁,ξ.∂₁n₂,ξ.∂₂n₁,ξ.∂₂n₂))
        ∂ᵦsₐ_ = Tensor{2,2}((ξ.∂₁s₁,ξ.∂₁s₂,ξ.∂₂s₁,ξ.∂₂s₂))
        # println(nₐ_)
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n¹ = ξ.n¹
        n² = ξ.n²
        s₁ = ξ.s₁
        s₂ = ξ.s₂
        s¹ = ξ.s¹
        s² = ξ.s²
        Δ = ξ.Δ
        n = n₁*n¹+n₂*n²
        s = s₁*s¹+s₂*s²
        # println("n: $n, s: $s")
        global err -= 𝑓(x_)*𝑣ᵃ(x_)⋅nₐ_*𝑤
        # global temp2 -= 𝑓(x_)*𝑣ᵃ(x_)⋅nₐ_*𝑤
        # global temp2 -= ξ.n₁*𝑤
        # global temp2 -= ξ.n₂*𝑤
        # global err -= 𝑓(x_)*𝑔(x_)*nₐ_[1]*𝑤
        # global err -= ∂ₐ𝑓(x_)⋅nᵃ_*nₐ_⋅Aᵃᵇ(x_)⋅nₐ_*𝑤 - 𝑓(x_)*(sₐ_⋅∂Aᵃᵇ(x_) ⊡ (nₐ_⊗sᵃ_) + nₐ_⋅Aᵃᵇ(x_)⋅∂ᵦsₐ_⋅sᵃ_ + sₐ_⋅Aᵃᵇ(x_)⋅∂ᵦnₐ_⋅sᵃ_)*𝑤 + Δ*𝑓(x_)*sₐ_⋅Aᵃᵇ(x_)⋅nₐ_ - Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅nₐ_*𝑤
        # global temp2 -= - 𝑓(x_)*(nₐ_⋅Aᵃᵇ(x_)⋅∂ᵦsₐ_⋅sᵃ_ + sₐ_⋅Aᵃᵇ(x_)⋅∂ᵦnₐ_⋅sᵃ_)*𝑤 + Δ*𝑓(x_)*sₐ_⋅Aᵃᵇ(x_)⋅nₐ_ - Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅nₐ_*𝑤
        # global temp2 -= - 𝑓(x_)*(nₐ_⋅Aᵃᵇ(x_)⋅∂ᵦsₐ_⋅sᵃ_ + sₐ_⋅Aᵃᵇ(x_)⋅∂ᵦnₐ_⋅sᵃ_)*𝑤 + Δ*𝑓(x_)*sₐ_⋅Aᵃᵇ(x_)⋅nₐ_
    end
end

println(err)
# println(temp2)

temp1 = 0.0
temp2 = 0.0
for a in elements["Γ"]
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    𝐿 = ((x₁-x₂)^2 + (y₁-y₂)^2)^0.5
    t¹ = (x₂-x₁)/𝐿
    t² = (y₂-y₁)/𝐿
    𝐿 ≠ a.𝐿 ? println("𝐿 wrong!") : nothing
    t₁(x) = cs.a₁₁(x)*t¹ + cs.a₁₂(x)*t²
    t₂(x) = cs.a₁₂(x)*t¹ + cs.a₂₂(x)*t²
    t(x) = (t₁(x)*t¹ + t₂(x)*t²)^0.5
    s¹(x) = t¹/t(x)
    s²(x) = t²/t(x)
    s₁(x) = t₁(x)/t(x)
    s₂(x) = t₂(x)/t(x)
    deta(x) = (cs.a₁₁(x)*cs.a₂₂(x) - cs.a₁₂(x)^2)^0.5
    n₁(x) = s²(x)*deta(x)
    n₂(x) =-s¹(x)*deta(x)
    n¹(x) = cs.a¹¹(x)*n₁(x) + cs.a¹²(x)*n₂(x)
    n²(x) = cs.a¹²(x)*n₁(x) + cs.a²²(x)*n₂(x)
    ∂₁n₁(x) = gradient(n₁,x)[1]
    ∂₂n₁(x) = gradient(n₁,x)[2]
    ∂₁n₂(x) = gradient(n₂,x)[1]
    ∂₂n₂(x) = gradient(n₂,x)[2]
    ∂₁s₁(x) = gradient(s₁,x)[1]
    ∂₂s₁(x) = gradient(s₁,x)[2]
    ∂₁s₂(x) = gradient(s₂,x)[1]
    ∂₂s₂(x) = gradient(s₂,x)[2]
    ∂₁s¹(x) = gradient(s¹,x)[1]
    ∂₂s¹(x) = gradient(s¹,x)[2]
    ∂₁s²(x) = gradient(s²,x)[1]
    ∂₂s²(x) = gradient(s²,x)[2]
    𝒏ₐ(x) = Vec{2}((n₁(x),n₂(x)))
    𝒏ᵃ(x) = Vec{2}((n¹(x),n²(x)))
    𝒔ₐ(x) = Vec{2}((s₁(x),s₂(x)))
    𝒔ᵃ(x) = Vec{2}((s¹(x),s²(x)))
    Mₛₙ(x) = 𝒔ₐ(x)⋅Aᵃᵇ(x)⋅𝒏ₐ(x)
    ∂Mₛₙ(x) = gradient(Mₛₙ,x)[1:2]
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        Δ = ξ.Δ
        x_ = Vec{3}((ξ.x,ξ.y,ξ.z))
        # check = s₁(x_) - ξ.s₁; check ≉ 0.0 ? println("s₁ wrong! $check") : nothing
        # check = s₂(x_) - ξ.s₂; check ≉ 0.0 ? println("s₂ wrong! $check") : nothing
        # check = n₁(x_) - ξ.n₁; check ≉ 0.0 ? println("n₁ wrong! $check") : nothing
        # check = n₂(x_) - ξ.n₂; check ≉ 0.0 ? println("n₂ wrong! $check") : nothing
        # check = t(x_) - ξ.t; check ≉ 0.0 ? println("t wrong! $check") : nothing
        # check = ∂₁n₁(x_) - ξ.∂₁n₁; check ≉ 0.0 ? println("∂₁n₁ wrong! $check") : nothing
        # check = ∂₂n₁(x_) - ξ.∂₂n₁; check ≉ 0.0 ? println("∂₂n₁ wrong! $check") : nothing
        # check = ∂₁n₂(x_) - ξ.∂₁n₂; check ≉ 0.0 ? println("∂₁n₂ wrong! $check") : nothing
        # check = ∂₂n₂(x_) - ξ.∂₂n₂; check ≉ 0.0 ? println("∂₂n₂ wrong! $check") : nothing
        # check = ∂₁s₁(x_) - ξ.∂₁s₁; check ≉ 0.0 ? println("∂₁s₁ wrong! $check") : nothing
        # check = ∂₂s₁(x_) - ξ.∂₂s₁; check ≉ 0.0 ? println("∂₂s₁ wrong! $check") : nothing
        # check = ∂₁s₂(x_) - ξ.∂₁s₂; check ≉ 0.0 ? println("∂₁s₂ wrong! $check") : nothing
        # check = ∂₂s₂(x_) - ξ.∂₂s₂; check ≉ 0.0 ? println("∂₂s₂ wrong! $check") : nothing
        # println(gradient(n₁,x_))
        # println(ξ.∂₂n₁)
        ∂𝒔ₐ = Vec{2}((ξ.∂₁s₁*ξ.s¹ + ξ.∂₂s₁*ξ.s²,ξ.∂₁s₂*ξ.s¹ + ξ.∂₂s₂*ξ.s²))
        ∂𝒏ₐ = Vec{2}((ξ.∂₁n₁*ξ.s¹ + ξ.∂₂n₁*ξ.s²,ξ.∂₁n₂*ξ.s¹ + ξ.∂₂n₂*ξ.s²))
        # global temp1 += - 𝑓(x_)*∂Mₛₙ(x_)⋅𝒔ᵃ(x_)*𝑤 + 𝑓(x_)*Mₛₙ(x_)*Δ
        # global temp1 += - 𝑓(x_)*∂Mₛₙ(x_)⋅𝒔ᵃ(x_)*𝑤
        # global temp1 += - 𝑓(x_)*∂Mₛₙ(x_)⋅𝒔ᵃ(x_)*𝑤 + 𝑓(x_)*Mₛₙ(x_)*Δ - 𝑓(x_)*Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤
        # global temp1 -= 𝑓(x_)*Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤 + 𝑓(x_)*divAᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤
        # global temp1 += - 𝑓(x_)*∂Mₛₙ(x_)⋅𝒔ᵃ(x_)*𝑤 + 𝑓(x_)*Mₛₙ(x_)*Δ + ∂ₐ𝑓(x_)⋅𝒏ᵃ(x_)*(𝒏ₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_))*𝑤 - 𝑓(x_)*Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤 - 𝑓(x_)*divAᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤
        # println(∂Aᵃᵇ(x_)[:,:,1])
        # println(𝒔ₐ(x_))
        global temp1 += - 𝑓(x_)*(𝒔ₐ(x_)⋅(Tensor{2,2}(∂Aᵃᵇ(x_)[:,:,1].*ξ.s¹+∂Aᵃᵇ(x_)[:,:,2].*ξ.s²))⋅𝒏ₐ(x_) + ∂𝒔ₐ⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_) + 𝒔ₐ(x_)⋅Aᵃᵇ(x_)⋅∂𝒏ₐ)*𝑤 + 𝑓(x_)*Mₛₙ(x_)*Δ + ∂ₐ𝑓(x_)⋅𝒏ᵃ(x_)*(𝒏ₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_))*𝑤 - 𝑓(x_)*Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤 - 𝑓(x_)*divAᵃᵇ(x_)⋅𝒏ₐ(x_)*𝑤
        # println(∂Mₛₙ(x_))
        # println(Aᵃᵇ(x_))
        # println(Δ)
        # println(𝑓(x_))
        # println(Mₛₙ(x_))
        # println(𝑤)
        # global temp2 += - 𝑓(x_)*(∂𝒔ₐ⋅Aᵃᵇ(x_)⋅𝒏ₐ(x_) + 𝒔ₐ(x_)⋅Aᵃᵇ(x_)⋅∂𝒏ₐ)*𝑤 + 𝑓(x_)*Mₛₙ(x_)*Δ
        # global temp2 += 𝑓(x_)*Mₛₙ(x_)*Δ
    end
end

for a in elements["Ω"]
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        x_ = Vec{3}((ξ.x,ξ.y,ξ.z))
        # global temp1 += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝑤
        # println(Γᵞᵧₐ(x_))
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*𝑣ᵃₐ(x_))*𝑤
        # global temp += 𝑓(x_)*Γᵞᵧₐ(x_)⋅𝑣ᵃ(x_)*𝑤
        # println(temp)
        
        # global err += (∇𝑓(x_)⋅𝒗(x_) + 𝑓(x_)*div𝒗(x_))*𝐽(x_)*𝑤
        # println(𝑓(x_)*div𝒗(x_))

        # global err += (∂ₐ𝑓(x_)[1]*𝑔(x_) + 𝑓(x_)*∂ₐ𝑔(x_)[1])*𝑤
        # global err += ∂ₐᵦ𝑓(x_) ⊡ Aᵃᵇ(x_)*𝑤 - 𝑓(x_)*(∂ₐᵦAᵃᵇ(x_) + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_) + 2*Γᵞᵧₐ(x_) ⋅ divAᵃᵇ(x_) + (Γᵞᵧₐ(x_)⊗Γᵞᵧₐ(x_)) ⊡ Aᵃᵇ(x_))*𝑤
        global temp1 -= ∂ₐᵦ𝑓(x_) ⊡ Aᵃᵇ(x_)*𝑤
        # println(Γ¹₁₁₁(x_))
        # println(Γ¹₁₁₂(x_))
        # println(Γ¹₁₂₁(x_))
        # println(Γ¹₁₂₂(x_))
        # println(Γ²₂₁₁(x_))
        # println(Γ²₂₁₂(x_))
        # println(Γ²₂₂₁(x_))
        # println(Γ²₂₂₂(x_))
        # println(temp1)
        global temp1 += 𝑓(x_)*(Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅Γᵞᵧₐ(x_) + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_) + 2*Γᵞᵧₐ(x_)⋅divAᵃᵇ(x_) + ∂ₐᵦAᵃᵇ(x_))*𝑤
        # global temp1 += Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅Γᵞᵧₐ(x_)*𝑤 + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_)*𝑤 
        # global temp2 -= Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅Γᵞᵧₐ(x_)*𝑤 + Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_)*𝑤 
        # global temp2 -= Γᵞᵧₐ(x_)⋅Aᵃᵇ(x_)⋅Γᵞᵧₐ(x_)*𝑤
        # global temp2 -=  Γᵞᵧᵦₐ(x_) ⊡ Aᵃᵇ(x_)*𝑤 
    end
    # println(temp)
end
println(temp1)
# println(temp2)