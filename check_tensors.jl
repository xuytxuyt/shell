using Tensors

κ₁₁ = rand()
κ₂₂ = rand()
κ₁₂ = rand()
𝜿ten = SymmetricTensor{2,2}((κ₁₁,κ₁₂,κ₂₂))
𝜿mat = [κ₁₁,κ₂₂,2*κ₁₂]
C¹¹¹¹ = rand()
C²²²² = rand()
C¹¹²² = rand()
C¹¹¹² = rand()
C²²¹² = rand()
C¹²¹² = rand()
Cmat = [C¹¹¹¹ C¹¹²² C¹¹¹²;C¹¹²² C²²²² C²²¹²;C¹¹¹² C²²¹² C¹²¹²]
Cten = SymmetricTensor{4,2}((C¹¹¹¹,C¹¹¹²,C¹¹²²,C¹¹¹²,C¹²¹²,C²²¹²,C¹¹²²,C²²¹²,C²²²²))

𝝈mat = Cmat*𝜿mat
𝝈ten = Cten ⊡ 𝜿ten
check = 𝝈mat[1] - 𝝈ten[1,1];check ≈ 0.0 ? println("pass") : println("𝝈 wrong: $check") 
check = 𝝈mat[2] - 𝝈ten[2,2];check ≈ 0.0 ? println("pass") : println("𝝈 wrong: $check") 
check = 𝝈mat[3] - 𝝈ten[1,2];check ≈ 0.0 ? println("pass") : println("𝝈 wrong: $check") 