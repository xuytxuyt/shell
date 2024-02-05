
opsGauss = quote
    op = Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h)
end

opsMix = quote
    ops = [
        Operator{:∫δεCεdΩ}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫δNεdΩ}(),
        Operator{:∫𝒏𝑵𝒗dΓ}(),
        Operator{:∫∇𝑵𝒗dΩ}(),
        Operator{:∫δκCκdΩ}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫δMκdΩ}(),
        Operator{:∫∇𝑴𝒏𝒂₃𝒗dΓ}(),
        Operator{:∫𝑴ₙₙ𝜽ₙdΓ}(),
        Operator{:ΔMₙₛ𝒂₃𝒗}(),
        Operator{:∫∇𝑴∇𝒂₃𝒗dΩ}(),
    ]
end

opsPenalty = quote
    opsα = [
        Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ*E),
        Operator{:∫δθθdΓ}(:α=>αᵣ*E),
    ]
end

opsNitsche = quote
    opsv = [
        Operator{:∫𝒏𝑵𝒈dΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫∇𝑴𝒏𝒂₃𝒈dΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
        Operator{:∫MₙₙθₙdΓ_Nitsche}(:E=>E,:ν=>ν,:h=>h),
    ]
end

opsHR = quote
    opsh = [
        Operator{:∫𝒏𝑵𝒈dΓ_HR}(),
        Operator{:∫∇𝑴𝒏𝒂₃𝒈dΓ_HR}(),
        Operator{:∫MₙₙθₙdΓ_HR}(),
        Operator{:ΔMₙₛ𝒂₃𝒈_HR}(),
    ]
end