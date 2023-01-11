using QuantumOptics
using CairoMakie
using Colors, ColorSchemes



φ = 2*π/3
α₁ = 2.3im
α₂ = α₁*exp(1im*2*π/3)
α₃ = α₁*exp(-1im*2*π/3)

Nfock = 50
fock_b = FockBasis(Nfock)
cat_1 = normalize(coherentstate(fock_b, α₁) + coherentstate(fock_b, α₂))
cat_2 = normalize(coherentstate(fock_b, α₂) + coherentstate(fock_b, α₃))
cat_3 = normalize(coherentstate(fock_b, α₃) + coherentstate(fock_b, α₁))

jgreen = HSLA(111, 0.60, 0.373, 1.0)
jpurple = HSLA(281, 0.369, 0.522, 1.0)
jred = HSLA(4, 0.598, 0.498, 1.0)
trans = HSLA(0, 0, 0, 0.0)

myjuliacolor1 = cgrad([jgreen, trans, jpurple], 3, categorical=true)
myjuliacolor2 = cgrad([jpurple, trans, jred], 3, categorical=true)
myjuliacolor3 = cgrad([jred, trans, jgreen], 3, categorical=true)

f2 = Figure(resolution = (800, 800))
ax = Axis(f2[1,1],
aspect = 1
)
W = QuantumOptics.wigner(cat_1, -4:0.01:4, -4:0.01:4)
w_lim = maximum(x->isnan(x) ? -Inf : x, W)
xpoints = reduce(hcat, [collect(-4:0.01:4) for i in 1:length(-4:0.01:4)])
contourf!(ax, vcat(xpoints'...), vcat(xpoints...), vcat(W...); interpolation=false, transparency=true, levels=3, colormap = myjuliacolor1, colorrange = (-w_lim, w_lim))
W = QuantumOptics.wigner(cat_2, -4:0.01:4, -4:0.01:4)
w_lim = maximum(x->isnan(x) ? -Inf : x, W)
contourf!(ax, vcat(xpoints'...), vcat(xpoints...), vcat(W...); interpolation=false, transparency=true, levels=3, colormap = myjuliacolor2, colorrange = (-w_lim, w_lim))
W = QuantumOptics.wigner(cat_3, -4:0.01:4, -4:0.01:4)
w_lim = maximum(x->isnan(x) ? -Inf : x, W)
contourf!(ax, vcat(xpoints'...), vcat(xpoints...), vcat(W...); interpolation=false, transparency=true, levels=3, colormap = myjuliacolor3, colorrange = (-w_lim, w_lim))
#W = QuantumOptics.wigner(cat_2, -6:0.01:6, -6:0.01:6)
#w_lim = maximum(x->isnan(x) ? -Inf : x, W)
#heatmap!(ax, -4:0.1:4, -4:0.1:4, W; interpolation=true, colormap = myjuliacolor2, colorrange = (-w_lim, w_lim))
#W = QuantumOptics.wigner(cat_3, -6:0.01:6, -6:0.01:6)
#w_lim = maximum(x->isnan(x) ? -Inf : x, W)
#heatmap!(ax, -4:0.1:4, -4:0.1:4, W; interpolation=true, colormap = myjuliacolor3, colorrange = (-w_lim, w_lim))
f2