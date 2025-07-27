```@example main
using OptimalControl
using NLPModelsIpopt
using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Roots
```

We consider a piecewise constant transmission rate that is given by

```math
\beta(t) = \begin{cases}
\beta_1, \quad \text{if } t < T_c \\[5pt]
\beta_2, \quad \text{if } t > T_c,
\end{cases}
```

Furthermore, we consider its regularization using a sigmoid function:

```@example main
function beta(t, τ, β1, β2, k)
    return β1 + (β2 - β1) / (1 + exp(-k * (t - τ))) 
end
```

```@example main
function beta_(t)
    if t<=150
        return 0.8
    else 
        return 0.4
    end
end

k = 300 # regularization parameter

plot(beta_ , 0, 300, label = L"$\beta$", color = :black, lw = 5,
    grid = false,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(t -> beta(t, 150, 0.8, 0.4, k), label = L"$\beta_k$", 
color = :red, linestyle = :dash, lw = 5, grid = false)
```

```@example main
tf_   = 300
τ     = 50
γ_    = 0.15
β11   = 0.4
β12   = 0.2
S0_   = 0.999
I0_   = 0.001
Q_    = 35
β21   = 0.7
β22   = 0.4
umax_ = 0.8
```

```@example main
function SIRocp2beta(T, S0, I0, Q, τ, β1, β2, γ, umax,k)
	
	ocp = @def begin
	        t ∈ [0, T], time
	        x = (S, I, C) ∈ R^3, state
	        u ∈ R, control
	
	        S(0) == S0
	        I(0) == I0
	        C(0) == Q
	
	        C(T) ≥ 0.0
	
	        0 ≤ u(t) ≤ umax
	
	        0.0 ≤ I(t) ≤ 1.0
	        0.0 ≤ S(t) ≤ 1.0
	
	        ẋ(t) == [-(1 - u(t)) * beta(t, τ, β1, β2, k) * S(t) * I(t),
	                  (1 - u(t)) * beta(t, τ, β1, β2, k) * S(t) * I(t) - γ * I(t),
	                -u(t)]
	
	        -S(T) → min
	end
	
	sol = nothing
	for N in [50, 100, 200,300, 400, 500,600,700,800,900,1000]
	    sol = solve(ocp, :direct, :adnlp, :ipopt; 
	                    disc_method = :gauss_legendre_3, 
	                    grid_size=N, 
	                    init=sol,
                        tol=1e-8, 
	                    display=true)
	end
	return sol
end
```

```@example main
sol1_ = SIRocp2beta(tf_, S0_, I0_, Q_, τ, β11, β12, γ_, umax_, k)
sol2_ = SIRocp2beta(tf_, S0_, I0_, Q_, τ, β21, β22, γ_, umax_, k)
sol3_ = SIRocp2beta(tf_, S0_, I0_, Q_, 150, 0.4, 0.8, γ_, umax_, k)
sol4_ = SIRocp2beta(tf_, S0_, I0_, Q_, 100, 0.2, 0.8, γ_, umax_, k)
```

## case 1: $\beta_1 > \beta_2$ (two interventions)


```@example main
gr()

plot(
    t -> state(sol1_)(t)[1], 0, tf_,
    label = false, color = :darkblue, lw = 5,
    grid = false,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol1_)(t)[2], 0, tf_,
    label = false, color = :darkgreen, lw = 5
)

plot!(
    [0, tf_], [γ_/β11, γ_/β11],
    label = false, lw = 3, ls = :dash, color = :black
)

plot!(
    [0, tf_], [γ_/β12, γ_/β12],
    label = false, lw = 3, ls = :dash, color = :purple
)

plot!(
    t -> control(sol1_)(t)[1], 0, tf_,
    label = false, color = :darkred, lw = 5
)

vline!(
    [τ], lw = 3, color = :black, linestyle = :dashdot, label = false
)

plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h^1$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$S_h^2$",   lw = 1.5, ls = :dash, color = :purple)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
plot!([NaN], [NaN], label = L"$\tau$",    lw = 1.5, ls = :dashdot, color = :black)
```

## case 2: $\beta_1 > \beta_2$ (one interventions)

```@example main
gr()

plot(
    t -> state(sol2_)(t)[1], 0, tf_,
    label = false, color = :darkblue, lw = 5,
    grid = false,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol2_)(t)[2], 0, tf_,
    label = false, color = :darkgreen, lw = 5
)

plot!(
    [0, tf_], [γ_/β21, γ_/β21],
    label = false, lw = 3, ls = :dash, color = :black
)

plot!(
    [0, tf_], [γ_/β22, γ_/β22],
    label = false, lw = 3, ls = :dash, color = :purple
)

plot!(
    t -> control(sol2_)(t)[1], 0, tf_,
    label = false, color = :darkred, lw = 5
)

vline!(
    [τ], lw = 3, color = :black, linestyle = :dashdot, label = false
)

plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h^1$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$S_h^2$",   lw = 1.5, ls = :dash, color = :purple)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
plot!([NaN], [NaN], label = L"$\tau$",    lw = 1.5, ls = :dashdot, color = :black)
```


## case 3: $\beta_1 < \beta_2$ (two interventions)


```@example main
gr()
plot(
    t -> state(sol3_)(t)[1], 0, tf_,
    label = false, color = :darkblue, lw = 5,
    grid = false,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol3_)(t)[2], 0, tf_,
    label = false, color = :darkgreen, lw = 5
)

plot!(
    [0, tf_], [γ_/0.4, γ_/0.4],
    label = false, lw = 3, ls = :dash, color = :black
)

plot!(
    [0, tf_], [γ_/0.8, γ_/0.8],
    label = false, lw = 3, ls = :dash, color = :purple
)

plot!(
    t -> control(sol3_)(t)[1], 0, tf_,
    label = false, color = :darkred, lw = 5
)

vline!(
    [150], lw = 3, color = :black, linestyle = :dashdot, label = false
)

# Dummy legend entries
plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h^1$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$S_h^2$",   lw = 1.5, ls = :dash, color = :purple)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
plot!([NaN], [NaN], label = L"$\tau$",    lw = 1.5, ls = :dashdot, color = :black)
```

## case 4: $\beta_1 < \beta_2$ (one interventions)

```@example main
gr()

plot(
    t -> state(sol4_)(t)[1], 0, tf_,
    label = false, color = :darkblue, lw = 5,
    grid = false,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol4_)(t)[2], 0, tf_,
    label = false, color = :darkgreen, lw = 5
)

plot!(
    [0, tf_], [γ_/0.2, γ_/0.2],
    label = false, lw = 3, ls = :dash, color = :black
)

plot!(
    [0, tf_], [γ_/0.8, γ_/0.8],
    label = false, lw = 3, ls = :dash, color = :purple
)

plot!(
    t -> control(sol4_)(t)[1], 0, tf_,
    label = false, color = :darkred, lw = 5
)

vline!(
    [100], lw = 3, color = :black, linestyle = :dashdot, label = false
)

plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h^1$",   lw = 2, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$S_h^2$",   lw = 2, ls = :dash, color = :purple)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
plot!([NaN], [NaN], label = L"$\tau$",    lw = 2, ls = :dashdot, color = :black)
```