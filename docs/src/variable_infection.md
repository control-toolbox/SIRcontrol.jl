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

plot(beta_ , 0, 300, label = L"$\beta$", color = :black, lw = 5,
    grid = true,
    gridlinewidth = 1.5,
        size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)
plot!(t -> beta(t, 150, 0.8, 0.4, 3), label = L"$\beta_1$", 
color = :orange, linestyle = :dash, lw = 5)
plot!(t -> beta(t, 150, 0.8, 0.4, 300), label = L"$\beta_{300}$", 
color = :red, linestyle = :dash, lw = 5)
xlims!(139, 161)
```

```@example main
tf   = 700
τ     = 50
γ    = 0.15
Q     = 35
umax = 0.8

β1   = [0.4, 0.7, 0.4, 0.2]
β2   = [0.2, 0.4, 0.8, 0.8]

S0   = 0.999
I0   = 0.001
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
	for N in [50, 100 ,600, 4000]
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
k=300
sol1 = SIRocp2beta(tf, S0, I0, Q, τ, β1[1], β2[1], γ, umax, k)
sol2 = SIRocp2beta(tf, S0, I0, Q, τ, β1[2], β2[2], γ, umax, k)
sol3 = SIRocp2beta(tf, S0, I0, Q, τ, β1[3], β2[3], γ, umax, k)
sol4 = SIRocp2beta(tf, S0, I0, Q, τ, β1[4], β2[4], γ, umax, k)
```

## case 1: $\beta_1 > \beta_2$ (two interventions)


```@example main

function Sh1(t)
    if t<=τ
        return γ/β1[1]
    else 
        return γ/β2[1]
    end
end

gr()


plot(
    t -> state(sol1)(t)[1], 0, tf,
    label = false, color = :darkblue, lw = 5,
    grid = true,
    gridlinewidth = 1.5,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol1)(t)[2], 0, tf,
    label = false, color = :darkgreen, lw = 5
)


plot!(
    t -> control(sol1)(t)[1], 0, tf,
    label = false, color = :darkred, lw = 5
)

plot!(t -> Sh1(t), 0, tf,     
    label = false, color = :black, ls = :dash, lw = 3
)

plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
xlims!(0, 126)
```

## case 2: $\beta_1 > \beta_2$ (one interventions)

```@example main

function Sh2(t)
    if t<=τ
        return γ/β1[2]
    else 
        return γ/β2[2]
    end
end
gr()

plot(
    t -> state(sol2)(t)[1], 0, tf,
    label = false, color = :darkblue, lw = 5,
    grid = true,
    gridlinewidth = 1.5,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol2)(t)[2], 0, tf,
    label = false, color = :darkgreen, lw = 5
)


plot!(
    t -> control(sol2)(t)[1], 0, tf,
    label = false, color = :darkred, lw = 5
)

plot!(t -> Sh2(t), 0, tf,     
    label = false, color = :black, ls = :dash, lw = 3 
)

#vline!(
#    [τ], lw = 3, color = :black, linestyle = :dashdot, label = false
#)

plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
xlims!(0,106)
```


## case 3: $\beta_1 < \beta_2$ (two interventions)


```@example main

function Sh3(t)
    if t<=τ
        return γ/β1[3]
    else 
        return γ/β2[3]
    end
end
gr()
plot(
    t -> state(sol3)(t)[1], 0, tf,
    label = false, color = :darkblue, lw = 5,
    grid = true,
    gridlinewidth = 1.5,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol3)(t)[2], 0, tf,
    label = false, color = :darkgreen, lw = 5
)

plot!(
    t -> control(sol3)(t)[1], 0, tf,
    label = false, color = :darkred, lw = 5
)
plot!(t -> Sh3(t), 0, tf,     
    label = false, color = :black, ls = :dash, lw = 3
)

# Dummy legend entries
plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h$",   lw = 1.5, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
xlims!(0, 126)
```

## case 4: $\beta_1 < \beta_2$ (one interventions)

```@example main

function Sh4(t)
    if t<=τ
        return γ/β1[4]
    else 
        return γ/β2[4]
    end
end

gr()

plot(
    t -> state(sol4)(t)[1], 0, tf,
    label = false, color = :darkblue, lw = 5,
    grid = true,
    gridlinewidth = 1.5,
    size = (1600, 1000),
    legendfontsize = 24,
    tickfontsize = 20,
    guidefontsize = 26,
    framestyle = :box,
    foreground_color_subplot = :black
)

plot!(
    t -> state(sol4)(t)[2], 0, tf,
    label = false, color = :darkgreen, lw = 5
)


plot!(
    t -> control(sol4)(t)[1], 0, tf,
    label = false, color = :darkred, lw = 5
)

plot!(t -> Sh4(t), 0, tf,     
    label = false, color = :black, ls = :dash, lw = 3 
)


plot!([NaN], [NaN], label = L"$S$",       lw = 2, color = :darkblue)
plot!([NaN], [NaN], label = L"$I$",       lw = 2, color = :darkgreen)
plot!([NaN], [NaN], label = L"$S_h$",   lw = 2, ls = :dash, color = :black)
plot!([NaN], [NaN], label = L"$u$",       lw = 2, color = :darkred)
xlims!(0, 126)
```