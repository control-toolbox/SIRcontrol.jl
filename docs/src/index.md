# On the problem of minimizing the epidemic final size for SIR model via social distancing

## 📘 Introduction

### Context and Motivation

We revisit the problem of minimizing the epidemic final size in the **SIR model** through **social distancing interventions**.
Traditionally, this problem assumes a fixed **interval structure** for the timing of interventions. In contrast, we investigate a more flexible approach: **L¹-constrained controls** that limit total intervention effort, rather than when and how it's applied.

**Key insight:**
> Even with this generalization, the optimal control still occurs over a **single time interval** when the transmission rate is constant, and over at most **two (separate) time intervals** if the transmission rate changes once.

---

## The SIR Model with Control

### System Dynamics

We consider the **SIR model** with a control variable _u_(_t_) (representing social distancing):

```math
\begin{cases}
\dot{S}(t) = -(1-u(t))\beta(t) S(t)I(t) \\
\dot{I}(t) = (1-u(t))\beta(t) S(t)I(t) - \gamma I(t)
\end{cases}
```

where _S_, _I_ are the susceptible and infected population fractions, _u_ ∈ [0, _ū_] is the intervention intensity (lockdown), _β_ is the transmission rate, and _γ_ is the recovery rate.

---

### Transmission Scenarios

We analyze two transmission cases:

**1. Constant _β_:**

```math
\beta(t) = \beta_0 \quad \forall t \geq 0
```

**2. Piecewise Constant _β_:**

```math
\beta(t) = \begin{cases}
\beta_1 & \text{if } 0 \leq t < T_c \\
\beta_2 & \text{if } t \geq T_c
\end{cases}
```

where _T_c is the time when the transmission rate _β_ changes.

---

## Optimal Control Formulation

### Objective

Our goal: **Minimize the final epidemic size** = maximize _S_(∞).

```math
J(u) = S(\infty) = \lim_{t \to \infty} S(t)
```

> **Note:** We deal with a non-standard cost function.

---

### Constraints

We impose an **L¹ budget constraint** on the intervention:

```math
\|u(\cdot)\|_{L^1} = \int_0^{\infty} u(t) \, dt \leq K
```

where _K_ > 0 is the total budget and _ū_ is the control upper bound.

> This approach extends optimal control problems by requiring that each intervention occur within an interval of the form [_t_, _t_ + _δ_], where _t_ is a decision variable.

Precisely, we consider the following set of admissible controls:

```math
\mathcal{U} = \left\{u : [0,\infty) \to [0,\bar{u}] \,\Big|\, u \text{ is measurable and } \|u\|_{L^1} \leq K \right\}
```

---

## Reproducibility

```@raw html
<details><summary>The documentation of this package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/control-toolbox/" *
                name *
                ".jl/tree/gh-pages/v" *
                version *
                "/assets/Manifest.toml"
link_project = "https://github.com/control-toolbox/" *
               name *
               ".jl/tree/gh-pages/v" *
               version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```