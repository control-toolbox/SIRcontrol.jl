# On the problem of minimizing the epidemic final size for SIR model via social distancing

## 📘 Introduction

### **Context and Motivation**

We revisit the problem of minimizing the epidemic final size in the **SIR model** through **social distancing interventions**.

Traditionally, this problem assumes a fixed **interval structure** for the timing of interventions. In contrast, we investigate a more flexible approach: **$L^1$-constrained controls** that limit total intervention effort, rather than when and how it's applied.

Key insight:

> Even with this generalization, the optimal control still occurs over a **single time interval** when the transmission rate is constant, and over at most **two (separate) time intervals** if the transmission rate changes once.
---


## The SIR Model with Control 

### **System Dynamics**

We consider the **SIR model** with a control variable $u(t)$ (representing social distancing):

```math
\begin{cases}
\dot{S}(t) = -(1-u(t))\beta(t) S(t)I(t) \\[5pt]
\dot{I}(t) = (1-u(t))\beta(t) S(t)I(t) - \gamma I(t)
\end{cases}
```

* \(S(t), I(t)\): Susceptible and infected population fractions
* \(u(t) \in [0, \bar u]\): Intervention intensity (lockdown)
* \(\beta(t)\): Transmission rate
* \(\gamma\): Recovery rate
---

### **Transmission Scenarios**

We analyze two transmission models:

1. **Constant $\beta$**  

```math
   \beta(t) = \beta_0 \quad \forall t \geq 0
```

2. **Piecewise Constant $\beta$**

```math
   \beta(t) = \begin{cases}
   \beta_1 & \text{if } 0 \leq t < T_c \\
   \beta_2 & \text{if } t \geq T_c
   \end{cases}
```

where $T_c$ is the time when the transmission rate $\beta$ changes.
---

## Optimal Control Formulation

### **Objective**

Our goal: **Minimize the final epidemic size** = maximize $S(\infty)$.

```math
J(u) = S(\infty) = \lim_{t \to \infty} S(t)
```

> Note that we deal with a non-standard cost function.
---


### **Constraints**

We impose an **$L^1$ budget constraint** on the intervention:

```math
\|u(\cdot)\|_{L^1} = \int_0^{\infty} u(t) \, dt \leq K
```

* \(K > 0\): Total budget
* \(\bar{u}\): Control upper bound

> This approach extends optimal control problems by requiring that each intervention occur within an interval of the form $[t, t + \delta]$, where $t$ is a decision variable.

---


Precisely, we consider the following set of admissible controls:

```math
\mathcal{U} = \{u : [0,\infty) \to [0,\bar{u}] \,|\, u \text{ is measurable and } \|u\|_{L^1} \leq K \}
```
