# Controls
This is meant to be a module to facilitate controlling State Space models such as the following

$$
\begin{align}
\dot x &= Ax + Bu \\
y &= Cx + Du
\end{align}
$$

This uses my `Modeling.jl` package which contains the `StateSpace` type which is defined as follows
```julia
struct StateSpace
    A::Matrix{<:Real}
    B::Vector{<:Real}
    C::Matrix{<:Real}
    D::Real
end
```

This package introduces the `Controller` type, which is defined as follows
```julia
mutable struct Controller
    C::Matrix{<:Real}
    a::Matrix{<:Real}
    A::Matrix{<:Real}
    α::Matrix{<:Real}
    K::Matrix{<:Real}
    kr::Real
end
```
defining the controllability matrix `C` along with the feedback control gains `K` and the feedworward gain `kr`
