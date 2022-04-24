export AbstractAnsatz, Ball, Asymptotic, AsymptoticExpansion, Symbolic

abstract type AbstractAnsatz{T} end

abstract type EvalType end

struct Ball <: EvalType end

struct Asymptotic <: EvalType end

struct AsymptoticExpansion <: EvalType end

struct Symbolic <: EvalType end
