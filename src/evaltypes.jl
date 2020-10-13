abstract type EvalType end

struct Ball <: EvalType end

struct Taylor <: EvalType end

struct Asymptotic <: EvalType end

struct AsymptoticExpansion <: EvalType end

struct Symbolic <: EvalType end

export Ball, Asymptotic, AsymptoticExpansion, Symbolic
