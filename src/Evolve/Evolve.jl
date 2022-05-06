module Evolve

using BioSequences
using StatsBase
using SubstitutionModels
using TreeTools

export evolve!
export compute_mutations!

include("main.jl")

end
