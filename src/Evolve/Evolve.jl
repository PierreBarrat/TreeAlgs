module Evolve

using BioSequences
using Chain
using StatsBase
using SubstitutionModels
using TreeTools

export evolve!
export compute_mutations!

include("main.jl")

end
