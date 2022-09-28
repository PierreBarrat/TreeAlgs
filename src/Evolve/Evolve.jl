module Evolve

using BioSequences
using StatsBase
using SubstitutionModels
using TreeTools
using FASTX

export evolve!
export compute_mutations!
export write_seq2fasta

include("main.jl")

end
