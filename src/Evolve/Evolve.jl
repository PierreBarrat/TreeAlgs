module Evolve

using BioSequences
using StatsBase
using SubstitutionModels
using TreeTools
using FASTX
using Dates

export evolve!
export compute_mutations!
export write_seq2fasta

include("main.jl")

end
