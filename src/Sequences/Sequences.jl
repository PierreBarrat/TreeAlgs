module Sequences

using FASTX
using TreeTools

include("main.jl")
export fasta_to_tree!, write_fasta

include("mutations.jl")
import Base: ==, hash, reverse
export Mutation

end # module
