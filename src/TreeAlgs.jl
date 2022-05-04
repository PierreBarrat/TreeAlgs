module TreeAlgs

include("sequences.jl")
export fasta2tree!

include("CompatibilityTree/CompatibilityTree.jl")
include("Evolve/Evolve.jl")
include("Fitch/Fitch.jl")
include("LBI/LBI.jl")

end # module
