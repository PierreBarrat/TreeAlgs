module CompatibilityTree

using LightGraphs
using TreeTools

export max_compatibility_splits, compatibility_graph

include("objects.jl")
include("compatibility_graph.jl")
include("main.jl")




end # module
