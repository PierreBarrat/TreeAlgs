module CompatibilityTree

using LightGraphs
using TreeTools

include("objects.jl")
include("compatibility_graph.jl")
include("main.jl")

export max_compatibility_splits, compatibility_graph



end # module
