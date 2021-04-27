module Fitch

using TreeTools
using BioSequences

import Base.isempty, Base.length, Base.intersect!, Base.union!

include("main.jl")
include("tools.jl")

export fitch!

end