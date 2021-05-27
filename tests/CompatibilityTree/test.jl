using Test
using TreeAlgs.CompatibilityTree
using LightGraphs
using TreeTools

 
#=
Example from Felsenstein's book
=#

dat = Array{Bool,2}(
		[1 0 0 1 1 0;
		0 0 1 0 0 0;
		1 1 0 0 0 0;
		1 1 0 1 1 1;
		0 0 1 1 1 0;
		0 0 0 0 0 0])

aln = [1 0 0 1 1 0;
		0 0 1 0 0 0;
		1 1 0 0 0 0;
		1 1 0 1 1 1;
		0 0 1 1 1 0;
		0 0 0 0 0 0]
labels = ["α", "β", "γ", "δ", "ε", "ω"]
tab = CompatibilityTree.CharTable(aln, labels)

gref1 = Graph(6) # Non-unique cols: 4 and 5 are interchangeable here (same col)
add_edge!(gref1, 1, 2)
add_edge!(gref1, 1, 3)
add_edge!(gref1, 1, 6)
add_edge!(gref1, 2, 3)
add_edge!(gref1, 2, 6)
add_edge!(gref1, 3, 6)
add_edge!(gref1, 4, 5)
add_edge!(gref1, 4, 6)
add_edge!(gref1, 5, 6) 

gref2 = Graph(5) # Unique columns
add_edge!(gref2, 1, 2)
add_edge!(gref2, 1, 3)
add_edge!(gref2, 1, 5)
add_edge!(gref2, 2, 3)
add_edge!(gref2, 2, 5)
add_edge!(gref2, 3, 5)
add_edge!(gref2, 4, 5)

Sref = TreeTools.SplitList(labels, 
	[Split([1,0,1,1,0,0]), Split([0,0,1,1,0,0]), Split([0,1,0,0,1,0]), Split([0,0,0,1,0,0])],
	ones(Bool, length(labels)), Dict{String,Split}())
sort!(Sref.splits, by=x->x.dat)

@testset "Felsenstein example (p93)" begin
	@test CompatibilityTree._compatibility_graph(Array{Bool,2}(dat')) == gref1
	@test CompatibilityTree.compatibility_graph(dat) == gref1
	@test CompatibilityTree.compatibility_graph(tab) == gref1
	@test max_compatibility_splits(tab) == Sref
end


#=
Test from a flu tree
=#
aln = [ 0  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  0  1  1  1  1
 1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  0  1  1  1  1  1  1  0  1  1  1  1  1  1  0  1  0  1  1  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  0  1  1  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
]
labels = ["A/Victoria/505/2013-egg", "A/PortoAlegre/LACENRS-275/2013", "A/Srilanka/68/2013", "NODE_58", "NODE_61", "NODE_54", "NODE_51"]
S = max_compatibility_splits(aln, labels)
@testset "Flu tree example" begin
	@test S.leaves[S.splits[1].dat] == ["A/Victoria/505/2013-egg", "NODE_54"]
	@test S.leaves[S.splits[2].dat] == ["A/Srilanka/68/2013", "A/Victoria/505/2013-egg", "NODE_54"]
	@test S.leaves[S.splits[3].dat] == ["A/PortoAlegre/LACENRS-275/2013", "NODE_61"]
end

