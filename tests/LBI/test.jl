using Test
using TreeAlgs.LBI
using TreeTools



nwk = "((A:2,B:2)AB:2,(C:2,D:2)CD:2)ABCD"
τ = 0.5

@testset "LBI" begin
	t1 = copy(node2tree(TreeTools.parse_newick(nwk)), LBIData)
	lbi!(t1, τ; normalize = false)
	e = exp(-t1.lnodes["A"].tau/τ)
	@test isapprox(t1.lnodes["ABCD"].data.lbi, 2τ*((1-e) + 2e*(1-e)), rtol=0.1 )
	t1.lnodes["A"].data.alive = false
	@test isapprox(t1.lnodes["ABCD"].data.lbi, 2τ*((1-e) + 3/2*e*(1-e)), rtol=0.1 )
	@test isapprox(t1.lnodes["CD"].data.lbi, τ*(2(1-e) + 1-e^3), rtol=0.1)
end
