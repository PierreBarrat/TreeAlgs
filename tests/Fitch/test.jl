using Chain
using Test
using TreeAlgs.Fitch
using TreeTools


treefile = "Fitch/tree.nwk"
alnfile = "Fitch/aln.fasta"
t = @chain treefile read_tree(_) Fitch.fasta_to_tree(_, alnfile)

mich = t.lnodes["A/Michigan/41/2015"]
sing = t.lnodes["A/Singapore/H2013721a/2013"]
penn = t.lnodes["A/Pennsylvania/28/2014"]
thai = t.lnodes["A/Thailand/CU-CB166/2014"]
flo = t.lnodes["A/Florida/62/2015"]
n1 = lca(flo, thai)
n2 = lca(sing, penn)
r = t.root

@testset "Up" begin
    for i in 1:5
        Fitch.reset_fitch_states!(t)
        Fitch.set_pos!(i)

        #Init
        foreach(Fitch.init_fitchstate!, nodes(t))
        @test mich.data.state == [Set("G"), Set("C"), Set("C"), Set("G"), Set("-")][i]

        # up
        Fitch.fitch_up!(t.root)
        @testset "Fitch up" begin
            @test mich.data.state == [Set('G'), Set('C'), Set('C'), Set('G'), Set('-')][i]
            @test n1.data.state == [Set(['C', 'A']), Set('C'), Set(['A', 'G']), Set('G'), Set(['A', '-'])][i]
            @test n2.data.state == [Set(['C', 'A']), Set('C'), Set(['C', 'A']), Set('G'), Set(['A', '-'])][i]
            @test mich.anc.data.state == [Set(['G', 'C', 'A']), Set('C'), Set(['C', 'A', 'G']), Set('G'), Set(['-'])][i]
            @test t.root.data.state == [Set(['C', 'A']), Set(['C']), Set(['C', 'A']), Set(['G']), Set(['-'])][i]
        end
    end
end

begin
    Fitch.reset_fitch_states!(t)
    Fitch.fitch!(t)
    @testset "result" begin
        @test string(t.root.data.reconstructed_sequence...) in ["ACAG-", "ACCG-", "CCAG-", "CCCG-"]
        @test string(ancestor(mich).data.reconstructed_sequence...) in ["ACAG-", "ACCG-", "CCAG-", "CCCG-"]
        @test string(n2.data.reconstructed_sequence...) in ["ACAG-", "ACCG-", "CCAG-", "CCCG-"]
        @test in(
            string(n1.data.reconstructed_sequence...),
            ["ACAG-", "ACAG-", "ACCG-", "ACGG-", "CCAG-", "CCAG-", "CCCG-", "CCGG-"]
        )
    end
end
