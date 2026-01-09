using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using LnFinder
using Test
using BioSequences

@testset "LnFinder Tests" begin

    @testset "Motif Search (PROSITE Logic)" begin
        # Lanthanide Case (D-x-D)
        # aa"" macro creates a LongAA type
        ln_seq = aa"MKTADLDG" 
        @test classify_pqq_adh(ln_seq) == LnFinder.Lanthanide

        # Calcium Case (D-x-[ATS])
        ca_seq = aa"MKTADLAG"
        @test classify_pqq_adh(ca_seq) == LnFinder.Calcium
        
        # String Input Support (Wrapper test)
        str_seq = "MKTADLDG"
        @test classify_pqq_adh(str_seq) == LnFinder.Lanthanide
    end


    @testset "Alignment Analysis" begin
        # Create a dummy "Aligned FASTA"
        # Notice the gaps (-) and that the motif aligns at column 6 (1-based)
        # Seq1: Ln (D-L-D) at col 6
        # Seq2: Ca (D-L-A) at col 6
        # Seq3: Noise (Has a "D-L-D" at col 1, but it doesn't align with the others!)
        aligned_fasta = """
        >Ln_Seq
        -----DLD-------
        >Ca_Seq
        -----DLA-------
        >Noise_Seq
        DLD--X-X-------
        """
        
        temp_file = "test_aligned.fasta"
        open(temp_file, "w") do io; write(io, aligned_fasta); end

        # Run analysis
        df, col_idx = analyze_alignment(temp_file,sample_size=500)
        
        # Check if it found the correct column (Column 6 is the 'D')
        @test col_idx == 6
        
        # Check Classifications
        # The Noise_Seq has "DLD" at col 1, but the CONSENSUS is at col 6.
        # At col 6, Noise_Seq has 'X', so it should be Unknown.
        @test df.motif[1] == LnFinder.Lanthanide
        @test df.motif[2] == LnFinder.Calcium
        @test df.motif[3] == LnFinder.Unknown # Correctly rejected the false positive!

        rm(temp_file)
    end
end