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

    @testset "FASTA Workflow (HMMER Web Integration)" begin
        # Create a dummy FASTA file representing HMMER output
        # Seq1: Ln-dependent (D-x-D)
        # Seq2: Ca-dependent (D-x-A)
        # Seq3: Garbage/Unknown
        dummy_fasta = """
        >seq_Ln
        MKTADLDG
        >seq_Ca
        MKTADLAG
        >seq_Unknown
        AAAAAA
        """
        
        # Write to a generic temp file
        temp_file = "test_hits.fasta"
        open(temp_file, "w") do io
            write(io, dummy_fasta)
        end

        # Run function
        df = classify_hits_from_fasta(temp_file)

        # Validate 
        @test size(df, 1) == 3
        
        # Check IDs
        @test df.id[1] == "seq_Ln"
        
        # Check Classifications
        @test df.motif[1] == LnFinder.Lanthanide    # Seq1 matches Ln
        @test df.motif[2] == LnFinder.Calcium       # Seq2 matches Ca
        @test df.motif[3] == LnFinder.Unknown       # Seq3 matches nothing

        # Cleanup
        rm(temp_file)
    end
end