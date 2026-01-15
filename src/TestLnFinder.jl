using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using LnFinder
using Test
using BioSequences
using DataFrames

@testset "LnFinder Tests" begin

    @testset "Motif Search (PROSITE Logic)" begin
        ln_seq = aa"MKTADLDG" 
        @test classify_pqq_adh(ln_seq) == LnFinder.Lanthanide

        ca_seq = aa"MKTADLAG"
        @test classify_pqq_adh(ca_seq) == LnFinder.Calcium
        
        unknown_seq = aa"MKTADLQG"
        @test classify_pqq_adh(unknown_seq) == LnFinder.Unknown

        str_seq = "MKTADLDG"
        @test classify_pqq_adh(str_seq) == LnFinder.Lanthanide
    end

    @testset "HMM Coordinate Logic" begin
        # Case A: Clean Sequence "MKDLD" -> D is at 3
        seq_clean = "MKDLD"
        @test identify_model_coordinate(seq_clean, "DLD") == 3

        # Case B: With inserts "MxxKxDLD" -> D is at 3 (inserts stripped)
        seq_messy = "MxxKxDLD" 
        @test identify_model_coordinate(seq_messy, "DLD") == 3
    end

    @testset "HMM Workflow Integration" begin
        # MOCK HMM ALIGNMENT
        # Ref 1: MxxK-DLD -> Clean: M, K, -, D, L, D
        # Indices: 1  2  3  4  5  6
        # Active site D is at Index 4.
        
        hmm_content = """
        >Ref_Ln (D-x-D)
        MxxK-DLD
        >Ref_Ca (D-x-A)
        MxxK-DLA
        >Test_Ln_Insertions
        MxxK-DLD
        >Test_Ca_Deletions
        M--DLA
        >Test_Gap_In_ActiveSite
        MxxK--LD
        >Test_Unknown_Residue
        MxxK-DLQ
        """
        # FIX EXPLANATION:
        # Test_Ca_Deletions changed from "M---DLA" to "M--DLA"
        # Columns: 1(M) 2(-) 3(-) 4(D) ... matches the Ref D-position at 4.
        
        temp_file = "test_hmm.afa"
        open(temp_file, "w") do io
            write(io, hmm_content)
        end

        try
            # Step A: Calibration
            idx = calibrate_active_site(temp_file)
            @test idx == 4
            
            # Step B: Parse
            df = parse_hmmalign_results(temp_file, idx)
            
            @test nrow(df) == 6
            
            # Ref 1 (Ln)
            @test df.classification[1] == :Lanthanide
            @test df.active_site[1] == "D-x-D"
            
            # Ref 2 (Ca)
            @test df.classification[2] == :Calcium
            @test df.active_site[2] == "D-x-A"
            
            # Test_Ln (Inserts ignored)
            @test df.classification[3] == :Lanthanide
            
            # Test_Ca (Deletions upstream should be ignored)
            @test df.classification[4] == :Calcium
            
            # Test_Gap (Missing anchor D)
            @test df.classification[5] == :Gap
            
            # Unknown
            @test df.classification[6] == :Unknown

        finally
            rm(temp_file, force=true)
        end
    end
end