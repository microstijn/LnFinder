using Revise
using LnFinder
using Test
using BioSequences

@testset "LnFinder Tests" begin

    @testset "Motif Search (PROSITE Logic)" begin
        # 1. Lanthanide Case (D-x-D)
        # aa"" macro creates a LongAA type
        ln_seq = aa"MKTADLDG" 
        @test classify_pqq_adh(ln_seq) == LnFinder.Lanthanide

        # 2. Calcium Case (D-x-[ATS])
        ca_seq = aa"MKTADLAG"
        @test classify_pqq_adh(ca_seq) == LnFinder.Calcium
        
        # 3. String Input Support (Wrapper test)
        str_seq = "MKTADLDG"
        @test classify_pqq_adh(str_seq) == LnFinder.Lanthanide
    end
    
    @testset "HMMER Parser" begin
        # (Same parser test as before...)
    dummy_content = """
            #                                                               --- full sequence ---
            # target_name        accession  query_name           accession    E-value  score  bias
            # ------------------- ---------- -------------------- ---------- --------- ------ -----
            Gene_A             -            XoxF                 -            1.2e-50  150.5   0.1
            Gene_B             -            ExaF                 -            0.001    10.5   0.0
            #
            """
        temp_file = "test_hmmer.tblout"
        open(temp_file, "w") do io
            write(io, dummy_content)
        end
        
        df = parse_hmmer_tblout(temp_file)
        @test df.target[1] == "Gene_A"
        rm(temp_file)
    end
end