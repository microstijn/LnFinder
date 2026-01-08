using BioSequences

@enum MotifType begin
    Lanthanide
    Calcium
    Unknown
end

# We use PROSITE notation as recommended in the docs.
# D-x-D corresponds to Asp-Any-Asp
const LN_PATTERN = prosite"D-x-D"
# D-x-[ATS] corresponds to Asp-Any-{Ala,Thr,Ser}
const CA_PATTERN = prosite"D-x-[ATS]"

"""
    classify_pqq_adh(seq::BioSequence)

Classifies a biological sequence using native BioSequences.jl PROSITE search.
Accepts any subtype of BioSequence (e.g., LongAA).
"""
function classify_pqq_adh(seq::BioSequence)
    # BioSequences.occursin is optimized for bit-parallel search
    if occursin(LN_PATTERN, seq)
        return Lanthanide
    elseif occursin(CA_PATTERN, seq)
        return Calcium
    else
        return Unknown
    end
end

"""
    classify_pqq_adh(seq::AbstractString)

Wrapper for strings. Converts to LongAA for efficient biological searching.
"""
function classify_pqq_adh(seq::AbstractString)
    # Convert to LongAA (LongSequence{AminoAcidAlphabet}) as defined in docs
    # Note: This assumes the string contains valid amino acid characters.
    return classify_pqq_adh(LongAA(seq))
end