# 🧬 REE-Binding protein discovery report
*Report Date: 2026-01-15*

---

## Global Residue Conservation and Classification

### Methods 
The active site residues were identified by anchoring HMM match-states to structural coordinates of **XoxF1** (Lanthanide-dependent) and **MxaF** (Calcium-dependent) PQQ-dehydrogenases using hmmsearch and hmmalign. 

* hmmsearch to find all possible hits in the uniprot database
* hmmalign forces every discovery sequence to map directly onto the HMM’s match states. This creates a universal coordinate system ensuring that we are looking at the exact same physical pocket in the 3D structure across all sequences.
    * Calibration with high-resolution crystal structure to the HMM. By identifying which HMM match-state nodes correspond to the known metal-coordinating Aspartates ($D_{180}$ and $D_{182}$ in the reference), we define the "Search Window" for the rest of the dataset.
* Extract active site and classify based on known D-x-D (Lanthanide) and "D-x-[ATS]"

### Results
#### Classification table
| Classification | Count | % of Total |
| :--- | :--- | :--- |
| **Lanthanide (XoxF)** | 24959 | 42.2% |
| **Calcium (MxaF)** | 4108 | 6.9% |
| **Other / Truncated** | 30108 | 50.9% |

#### HMM match-state frequencies

| Rank | Amino Acid | Pos 1 (Anchor) | Pos 2 (Spacer) | Pos 3 (Switch) |
| :--- | :--- | :---: | :---: | :---: |
| 1 | **Y** | 0.06% | **44.2%** | 0.0% |
| 2 | **F** | 0.45% | **17.96%** | 0.01% |
| 3 | **M** | 0.0% | **13.4%** | 0.0% |
| 4 | **W** | 0.0% | **7.87%** | 0.01% |
| 5 | **L** | 0.01% | **5.86%** | 0.0% |
| 6 | **R** | 0.04% | **4.78%** | 0.02% |
| 7 | **A** | 0.07% | **1.75%** | 4.88% |
| 8 | **H** | 0.02% | **1.41%** | 0.01% |
| 9 | **V** | 0.11% | **0.76%** | 0.02% |
| 10 | **Q** | 0.02% | **0.47%** | 0.02% |
| 11 | **K** | 0.0% | **0.34%** | 0.25% |
| 12 | **N** | 4.2% | **0.23%** | 1.19% |
| 13 | **S** | 0.06% | **0.23%** | 2.85% |
| 14 | **I** | 0.02% | **0.18%** | 0.0% |
| 15 | **G** | 1.03% | **0.15%** | 0.39% |
| 16 | **T** | 0.34% | **0.11%** | 7.31% |
| 17 | **C** | 0.0% | **0.11%** | 0.0% |
| 18 | **D** | 93.19% | **0.07%** | 82.42% |
| 19 | **E** | 0.34% | **0.06%** | 0.57% |
| 20 | **P** | 0.02% | **0.02%** | 0.03% |


---

## 2. Taxonomic  // Environment Distribution

### Methods

* Uniprot identifyer --> ncbi taxid
    * fallback uniprot --> uniparc --> ncbi taxid 
* ncbi taxid --> taxonomy trough local taxonomy database. 
* ncbi taxid --> Omnicrobe semantic database
* To prevent pollution with identical or very similar sequences, used usearch to cluster to centroids (0.9 threshold). 

### Results

### Taxonomy

| Kingdom | Centroids | Percentage |
| :--- | :--- | :--- |
| **Unknown** | 7920 | 65.32% |
| **Bacteria** | 3946 | 32.54% |
| **Eukaryota** | 216 | 1.78% |
| **Archaea** | 43 | 0.35% |

**Archaea Highlights (43 centroids):**
- Halococcus dombrowskii
- Halolamina salifodinae
- Haloplanus vescus
- Halalkalicoccus paucihalophilus
- Natronococcus pandeyae
- Halopiger xanaduensis
- Natrialba asiatica
- Natrinema thermotolerans
- Halogranum amylolyticum
- Natronococcus occultus
- Haloterrigena salina
- Halorussus limi
- Halorussus gelatinilyticus
- Natrialba aegyptia
- Natrinema halophilum
- Natrinema ejinorense
- Natrinema versiforme
- Natrialba chahannaoensis
- Natrinema zhouii
- Methanoculleus thermophilus
- Methanoculleus bourgensis
- Methanothrix harundinacea
- Candidatus Caldarchaeum subterraneum
- Methanoculleus taiwanensis
- uncultured marine thaumarchaeote KM3_90_H07
- Methanospirillum lacunae
- Candidatus Nitrosocosmicus arcticus
- Candidatus Nitrosocosmicus franklandus
- Candidatus Nitrosocosmicus oleophilus

**Eukaryota Highlights (216 centroids):**
- Geodia barretti
- Plasmodium yoelii
- Cyprideis torosa
- Potamilus streckersoni
- Arachis hypogaea
- Gossypium stocksii
- Oppiella nova
- Symbiodinium pilosum
- Symbiodinium necroappetens
- Owenia fusiformis
- Gossypium klotzschianum
- Mikania micrantha
- Vitis vinifera
- Dillenia turbinata
- Dipteronia dyeriana
- Herrania umbratica
- Jatropha curcas
- Perilla frutescens
- Hibiscus syriacus
- Arachis duranensis
- Gossypium anomalum
- Canavalia gladiata
- Nepenthes gracilis
- Solanum tuberosum
- Vigna unguiculata
- Theobroma cacao
- Cicer arietinum
- Hibiscus trionum
- Lupinus angustifolius
- Fagus sylvatica
- Carpinus fangiana
- Saponaria officinalis
- Manihot esculenta
- Camellia sinensis
- Kingdonia uniflora
- Citrus x changshan-huyou
- Dovyalis caffra
- Carnegiea gigantea
- Juglans regia
- Morella rubra
- Erythroxylum novogranatense
- Vigna radiata
- Pisum sativum
- Liquidambar formosana
- Actinidia rufa
- Thalictrum thalictroides
- Prunus avium
- Quercus lobata
- Aquilegia coerulea
- Nyssa sinensis
- Tetracentron sinense
- Gossypium barbadense
- Psophocarpus tetragonolobus
- Turnera subulata
- Mucuna pruriens
- Coptis chinensis
- Sphenostylis stenocarpa
- Quillaja saponaria
- Nelumbo nucifera
- Glycine soja
- Carya illinoinensis
- Glycine max
- Phaseolus vulgaris
- Daucus carota
- Citrus unshiu
- Castanea mollissima
- Medicago truncatula
- Tripterygium wilfordii
- Prunus armeniaca
- Coffea canephora
- Prunus dulcis
- Cucumis melo
- Linum trigynum
- Helianthus annuus
- Ambrosia artemisiifolia
- Rhamnella rubrinervis
- Acer yangbiense
- Gossypium laxum
- Cinnamomum micranthum
- Pyrus ussuriensis x Pyrus communis
- Vigna angularis
- Populus deltoides
- Drosophila rhopaloa
- Lactuca sativa
- Malus domestica
- Centaurea solstitialis
- Ziziphus jujuba
- Olea europaea
- Rubus argutus
- Sesamum radiatum
- Deinandra increscens
- Oldenlandia corymbosa
- Ficus carica
- Nicotiana sylvestris
- Chenopodium quinoa
- Ceratopteris richardii
- Fraxinus pennsylvanica
- Rosa chinensis
- Eucalyptus grandis
- Trifolium subterraneum
- Tagetes erecta
- Anisodus tanguticus
- Buddleja alternifolia
- Citrus x clementina
- Callicarpa americana
- Actinidia chinensis
- Cajanus cajan
- Solanum verrucosum
- Salvia splendens
- Adiantum capillus-veneris
- Trema orientale
- Acer saccharum
- Physcomitrium patens
- Rhododendron griersonianum
- Sesamum indicum
- Acer negundo
- Trifolium pratense
- Erythranthe guttata
- Dorcoceras hygrometricum
- Sesamum calycinum
- Rubroshorea leprosula
- Morus notabilis
- Sesamum angustifolium
- Artemisia annua
- Prunus yedoensis
- Capsicum baccatum
- Lithospermum erythrorhizon
- Genlisea aurea
- Quercus rubra
- Podospora australis

**Conclusions:**
* Most eukaryote sequences are likely contaminations.
* Interesting set of Archaea, food for a perspective paper.

### Semantic link from taxonomy to environment

The discovered sequences were mapped to broad ecological categories using OntoBiotope (OBT) leaf nodes.

| Habitat | Count | Proportion |
| :--- | :--- | :--- |
| _Not Mapped_ | 39749 | 67.2% |
| Other | 5705 | 9.6% |
| Host-associated | 4965 | 8.4% |
| Aquatic | 3754 | 6.3% |
| Terrestrial | 3344 | 5.7% |
| Anthropogenic | 1658 | 2.8% |
---

**Conclusion:**
* NOT just aquatic, but widespread.
* Caveat is that this is a SEMANTIC link. NOT based on where the sample was taken.
* Further analysis of sample origin could strengthen / weaken this result. 

## 3. Niche Enrichment Analysis (Z-Scores)

### Methods

* Deduplicated hits with useach (usearch11 -sortbylength A.fa -fastaout B.fa) && (usearch11 -cluster_fast B.fa -id 0.9 -centroids centroids.fa -uc clusters.uc)
* Pearson’s Chi-squared test on the contingency of habitat vs. amino acid spacer. 
* "Expected" frequency of each residue assuming a random distribution across environments.
* Calculated Z-scores ($Z = \frac{O-E}{\sqrt{E}}$) to identify residues that were statistically over-represented (🔥) or under-represented (❄️) in specific niches.  

### Results

**Chi-Squared p-value:** $3.2e^{-9}$ suggests that the choice of (Y/F/W) is not random but an evolutionary adaptation to specific ecological requirements.

*Values represent Z-scores: (Observed - Expected) / sqrt(Expected). 🔥 indicates enrichment (>2.0).*

| AA (Pos 2) | Aquatic | Host-associated | Terrestrial | Anthropogenic |
| :--- | :---: | :---: | :---: | :---: |
| **Y** | 0.1 | 0.0 | -0.1 | -0.1 | 
| **F** | 1.3 | 1.1 | ❄️ -2.7 | 1.0 | 
| **W** | ❄️ -2.9 | -1.4 | 🔥 4.9 | -1.8 | 
| **A** | 1.0 | -1.0 | 0.4 | -0.9 | 
| **G** | 0.1 | -1.2 | -0.7 | 🔥 3.3 | 


## Archaea vs. Bacteria second AA 
Een vergelijking van de 43 archeale centroids met de bacteriële meerderheid laat een fundamenteel verschil zien in de architectuur van de actieve site.

**Statistische Validatie:**
- **Chi-Squared p-waarde:** 1.02e-6

| Kenmerk | Bacteria | Archaea | Implicatie |
| :--- | :--- | :--- | :--- |
| **Tyrosine (Y) Dominantie** | 46.5% | **69.8%** | Sterkere focus op cofactor-stabiliteit? |
| **Fenylalanine (F)** | 12.8% | **0.0%** | Negatieve selectie in Archaea |
| **Histidine (H) Verrijking** | 1.3% | **9.3%** | ? |

**Conclusie:** De Lanthanide-afhankelijkheid in Archaea is geen kopie van het bacteriële systeem, maar een onafhankelijke evolutionaire lijn. Dat maakt de cultivatie van de eerdere Halo Archaea nog interessanter. 

## Phylogeny

### Methods

* Used usearch 0.9 centroids 
* Only selected those assigned a taxonomy 
* FAMSA to align
* FastTree to build tree

### Results

<img width="8928" height="6397" alt="_ff00ff" src="https://github.com/user-attachments/assets/0ff000b2-c73f-4b38-a81c-40d072aaae44" />

**Conclusion:**
* Eukarya not diffuse but focussed. 2 explana possible:
     * There is a eukarya specific Ln dependent gene.
     * The eukaryote cluster is a bacterial lineage living in and on Eukarya, and thus frequently contaminates Eukaryotic genomes.
* Archaea form 2 clusters.  



