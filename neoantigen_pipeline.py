#!/usr/bin/env python3
"""
Personalised Cancer Neoantigen Vaccine Designer - Simulation Engine Core
Based on the 9-Stage Computational Pipeline Technical Reference
"""

import math
import random
import logging
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

# Setup high-fidelity logging output
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("PipelineEngine")

# ==========================================
# DATA STRUCTURES & CONFIGURATION (Stages 2, 4, 5, 6, 7, 8)
# ==========================================

@dataclass
class SomaticVariant:
    gene: str
    transcript: str
    aa_change: str         # e.g., 'p.Val600Glu' (BRAF V600E) or 'p.Gly12Val' (KRAS G12V)
    consequence: str       # 'missense_variant', 'frameshift_variant'
    tumor_lod: float       # Mutect2 tumor Log-Odds score
    tpm_expression: float  # RNA-seq normalized expression (Transcripts Per Million)
    ccf: float             # Cancer Cell Fraction (PyClone-VI clonal estimate)

@dataclass
class PeptidePair:
    mutation_id: str
    gene: str
    length: int
    mutant_peptide: str
    wt_peptide: str
    mut_position: int      # 0-based index of mutation in the peptide
    consequence: str
    tpm_expression: float
    ccf: float

@dataclass
class ScoredNeoantigen:
    gene: str
    mutation_id: str
    allele: str
    mhc_class: str          # 'MHC-I' or 'MHC-II'
    mutant_peptide: str
    wt_peptide: str
    ic50_mutant: float
    ic50_wt: float
    rank_pct: float
    dai: float
    tpm: float
    ccf: float
    binder_class: str       # 'Strong' or 'Weak'
    total_score: float

# ==========================================
# PARSING & GENERATION MODULES (Stage 4)
# ==========================================

class PeptideGenerator:
    """Simulates Stage 4: Extracting flanking windows and creating Mutant/WT pairs."""
    
    @staticmethod
    def parse_aa_change(aa_change: str) -> Tuple[int, str, str]:
        """Parses Ensembl HGVSp mutations (e.g. p.Val600Glu -> (600, 'V', 'E'))"""
        aa3_to_1 = {
            'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E',
            'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
            'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*'
        }
        try:
            # Extract 3-letter codes and absolute amino acid position
            import re
            match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)', aa_change)
            if not match:
                raise ValueError("Format mismatch")
            pos = int(match.group(2))
            ref_aa = aa3_to_1.get(match.group(1), 'X')
            alt_aa = aa3_to_1.get(match.group(3), 'X')
            return pos, ref_aa, alt_aa
        except Exception:
            # Fallback for shorthand styles (e.g., V600E)
            import re
            match = re.match(r'p?\.?([A-Z])(\d+)([A-Z\*])', aa_change)
            if match:
                return int(match.group(2)), match.group(1), match.group(3)
            raise ValueError(f"Unable to safely parse AA change sequence: {aa_change}")

    def generate_pairs(self, variant: SomaticVariant) -> List[PeptidePair]:
        pairs = []
        try:
            pos, ref_aa, alt_aa = self.parse_aa_change(variant.aa_change)
        except ValueError as e:
            logger.warning(f"Skipping variant on gene {variant.gene}: {e}")
            return []

        # Simulated baseline protein context generation (representing local surrounding fasta database)
        # Creating a reproducible pseudo-sequence background unique to each gene string sequence
        random.seed(sum(ord(c) for c in variant.gene))
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        upstream = "".join(random.choices(amino_acids, k=30))
        downstream = "".join(random.choices(amino_acids, k=30))
        
        wt_protein = upstream + ref_aa + downstream
        mut_idx = len(upstream)  # Relative index inside our focal structural window

        # [span_7](start_span)Process lengths for MHC-I (8-11) and MHC-II (13-21)[span_7](end_span)
        mhc1_lengths = range(8, 12)
        mhc2_lengths = range(13, 22)
        all_lengths = list(mhc1_lengths) + list(mhc2_lengths)

        for L in all_lengths:
            # [span_8](start_span)Shift windows of length L over mutation index site[span_8](end_span)
            start_min = max(0, mut_idx - L + 1)
            start_max = min(mut_idx, len(wt_protein) - L)
            
            for start in range(start_min, start_max + 1):
                end = start + L
                wt_pep = wt_protein[start:end]
                
                # Insert the mutant residue variant
                mut_pep_list = list(wt_pep)
                local_mut_idx = mut_idx - start
                if local_mut_idx < len(mut_pep_list):
                    mut_pep_list[local_mut_idx] = alt_aa
                mut_pep = "".join(mut_pep_list)
                
                if wt_pep == mut_pep or '*' in mut_pep:
                    continue
                    
                pairs.append(PeptidePair(
                    mutation_id=f"p.{ref_aa}{pos}{alt_aa}",
                    gene=variant.gene,
                    length=L,
                    mutant_peptide=mut_pep,
                    wt_peptide=wt_pep,
                    mut_position=local_mut_idx,
                    consequence=variant.consequence,
                    tpm_expression=variant.tpm_expression,
                    ccf=variant.ccf
                ))
        return pairs

# ==========================================
# MATHEMATICAL PREDICTION ENGINES (Stages 5, 6)
# ==========================================

class BindingPredictor:
    """Simulates NetMHCpan-4.1 (MHC-I) & NetMHCIIpan-4.0 (MHC-II) analytical output profiles."""
    
    @staticmethod
    def compute_dai(mut_ic50: float, wt_ic50: float) -> float:
        [span_9](start_span)[span_10](start_span)"""DAI = log2(IC50_WT / IC50_mutant)[span_9](end_span)[span_10](end_span). Higher indicates enhanced mutant binding vs wildtype."""
        if mut_ic50 <= 0 or wt_ic50 <= 0:
            return 0.0
        return math.log2(wt_ic50 / mut_ic50)

    def predict_binding(self, pair: PeptidePair, allele: str, mhc_class: str) -> Optional[ScoredNeoantigen]:
        # Anchor matrix simulations based on common cancer variant features mentioned in documentation text:
        # e.g., KRAS G12V or IDH1 R132H driving highly specialized binding improvements
        random.seed(sum(ord(c) for c in pair.mutant_peptide + allele))
        
        if mhc_class == 'MHC-I':
            # [span_11](start_span)MHC-I Thresholds: Strong < 0.5% Rank, Weak < 2.0% Rank[span_11](end_span)
            rank_pct = random.uniform(0.01, 15.0)
            
            # Anchor structural optimization logic
            if "G12V" in pair.mutation_id or "R132H" in pair.mutation_id:
                rank_pct = random.uniform(0.02, 0.45)  # Force strong binder profile
            
            # Convert rank profile gracefully into IC50 configurations
            ic50_mutant = exp_ic50 = 50000.0 * (0.01 ** (1.0 - (rank_pct / 15.0)))
            ic50_mutant = max(1.0, min(ic50_mutant, 15000.0))
            
            # Back-calculate Wildtype baseline affinity to establish DAI variance
            if rank_pct < 0.5:
                binder_class = 'Strong'
                ic50_wt = ic50_mutant * random.uniform(5.0, 35.0) # Elevated WT IC50 = High DAI
            elif rank_pct < 2.0:
                binder_class = 'Weak'
                ic50_wt = ic50_mutant * random.uniform(1.1, 4.0)
            else:
                return None  # Discard non-binders directly during pipeline sweep
                
        else: # MHC-II
            # [span_12](start_span)[span_13](start_span)MHC-II Thresholds: Strong < 2.0% Rank, Weak < 10.0% Rank[span_12](end_span)[span_13](end_span)
            rank_pct = random.uniform(0.1, 25.0)
            if rank_pct < 2.0:
                binder_class = 'Strong'
                ic50_mutant = random.uniform(10.0, 150.0)
                ic50_wt = ic50_mutant * random.uniform(2.0, 10.0)
            elif rank_pct < 10.0:
                binder_class = 'Weak'
                ic50_mutant = random.uniform(151.0, 800.0)
                ic50_wt = ic50_mutant * random.uniform(1.0, 3.0)
            else:
                return None

        dai = self.compute_dai(ic50_mutant, ic50_wt)
        
        return ScoredNeoantigen(
            gene=pair.gene, mutation_id=pair.mutation_id, allele=allele, mhc_class=mhc_class,
            mutant_peptide=pair.mutant_peptide, wt_peptide=pair.wt_peptide,
            ic50_mutant=round(ic50_mutant, 2), ic50_wt=round(ic50_wt, 2),
            rank_pct=round(rank_pct, 3), dai=round(dai, 2),
            tpm=pair.tpm_expression, ccf=pair.ccf, binder_class=binder_class,
            total_score=0.0 # Will be derived inside the Stage 9 composite orchestrator
        )

# ==========================================
# EXPRESSION, CLONALITY & RANKING ENGINE (Stages 7, 8, 9)
# ==========================================

class ImmunogenicityScorer:
    """Calculates prioritized holistic rank metrics according to bio-expression parameters."""
    
    @staticmethod
    def filter_by_expression(variants: List[SomaticVariant]) -> List[SomaticVariant]:
        [span_14](start_span)[span_15](start_span)"""Stage 7 Biological Rationale: Filtering out genes with TPM <= 1[span_14](end_span)[span_15](end_span)."""
        return [v for v in variants if v.tpm_expression > 1.0]

    @staticmethod
    def calculate_vaccine_score(candidate: ScoredNeoantigen) -> float:
        """
        Calculates a matrix mathematical rank priority score:
        Combines Binding Affinity, Differential Agretopicity (DAI), RNA Expression (TPM), and Clonality (CCF).
        """
        # 1. Base affinity score component (Inversely proportional to IC50 mutant binding maximum)
        affinity_score = 100.0 * (1.0 - min(1.0, candidate.ic50_mutant / 500.0))
        
        # 2. Differential visibility score component (Higher DAI = newly visible to CTLs)
        dai_score = max(0.0, candidate.dai * 15.0)
        
        # 3. Logarithmic scalar optimization for expression level weightings
        expression_weight = math.log1p(candidate.tpm) / math.log1p(30.0) 
        
        # 4. Clonal weight linear modulation (CCF = 1.0 focuses attacks directly on corporate trunk lines)
        clonality_weight = candidate.ccf
        
        composite_score = (affinity_score + dai_score) * expression_weight * clonality_weight
        return round(composite_score, 2)

# ==========================================
# FULL INTEGRATED PIPELINE ORCHESTRATOR
# ==========================================

def execute_vaccine_design_pipeline(patient_variants: List[SomaticVariant], patient_hla: Dict[str, List[str]]) -> List[ScoredNeoantigen]:
    logger.info("Initializing Stage 1 & 2 Workflow Simulation: Reading patient data parameters...")
    [span_16](start_span)logger.info(f"Loaded {len(patient_variants)} raw somatic variant records called via Mutect2[span_16](end_span).")
    
    # Stage 7: Filter by verified RNA-seq trace visibility limits
    filtered_variants = ImmunogenicityScorer.filter_by_expression(patient_variants)
    [span_17](start_span)logger.info(f"Stage 7 Expression Filtering Complete: {len(filtered_variants)} / {len(patient_variants)} variants passed TPM > 1.0[span_17](end_span).")
    
    # Stage 4: Peptide Processing Loops
    pep_gen = PeptideGenerator()
    all_generated_pairs = []
    for var in filtered_variants:
        all_generated_pairs.extend(pep_gen.generate_pairs(var))
    [span_18](start_span)logger.info(f"Stage 4 Window Matrix Generation Complete: Extracted {len(all_generated_pairs)} overlapping candidate mutant/WT peptide pairs[span_18](end_span).")
    
    # Stage 5 & 6: Binding Predictions sweeps matching patient specific HLA complexes
    predictor = BindingPredictor()
    validated_candidates = []
    
    for pair in all_generated_pairs:
        # Check MHC-I alignments (Lengths 8-11)
        if 8 <= pair.length <= 11:
            for hla_allele in patient_hla.get('mhc1', []):
                res = predictor.predict_binding(pair, hla_allele, 'MHC-I')
                if res: validated_candidates.append(res)
                
        # Check MHC-II alignments (Lengths 13-21)
        elif 13 <= pair.length <= 21:
            for hla_allele in patient_hla.get('mhc2', []):
                res = predictor.predict_binding(pair, hla_allele, 'MHC-II')
                if res: validated_candidates.append(res)

    logger.info(f"Stage 5 & 6 Neural Affinity Scans Finished: Screened valid binder configurations.")
    
    # Stage 9: Integrated Multivariable Structural Ranking Sweep
    for candidate in validated_candidates:
        candidate.total_score = ImmunogenicityScorer.calculate_vaccine_score(candidate)
        
    # Order candidates cleanly descending based on programmatic mathematical score calculations
    final_ranked_vaccine_output = sorted(validated_candidates, key=lambda x: x.total_score, reverse=True)
    
    logger.info("Stage 9 Vaccine Matrix Compilation finalized successfully.")
    return final_ranked_vaccine_output

# ==========================================
# EXECUTION TEST RUN - PATIENT PROFILE PARAMETERS
# ==========================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("   PATIENT PERSONALIZED CANCER NEOANTIGEN VACCINE PIPELINE REPORT")
    print("="*70)
    
    # Shared input data values matching exact reference metrics specified in documentation
    patient_mutations_dataset = [
        SomaticVariant(gene="KRAS", transcript="ENST00000311936", aa_change="p.Gly12Val", consequence="missense_variant", tumor_lod=14.2, tpm_expression=45.2, ccf=1.00),
        SomaticVariant(gene="IDH1", transcript="ENST00000415913", aa_change="p.Arg132H", consequence="missense_variant", tumor_lod=12.1, tpm_expression=28.4, ccf=1.00),
        SomaticVariant(gene="TP53", transcript="ENST00000269305", aa_change="p.Arg248Trp", consequence="missense_variant", tumor_lod=9.8, tpm_expression=18.1, ccf=0.85),
        SomaticVariant(gene="BRAF", transcript="ENST00000288602", aa_change="p.Val600Glu", consequence="missense_variant", tumor_lod=15.5, tpm_expression=0.3, ccf=0.91), # Will get filtered (TPM < 1)
        SomaticVariant(gene="PIK3CA", transcript="ENST00000263967", aa_change="p.Glu545Lys", consequence="missense_variant", tumor_lod=8.4, tpm_expression=8.7, ccf=0.42)
    ]
    
    patient_typed_hlas = {
        'mhc1': ['HLA-A*02:01', 'HLA-B*44:02', 'HLA-C*07:01'], # Patient MHC-I alleles
        [span_19](start_span)'mhc2': ['HLA-DRB1*01:01', 'HLA-DRB1*15:01']          # Patient MHC-II alleles[span_19](end_span)
    }
    
    # Trigger full calculation run
    top_candidates = execute_vaccine_design_pipeline(patient_mutations_dataset, patient_typed_hlas)
    
    # Display the top 5 high-priority immunogenicity-ranked output selections
    print(f"\nSUCCESS: Pipeline finalized. Generated ranked top vaccine formulations:")
    print("-" * 115)
    print(f"{'Rank':<5}{'Gene':<8}{'Mutation':<10}{'HLA Allele':<14}{'Class':<8}{'Peptide (Mutant)':<16}{'IC50 (nM)':<11}{'DAI':<7}{'TPM':<7}{'CCF':<6}{'Score':<8}")
    print("-" * 115)
    
    for idx, cand in enumerate(top_candidates[:5], start=1):
        print(f"{idx:<5}{cand.gene:<8}{cand.mutation_id:<10}{cand.allele:<14}{cand.mhc_class:<8}"
              f"{cand.mutant_peptide[:14]+'...':<16}{cand.ic50_mutant:<11}{cand.dai:<7}{cand.tpm:<7}{cand.ccf:<6}{cand.total_score:<8}")
    print("-" * 115)
    print("Interpretation Guide: Rank 1 choices display strong binding features (low IC50), positive Differential Agretopicity (DAI), and high clonal presence (CCF=1.00).\n")
      
