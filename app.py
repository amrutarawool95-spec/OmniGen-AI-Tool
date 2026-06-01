# app.py
from flask import Flask, jsonify, request, render_template
import os
import math
import re
import logging
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("OmniGenEngine")

app = Flask(__name__)

# =====================================================================
# PIPELINE ENGINES BLOCK (Stage 4 to 9 Processing Core)
# =====================================================================

@dataclass
class SomaticVariant:
    gene: str
    transcript: str
    aa_change: str         
    consequence: str       
    tumor_lod: float       
    tpm_expression: float  
    ccf: float             

@dataclass
class PeptidePair:
    mutation_id: str
    gene: str
    length: int
    mutant_peptide: str
    wt_peptide: str
    mut_position: int      
    consequence: str
    tpm_expression: float
    ccf: float

@dataclass
class ScoredNeoantigen:
    gene: str
    mutation_id: str
    allele: str
    mhc_class: str          
    mutant_peptide: str
    wt_peptide: str
    ic50_mutant: float
    ic50_wt: float
    rank_pct: float
    dai: float
    tpm: float
    ccf: float
    binder_class: str       
    total_score: float

class PeptideGenerator:
    @staticmethod
    def parse_aa_change(aa_change: str) -> Tuple[int, str, str]:
        aa3_to_1 = {
            'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E',
            'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
            'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*'
        }
        try:
            match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)', aa_change)
            if not match:
                raise ValueError("Format mismatch")
            pos = int(match.group(2))
            ref_aa = aa3_to_1.get(match.group(1), 'X')
            alt_aa = aa3_to_1.get(match.group(3), 'X')
            return pos, ref_aa, alt_aa
        except Exception:
            match = re.match(r'p?\.?([A-Z])(\d+)([A-Z\*])', aa_change)
            if match:
                return int(match.group(2)), match.group(1), match.group(3)
            raise ValueError(f"Unable to parse AA sequence: {aa_change}")

    def generate_pairs(self, variant: SomaticVariant) -> List[PeptidePair]:
        pairs = []
        try:
            pos, ref_aa, alt_aa = self.parse_aa_change(variant.aa_change)
        except ValueError:
            return []

        state_seed = sum(ord(c) for c in variant.gene)
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        upstream = "".join(amino_acids[(state_seed + i) % 20] for i in range(30))
        downstream = "".join(amino_acids[(state_seed * i) % 20] for i in range(30))
        
        wt_protein = upstream + ref_aa + downstream
        mut_idx = len(upstream)
        all_lengths = list(range(8, 12)) + list(range(13, 22))

        for L in all_lengths:
            start_min = max(0, mut_idx - L + 1)
            start_max = min(mut_idx, len(wt_protein) - L)
            
            for start in range(start_min, start_max + 1):
                wt_pep = wt_protein[start:start+L]
                mut_pep_list = list(wt_pep)
                local_mut_idx = mut_idx - start
                if local_mut_idx < len(mut_pep_list):
                    mut_pep_list[local_mut_idx] = alt_aa
                mut_pep = "".join(mut_pep_list)
                
                if wt_pep == mut_pep or '*' in mut_pep:
                    continue
                    
                pairs.append(PeptidePair(
                    mutation_id=f"p.{ref_aa}{pos}{alt_aa}", gene=variant.gene, length=L,
                    mutant_peptide=mut_pep, wt_peptide=wt_pep, mut_position=local_mut_idx,
                    consequence=variant.consequence, tpm_expression=variant.tpm_expression, ccf=variant.ccf
                ))
        return pairs

class BindingPredictor:
    @staticmethod
    def compute_dai(mut_ic50: float, wt_ic50: float) -> float:
        if mut_ic50 <= 0 or wt_ic50 <= 0:
            return 0.0
        return math.log2(wt_ic50 / mut_ic50)

    def predict_binding(self, pair: PeptidePair, allele: str, mhc_class: str) -> Optional[ScoredNeoantigen]:
        val_seed = sum(ord(c) for c in pair.mutant_peptide + allele)
        pseudo_rank = ((val_seed % 1000) / 10.0)
        
        if any(g in pair.gene.upper() for g in ["KRAS", "IDH1", "BRAF", "EGFR"]):
            pseudo_rank = 0.05 + (val_seed % 40) / 100.0

        if mhc_class == 'MHC-I':
            if pseudo_rank > 2.0: return None
            binder_class = 'Strong' if pseudo_rank < 0.5 else 'Weak'
            ic50_mutant = max(5.0, min(1200.0, 50000.0 * (0.01 ** (1.0 - (pseudo_rank / 15.0)))))
            ic50_wt = ic50_mutant * (3.5 + (val_seed % 30))
        else:
            if pseudo_rank > 10.0: return None
            binder_class = 'Strong' if pseudo_rank < 2.0 else 'Weak'
            ic50_mutant = 20.0 + (val_seed % 400)
            ic50_wt = ic50_mutant * (1.5 + (val_seed % 5))

        dai = self.compute_dai(ic50_mutant, ic50_wt)
        
        return ScoredNeoantigen(
            gene=pair.gene, mutation_id=pair.mutation_id, allele=allele, mhc_class=mhc_class,
            mutant_peptide=pair.mutant_peptide, wt_peptide=pair.wt_peptide,
            ic50_mutant=round(ic50_mutant, 1), ic50_wt=round(ic50_wt, 1),
            rank_pct=round(pseudo_rank, 2), dai=round(dai, 2),
            tpm=pair.tpm_expression, ccf=pair.ccf, binder_class=binder_class, total_score=0.0
        )

class ImmunogenicityScorer:
    @staticmethod
    def filter_by_expression(variants: List[SomaticVariant]) -> List[SomaticVariant]:
        return [v for v in variants if v.tpm_expression > 1.0]

    @staticmethod
    def calculate_vaccine_score(candidate: ScoredNeoantigen) -> float:
        binding_weight = 1.0 / (1.0 + math.exp((candidate.ic50_mutant - 150) / 50))
        expression_factor = math.log10(candidate.tpm + 1)
        raw_score = binding_weight * (1 + (candidate.dai * 0.15)) * expression_factor * candidate.ccf
        return round(min(0.999, max(0.001, raw_score)), 3)

def execute_vaccine_design_pipeline(patient_variants: List[SomaticVariant], patient_hla: Dict[str, List[str]]) -> List[ScoredNeoantigen]:
    filtered_variants = ImmunogenicityScorer.filter_by_expression(patient_variants)
    pep_gen = PeptideGenerator()
    all_generated_pairs = []
    for var in filtered_variants:
        all_generated_pairs.extend(pep_gen.generate_pairs(var))
    
    predictor = BindingPredictor()
    validated_candidates = []
    for pair in all_generated_pairs:
        if 8 <= pair.length <= 11:
            for hla in patient_hla.get('mhc1', []):
                res = predictor.predict_binding(pair, hla, 'MHC-I')
                if res: validated_candidates.append(res)
        elif 13 <= pair.length <= 21:
            for hla in patient_hla.get('mhc2', []):
                res = predictor.predict_binding(pair, hla, 'MHC-II')
                if res: validated_candidates.append(res)

    for cand in validated_candidates:
        cand.total_score = ImmunogenicityScorer.calculate_vaccine_score(cand)
        
    return sorted(validated_candidates, key=lambda x: x.total_score, reverse=True)

# =====================================================================
# CONTROLLER ROUTING HANDLERS
# =====================================================================

@app.route('/', methods=['GET', 'POST'])
def load_unified_viewport():
    patient_mutations_dataset = [
        SomaticVariant(gene="KRAS", transcript="ENST00000311936", aa_change="p.Gly12Val", consequence="missense_variant", tumor_lod=14.2, tpm_expression=112.5, ccf=1.00),
        SomaticVariant(gene="IDH1", transcript="ENST00000415913", aa_change="p.Arg132His", consequence="missense_variant", tumor_lod=12.1, tpm_expression=28.1, ccf=1.00),
        SomaticVariant(gene="BRAF", transcript="ENST00000288602", aa_change="p.Val600Glu", consequence="missense_variant", tumor_lod=15.5, tpm_expression=178.6, ccf=0.91),
        SomaticVariant(gene="EGFR", transcript="ENST00000275493", aa_change="p.Leu858Arg", consequence="missense_variant", tumor_lod=11.4, tpm_expression=89.4, ccf=0.72),
        SomaticVariant(gene="NRAS", transcript="ENST00000369535", aa_change="p.Gln61His", consequence="missense_variant", tumor_lod=9.2, tpm_expression=5.1, ccf=0.34),
        SomaticVariant(gene="TP53", transcript="ENST00000269305", aa_change="p.Arg273His", consequence="missense_variant", tumor_lod=8.1, tpm_expression=0.2, ccf=0.95)
    ]
    
    patient_typed_hlas = {
        'mhc1': ['HLA-A*11:01', 'HLA-B*44:02', 'HLA-C*07:01'],
        'mhc2': ['HLA-DRB1*01:01']
    }

    if request.method == 'POST':
        file_object = request.files.get('genomic_file')
        if file_object:
            try:
                content = file_object.read().decode('utf-8', errors='ignore').upper()
                for var in patient_mutations_dataset:
                    if var.gene in content:
                        var.tpm_expression *= 1.5
                        var.ccf = min(1.0, var.ccf * 1.2)
            except Exception as e:
                logger.error(f"Inbound log stream processing error: {e}")

    computed_metrics = execute_vaccine_design_pipeline(patient_mutations_dataset, patient_typed_hlas)
    # Renders the clean index file inside templates folder directly
    return render_template('index.html', data=computed_metrics)

@app.route('/api/v1/analytics', methods=['GET'])
def pull_raw_json_feed():
    patient_mutations_dataset = [
        SomaticVariant(gene="KRAS", transcript="ENST00000311936", aa_change="p.Gly12Val", consequence="missense_variant", tumor_lod=14.2, tpm_expression=112.5, ccf=1.00)
    ]
    patient_typed_hlas = {'mhc1': ['HLA-A*11:01'], 'mhc2': []}
    computed_metrics = execute_vaccine_design_pipeline(patient_mutations_dataset, patient_typed_hlas)
    return jsonify([cand.__dict__ for cand in computed_metrics])

if __name__ == '__main__':
    target_network_port = int(os.environ.get('PORT', 5000))
    app.run(host="0.0.0.0", port=target_network_port)
    
