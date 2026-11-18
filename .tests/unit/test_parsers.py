"""Unit tests for workflow/scripts/parsers.py. Fixtures are real strings seen in L. nagelii runs (BAGEL4 session, TYGS table).
Run: python -m pytest .tests/unit -q"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))
from parsers import accession_from_busco_json, acc_from_upload, cas_type, db_hits, html_text  # noqa: E402


# --- BAGEL4 'function' cells -------------------------------------------------
def test_db_hits_glued_pairs_from_aga58():
    cell = "bacteriocinII; Bacteriocin_II; L_biotic_typeA; Bacteriocin_IIc; 74.1;Mutacin_II(MutacinH-29B)166.2;Plantaricin_423"
    assert db_hits(cell) == [(74.1, "Mutacin_II(MutacinH-29B)"), (166.2, "Plantaricin_423")]


def test_db_hits_single_pair():
    assert db_hits("70.3;Helveticin-J") == [(70.3, "Helveticin-J")]


def test_db_hits_motif_only_has_no_scored_hit():
    assert db_hits("Lactococcin; L_biotic_typeA; Bacteriocin_IIc; ") == []


def test_db_hits_name_ending_in_integer_is_not_split():
    # 'Plantaricin_423' ends in digits without a decimal point: it is a name, not a score
    assert db_hits("166.2;Plantaricin_423") == [(166.2, "Plantaricin_423")]


# --- accession recovery -------------------------------------------------------
def test_acc_from_upload_restores_underscore():
    assert acc_from_upload("GCA061326205.1.fna") == "GCA_061326205.1"


def test_acc_from_upload_keeps_correct_name():
    assert acc_from_upload("GCA_061326205.1.fna") == "GCA_061326205.1"


def test_accession_from_busco_json_uses_file_name_not_folder():
    # regression: the folder name ('busco') was taken as the accession, so the focal genome looked missing
    assert accession_from_busco_json("results/qc/busco/GCA_020124705.2.json") == "GCA_020124705.2"


# --- Cas typing ---------------------------------------------------------------
def test_cas_type_ii_a_alone():
    prods = ["type II CRISPR RNA-guided endonuclease Cas9", "type II CRISPR-associated endonuclease Cas1",
             "CRISPR-associated endonuclease Cas2", "type II-A CRISPR-associated protein Csn2"]
    assert cas_type(prods) == "II-A"


def test_cas_type_i_b():
    prods = ["CRISPR-associated helicase/endonuclease Cas3", "type I-B CRISPR-associated protein Cas8b1/Cst1"]
    assert cas_type(prods) == "I-B"


def test_cas_type_mixed_i_b_and_iii_a():
    assert cas_type(["CRISPR-associated helicase Cas3", "type III-A CRISPR-associated protein Cas10/Csm1"]) == "I-B + III-A"


def test_cas_type_adaptation_genes_only_is_fragment():
    assert cas_type(["CRISPR-associated endonuclease Cas1", "CRISPR-associated endonuclease Cas2"]) == "fragment"


# --- TYGS cells ---------------------------------------------------------------
def test_html_text_genome_cell():
    assert html_text("<b>'GCA_047717745.1.fna'</b> ") == "GCA_047717745.1"


def test_html_text_type_strain_cell():
    cell = '<I>Liquorilactobacillus</I> <I>nagelii</I> <a target="_blank" href="https://www.dsmz.de/x">DSM 13675</a>'
    assert html_text(cell) == "Liquorilactobacillus nagelii DSM 13675"
