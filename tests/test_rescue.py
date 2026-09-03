import pandas as pd
from hspc_response.rescue import compute_rescue_percent, classify_rescue


def test_compute_rescue_percent_full_rescue():
    # Gene fully rescued: GE moved it to log2FC=2, GE_ANAK brought it back to 0
    ge_results = pd.DataFrame({"log2FoldChange": [2.0]}, index=["gene_a"])
    ge_anak_results = pd.DataFrame({"log2FoldChange": [0.0]}, index=["gene_a"])

    rescue_pct = compute_rescue_percent(ge_results, ge_anak_results)
    assert rescue_pct["gene_a"] == 100.0


def test_compute_rescue_percent_no_rescue():
    # Gene not rescued at all: GE_ANAK stayed exactly where GE put it
    ge_results = pd.DataFrame({"log2FoldChange": [2.0]}, index=["gene_a"])
    ge_anak_results = pd.DataFrame({"log2FoldChange": [2.0]}, index=["gene_a"])

    rescue_pct = compute_rescue_percent(ge_results, ge_anak_results)
    assert rescue_pct["gene_a"] == 0.0


def test_compute_rescue_percent_partial_rescue():
    # Gene half-rescued: GE=2.0, GE_ANAK=1.0, so 50% rescue
    ge_results = pd.DataFrame({"log2FoldChange": [2.0]}, index=["gene_a"])
    ge_anak_results = pd.DataFrame({"log2FoldChange": [1.0]}, index=["gene_a"])

    rescue_pct = compute_rescue_percent(ge_results, ge_anak_results)
    assert rescue_pct["gene_a"] == 50.0


def test_classify_rescue_full_rescue():
    # GE_ANAK significantly differs from GE, and does NOT differ from RNP_NEG
    ge_anak_vs_ge = pd.DataFrame({"padj": [0.01]}, index=["gene_a"])
    ge_anak_vs_rnp = pd.DataFrame({"padj": [0.5]}, index=["gene_a"])

    classification = classify_rescue(ge_anak_vs_ge, ge_anak_vs_rnp, alpha=0.05)
    assert classification["gene_a"] == "full_rescue"


def test_classify_rescue_partial_rescue():
    # GE_ANAK significantly differs from GE, AND still differs from RNP_NEG
    ge_anak_vs_ge = pd.DataFrame({"padj": [0.01]}, index=["gene_a"])
    ge_anak_vs_rnp = pd.DataFrame({"padj": [0.01]}, index=["gene_a"])

    classification = classify_rescue(ge_anak_vs_ge, ge_anak_vs_rnp, alpha=0.05)
    assert classification["gene_a"] == "partial_rescue"


def test_classify_rescue_persistent():
    # GE_ANAK does NOT significantly differ from GE (no movement at all)
    ge_anak_vs_ge = pd.DataFrame({"padj": [0.5]}, index=["gene_a"])
    ge_anak_vs_rnp = pd.DataFrame({"padj": [0.01]}, index=["gene_a"])

    classification = classify_rescue(ge_anak_vs_ge, ge_anak_vs_rnp, alpha=0.05)
    assert classification["gene_a"] == "persistent"