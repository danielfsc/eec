import sys
from pathlib import Path
# from utils import parse_modma_filename, fix_integer_wraparound, jensen_shannon, cramers_v
import utils as ut
from vars import TASK_UNSPECIFIED, EXPECTED_SUBJECT_IDS_55
import numpy as np
import pandas as pd
from derive_temporal import derive_temporal_protocol
from temporal_protocol import TemporalProtocolConfig
from evidence_policy import EvidencePolicy
from analysis import AnalysisConfig
from cohort import CohortConfig
from qc import qc_distribution_report, assert_qc_thresholds_frozen
from qc_freeze import QCFreezeConfig
from run import run_main
import warnings
# Add the parent directory to the module search path
# sys.path.append(str(Path(__file__).resolve().parent.parent))

# Now import your module normally


from run_config import RunConfig
from typing import  Dict


def selftest_v24(verbose: bool = True) -> Dict[str, bool]:
    """Autotestes com respostas conhecidas, incluindo os novos testes B1-B11."""
    import tempfile
    res: Dict[str, bool] = {}
    cfg = RunConfig()
    fs = cfg.acquisition.fs

    # --- basicos herdados ---
    res["nome_com_task"] = ut.parse_modma_filename("02010002_still.txt") == ("02010002", "still")
    res["nome_sem_task"] = ut.parse_modma_filename("02010001.txt") == ("02010001", TASK_UNSPECIFIED)
    a, n = ut.fix_integer_wraparound(np.array([[10.0, 4294967295.0]]), 32)
    res["wraparound"] = bool(n == 1 and a[0, 1] == -1.0)
    res["js_identico_zero"] = abs(ut.jensen_shannon([.5, .5], [.5, .5])) < 1e-12
    res["js_disjunto_um"] = abs(ut.jensen_shannon([1, 1e-12], [1e-12, 1]) - 1.0) < 1e-3
    res["cramers_v_colinear"] = abs(ut.cramers_v(["a", "a", "b", "b"], [1, 1, 0, 0]) - 1.0) < 1e-9
    res["bh_fdr_monotono"] = bool(np.all(np.diff(np.sort(ut.bh_fdr([0.001, 0.02, 0.3, 0.7]))) >= 0))

    # # --- B2: verificacao empirica de fs ---
    t = np.arange(int(fs * 60)) / fs
    x50 = np.tile(1000 * (np.sin(2 * np.pi * 10 * t) + 3.0 * np.sin(2 * np.pi * 50 * t)), (3, 1))
    v50 = ut.verify_sampling_rate(x50, cfg.acquisition)
    res["B2_detecta_pico_50hz"] = bool(v50["fs_supported"]
                                       and abs(v50["line_peak_hz_under_assumed_fs"] - 50) <= 0.6)
    x43 = np.tile(1000 * (np.sin(2 * np.pi * 10 * t) + 3.0 * np.sin(2 * np.pi * 43 * t)), (3, 1))
    res["B2_rejeita_pico_inesperado"] = bool(not ut.verify_sampling_rate(x43,
                                                                     cfg.acquisition)["fs_supported"])
    res["B2_fs_declarado_suposicao"] = bool(cfg.acquisition.fs_is_assumption)

    # # --- A5/B3: acomodacao e referencia externa ---
    x = ut.synth_recording(fs, 700.0, settle_s=150.0, seed=1)
    filt, fdiag = ut.preprocess_continuous(x, cfg.acquisition, cfg.preproc)
    prof = ut.block_profile(filt, fs, cfg.temporal.block_seconds)
    ref = ut.cohort_terminal_profile([prof], cfg.temporal)
    res["B3_referencia_coorte_normalizada"] = bool(ref is not None
                                                   and abs(float(np.sum(ref)) - 1.0) < 1e-6)
    st = ut.subject_settling_time(prof, cfg.temporal, cohort_ref_profile=ref)
    res["A5_detecta_transitorio"] = bool(np.isfinite(st["settling_time_s"])
                                         and st["settling_time_s"] > 30.0)
    res["B3_reporta_referencia"] = st["reference_used"] == "cohort_terminal_blocks"
    res["B3_flag_deriva_global"] = "drift_spans_record" in st
    y = ut.synth_recording(fs, 700.0, settle_s=0.5, seed=2)
    fl2, _ = ut.preprocess_continuous(y, cfg.acquisition, cfg.preproc)
    st2 = ut.subject_settling_time(ut.block_profile(fl2, fs, cfg.temporal.block_seconds),
                                cfg.temporal, cohort_ref_profile=ref)
    res["A5_estavel_acomoda_cedo"] = bool(st2["settling_time_s"] <= st["settling_time_s"])

    # # --- B4: cap explicito ---
    sdf = pd.DataFrame({"subject_id": ["a", "b", "c"], "settling_time_s": [60.0, 100.0, 80.0],
                        "drift_spans_record": [False, False, True]})
    prot = derive_temporal_protocol(sdf, {"a": 1200.0, "b": 1200.0, "c": 700.0}, cfg.temporal)
    res["A5_janela_derivada"] = bool(prot["skip_seconds"] >= 80.0 and prot["window_seconds"] > 0)
    res["B4_expoe_cap"] = ("skip_capped" in prot and "window_capped" in prot
                           and "protocol_valid" in prot)
    res["B4_conta_deriva_global"] = bool(prot["n_drift_spans_record"] == 1)
    tcap = TemporalProtocolConfig(max_skip_seconds=40.0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        pcap = derive_temporal_protocol(sdf, {"a": 1200.0, "b": 1200.0, "c": 700.0}, tcap)
        res["B4_avisa_quando_trunca"] = bool(pcap["skip_capped"]
                                             and not pcap["protocol_valid"] and len(w) >= 1)
    ws = prot["sensitivity_windows"]
    res["A5_janelas_disjuntas"] = all(ws[i][0] + ws[i][1] <= ws[i + 1][0] + 1e-9
                                      for i in range(len(ws) - 1))
    res["A5_exclui_curto"] = bool(derive_temporal_protocol(
        sdf, {"a": 1200.0, "b": 1200.0, "c": 72.8}, cfg.temporal)["excluded_short"] == ["c"])

    # # --- B9: parser de sexo ---
    s_txt, r_txt = ut.parse_sex_column(pd.Series(["M", "F", "male", "Female"]))
    res["B9_texto"] = list(s_txt) == [1.0, 0.0, 1.0, 0.0]
    s_num, r_num = ut.parse_sex_column(pd.Series([1, 2, 1, 2]))
    res["B9_numerico_1M_2F"] = (list(s_num) == [1.0, 0.0, 1.0, 0.0]
                                and r_num["coding_detected"] == "numeric_1M_2F")
    res["B9_relata_nao_mapeado"] = ut.parse_sex_column(pd.Series(["?", "?"]))[1]["n_mapped"] == 0

    # # --- B7: concordancia com IC ---
    g1 = pd.Series({"a": 0.4, "b": -0.3, "c": 0.2, "d": -0.5})
    ag = ut.sign_agreement_with_ci(g1, g1, n_boot=200, seed=1)
    res["B7_concordancia_total"] = abs(ag["sign_agreement"] - 1.0) < 1e-12
    res["B7_tem_ic_e_p"] = ("ci95" in ag and np.isfinite(ag["p_binomial_vs_chance"]))
    res["B7_discordancia"] = abs(ut.sign_agreement_with_ci(g1, -g1, n_boot=200,seed=1)["sign_agreement"]) < 1e-12

    # # --- B1: estimadores pareados ---
    res["B1_reps_pareados"] = bool(AnalysisConfig().permutation_matched_repeats >= 1)
    res["B8_repeticoes_iguais"] = AnalysisConfig().sensitivity_outer_repeats is None

    # # --- B5: distribuicao de QC e trava de congelamento ---
    qdf = pd.DataFrame({"subject_id": list("abcde"),
                        "ocular_index": [0.1, 0.2, 0.3, 0.9, 0.95],
                        "muscle_ratio": [0.01, 0.02, 0.03, 0.04, 0.05],
                        "min_channel_corr": [0.9, 0.9, 0.9, 0.9, -0.9]})
    dist = qc_distribution_report(qdf, cfg.qc)
    res["B5_tabela_distribuicao"] = bool(len(dist) >= 3
                                         and "n_excluded_by_threshold" in dist.columns)
    row = dist.loc[dist["metric"] == "ocular_index"].iloc[0]
    res["B5_conta_exclusoes"] = int(row["n_excluded_by_threshold"]) == 2
    res["B5_trava_research"] = False
    try:
        assert_qc_thresholds_frozen(QCFreezeConfig(), EvidencePolicy(mode="research"))
    except RuntimeError:
        res["B5_trava_research"] = True
    res["B5_libera_congelado"] = True
    try:
        assert_qc_thresholds_frozen(QCFreezeConfig(thresholds_frozen=True,
                                                   freeze_date="2026-08-22",
                                                   freeze_source="qc_distribution.csv"),
                                    EvidencePolicy(mode="research"))
    except RuntimeError:
        res["B5_libera_congelado"] = False
    res["B5_exige_data"] = False
    try:
        QCFreezeConfig(thresholds_frozen=True)
    except ValueError:
        res["B5_exige_data"] = True

    # # --- B11: coorte canonica ---
    res["B11_55_ids_unicos"] = len(set(EXPECTED_SUBJECT_IDS_55)) == 55
    aud = ut.audit_expected_cohort(list(EXPECTED_SUBJECT_IDS_55[:54]), CohortConfig())
    res["B11_audita_faltantes"] = (len(aud["missing_ids"]) == 1 and not aud["unexpected_ids"])
    res["B11_threshold_configuravel"] = abs(AnalysisConfig().decision_threshold - 0.5) < 1e-12

    # # --- features e barreiras ---
    ep = ut.epoch_signal(ut.extract_window(filt, fs, 300.0, 240.0), fs, 4.0)
    f, psd = ut.compute_psd(ep, fs, cfg.features)
    good, diag = ut.epoch_rejection_mask(ep, f, psd, cfg.qc, cfg.features)
    res["qc_epoca_maioria_boa"] = bool(good.mean() > 0.5)
    arrs = ut.epoch_feature_arrays(ep[good][:6], fs, cfg.features, cfg.acquisition.channel_names)
    row_f = ut.aggregate_features(arrs, cfg.acquisition.channel_names, cfg.features)
    res["nove_features_primarias"] = sum(k.startswith("p__") for k in row_f) == 9
    res["alfa_detectado"] = bool(abs(np.nanmedian(arrs["alpha_peak_hz"]) - 10.0) <= 1.5)
    res["barreira_primaria"] = False
    try:
        ut.assert_primary_only(["p__rel_alpha", "e__faa"])
    except RuntimeError:
        res["barreira_primaria"] = True
    res["barreira_circular"] = False
    try:
        ut.assert_no_circular_features(["p__rel_alpha", "phq9"])
    except RuntimeError:
        res["barreira_circular"] = True

    # # --- execucao ponta a ponta em dados sinteticos ---
    with tempfile.TemporaryDirectory() as td:
        root = Path(td); 
        eegd = root / "eeg" 
        eegd.mkdir()
        # print(eegd);
        for i, sid in enumerate(["02010001", "02010002", "02020004", "02020007"]):
            np.savetxt(eegd / ("%s_still.txt" % sid),
                       ut.synth_recording(fs, 700.0, settle_s=120.0, seed=i).T, fmt="%d")
        np.savetxt(eegd / "02020018_still.txt", ut.synth_recording(fs, 72.8, seed=9).T, fmt="%d")
        pd.DataFrame({"subject id": [2010001, 2010002, 2020004, 2020007, 2020018],
                      "type": ["MDD", "MDD", "HC", "HC", "HC"],
                      "age": [30, 40, 35, 45, 50], "gender": [1, 2, 1, 2, 1],
                      "education（years）": [12, 16, 14, 18, 12]}).to_csv(root / "meta.csv",
                                                                        index=False)
        got = run_main(eegd, root / "meta.csv", root / "out",
                      selection_criterion="Autoteste sintetico com transitorio inicial.")
        s = got["summary"]
        res["run_gera_protocolo"] = s["temporal_protocol"]["rule"] == "derived_from_data"
        res["run_exclui_curto"] = s["consort"]["n_excluded_short_duration"] == 1
        res["run_bloqueia_integration"] = "ROC/AUC" in s["blocked_analyses"]
        res["run_grava_artefatos"] = (root / "out" / "manifest.json").exists()
        res["run_sensibilidade"] = len(got["features_by_window"]) >= 1
        res["B5_grava_distribuicao"] = (root / "out" / "02_quality"
                                        / "qc_distribution.csv").exists()
        res["B2_grava_verificacao_fs"] = (root / "out" / "02_quality"
                                          / "sampling_rate_verification.csv").exists()
        res["B6_falhas_por_estagio"] = "failures_by_stage" in s["consort"]
        res["B11_consort_audita_coorte"] = bool(
            s["consort"]["cohort_audit"]["expected_declared"])
        res["B9_sexo_numerico_no_pipeline"] = bool(
            got["features"]["sex"].notna().all()) if "sex" in got["features"] else False
    if verbose:
        for k, v in res.items():
            print("  [%s] %s" % ("OK " if v else "FALHA", k))
        print("RESULTADO: %d/%d" % (sum(res.values()), len(res)))
    return res
