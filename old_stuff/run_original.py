from typing import Dict, Any
from pathlib import Path
from old_stuff.evidence_policy import EvidencePolicy
from old_stuff.run_config import RunConfig
from old_stuff.subject_recording import SubjectRecording
from old_stuff.derive_temporal import derive_temporal_protocol
from vars  import __version__, OPEN_LIMITATIONS, BLOCKED_ANALYSES, CHANGELOG_V24

from dataclasses import asdict
import old_stuff.qc as qc;
import old_stuff.utils as ut
import pandas as pd
import numpy as np
import json

def run_main(eeg_dir: str | Path, metadata_path: str | Path, output_dir: str | Path,
            selection_criterion: str, policy: EvidencePolicy = EvidencePolicy(),
            cfg: RunConfig = RunConfig()) -> Dict[str, Any]:
    """Executa o pipeline completo: integracao -> protocolo temporal -> QC -> features -> analise.

    Artefatos em quatro pastas: 01_integration, 02_quality, 03_analysis, 04_sensitivity.
    Falhas por arquivo NAO abortam o lote; vao para failures.csv, agora segregadas
    por estagio (B6). O modo research so e liberado com limiares de QC congelados (B5).
    """
    if not isinstance(selection_criterion, str) or len(selection_criterion.strip()) < 10:
        raise ValueError("selection_criterion obrigatorio (>=10 caracteres).")
    qc.assert_qc_thresholds_frozen(cfg.qc_freeze, policy)  # B5
    out = Path(output_dir)
    d_int, d_qual = out / "01_integration", out / "02_quality"
    d_ana, d_sens = out / "03_analysis", out / "04_sensitivity"
    for d in (d_int, d_qual, d_ana, d_sens):
        d.mkdir(parents=True, exist_ok=True)

    files, failures = ut.discover_eeg_files(eeg_dir)
    metadata = ut.load_modma_metadata(metadata_path)

    # --- Etapa 1: leitura, esquema, filtragem, QC e perfil temporal (cego ao rotulo) ---
    recs: Dict[str, SubjectRecording] = {}
    filtered_map: Dict[str, np.ndarray] = {}
    qc_rows, block_rows, fs_rows = [], [], []
    for fp in files:
        # print(['Arquivo', fp])
        try:
            rec = ut.load_modma_txt(fp, cfg.acquisition, fp.name, cfg.schema)
            if rec.subject_id in recs:
                raise ValueError("Sujeito duplicado: %s" % rec.subject_id)
            fsv = ut.verify_sampling_rate(rec.data_counts, cfg.acquisition)  # B2
            fsv["subject_id"] = rec.subject_id
            filt, fdiag = ut.preprocess_continuous(rec.data_counts, cfg.acquisition, cfg.preproc)
            prof = ut.block_profile(filt, cfg.acquisition.fs, cfg.temporal.block_seconds)
            prof["subject_id"] = rec.subject_id
            recs[rec.subject_id] = rec
            block_rows.append(prof)
            filtered_map[rec.subject_id] = filt
            # print(['Entrada da função',rec, filt, cfg.acquisition, cfg.features, fdiag])
            # print(['resultado da subeject_qc_metrics',ut.subject_qc_metrics(rec, filt, cfg.acquisition, cfg.features, fdiag)])
            qc_rows.append(ut.subject_qc_metrics(rec, filt, cfg.acquisition, cfg.features, fdiag))
            # print(prof)
            # print(block_rows)
            fs_rows.append(fsv)
        except Exception as exc:
            # print(['Erro em alguma função', exc])
            failures.append({"file": fp.name, "subject_id": fp.stem.split("_")[0],
                             "stage": "read_or_profile",
                             "error": "%s: %s" % (type(exc).__name__, str(exc)[:200])})
    if not recs:
        raise RuntimeError("Nenhum registro utilizavel.")

    # B3: referencia espectral EXTERNA, mediana da coorte nos blocos terminais
    # print(block_rows)
    # print(['Block_rows-antes do cohort', block_rows])
    # print(['qc_rows-antes do cohort', qc_rows])
    cohort_ref = ut.cohort_terminal_profile(block_rows, cfg.temporal)
    settle_rows = []
    # print(['qc_rows-antes do for', qc_rows])
    for prof in block_rows:
        sid = str(prof["subject_id"].iloc[0])
        st = ut.subject_settling_time(prof, cfg.temporal, cohort_ref_profile=cohort_ref)
        st["subject_id"] = sid
        st["duration_s"] = recs[sid].parse_report["duration_s"]
        settle_rows.append(st)
    # print(['qc_rows - depois do for',qc_rows])
    qc_df = qc.cohort_qc_decision(pd.DataFrame(qc_rows), cfg.qc)
    qc_dist = qc.qc_distribution_report(qc_df, cfg.qc)                      # B5
    settle_df = pd.DataFrame(settle_rows)
    blocks_df = pd.concat(block_rows, ignore_index=True)
    fs_df = pd.DataFrame(fs_rows)
    print(['FS_DF',fs_df])

    # --- Etapa 2 (A5/B3/B4): derivacao do protocolo temporal ---
    durations = {s: r.parse_report["duration_s"] for s, r in recs.items()}
    protocol = derive_temporal_protocol(settle_df, durations, cfg.temporal)
    windows = protocol["sensitivity_windows"]
    primary_start, primary_len = windows[0]
    # --- Etapa 3: features por janela ---
    feat_by_window: Dict[str, pd.DataFrame] = {}
    for (s0, L) in windows:
        rows = []
        for sid, rec in recs.items():
            if sid not in protocol["eligible_subjects"]:
                continue
            try:
                rows.append(ut.extract_features_for_window(filtered_map[sid], rec, cfg, s0, L))
                # print(rows)
            except Exception as exc:
                failures.append({"file": rec.source_name, "subject_id": sid,
                                 "stage": "features_w%.0f" % s0,
                                 "error": "%s: %s" % (type(exc).__name__, str(exc)[:200])})
        feat_by_window["w_%.0f_%.0f" % (s0, L)] = pd.DataFrame(rows)

    key_primary = "w_%.0f_%.0f" % (primary_start, primary_len)
    features = feat_by_window[key_primary]
    print(['features',features])
    # --- Etapa 4: pareamento e coorte analitica ---
    merged = features.merge(metadata, on="subject_id", how="left")
    unmatched = merged.loc[merged["label"].isna(), "subject_id"].tolist()
    merged = merged.loc[merged["label"].notna()].copy()
    merged["label"] = merged["label"].astype(int)
    print(qc_df)
    qc_pass = set(qc_df.loc[qc_df["qc_pass"], "subject_id"])
    merged["qc_pass"] = merged["subject_id"].isin(qc_pass)
    drift_ok = ((pd.to_numeric(merged.get("drift_log2_rms"), errors="coerce")
                 <= cfg.temporal.max_within_window_drift_log2)
                & (pd.to_numeric(merged.get("drift_js"), errors="coerce")
                   <= cfg.temporal.max_within_window_js))
    merged["drift_pass"] = drift_ok.fillna(False)
    analytic = merged.loc[merged["qc_pass"] & merged["drift_pass"]].reset_index(drop=True)

    excl = merged.loc[~(merged["qc_pass"] & merged["drift_pass"])]
    tab = pd.crosstab(merged["label"], merged["qc_pass"] & merged["drift_pass"])
    fisher_p = (float(stats.fisher_exact(tab.to_numpy())[1])
                if tab.shape == (2, 2) else float("nan"))

    # B6: falhas segregadas por estagio; a v23 somava leitura e janela no mesmo total,
    # de modo que um sujeito podia ser contado ate 3 vezes como "arquivo que falhou".
    fail_df = pd.DataFrame(failures)
    print(['fail_df',fail_df])
    if not fail_df.empty and "stage" in fail_df.columns:
        by_stage = fail_df.groupby("stage").size().to_dict()
        subj_by_stage = {k: int(v["subject_id"].nunique())
                         for k, v in fail_df.groupby("stage")} if "subject_id" in fail_df else {}
        n_read_fail = int(fail_df.loc[fail_df["stage"].isin(
            ["discovery", "read_or_profile"])].shape[0])
    else:
        by_stage, subj_by_stage, n_read_fail = {}, {}, 0

    consort = {"n_files_found": len(files),
               "n_failures_total_events": int(len(failures)),
               "n_failures_read_stage": n_read_fail,
               "failures_by_stage": by_stage,
               "n_unique_subjects_failed_by_stage": subj_by_stage,
               "n_records_read": len(recs),
               "cohort_audit": ut.audit_expected_cohort(list(recs), cfg.cohort),
               "n_excluded_short_duration": len(protocol["excluded_short"]),
               "excluded_short_ids": protocol["excluded_short"],
               "n_with_features_primary_window": int(len(features)),
               "n_unmatched_metadata": len(unmatched), "unmatched_ids": unmatched,
               "n_after_qc_and_drift": int(len(analytic)),
               "n_excluded_qc_or_drift": int(len(excl)),
               "excluded_qc_ids": excl["subject_id"].tolist(),
               "differential_exclusion_fisher_p": fisher_p,
               "label_distribution_analytic": analytic["label"].value_counts().to_dict()}

    # --- Etapa 5: gravacao dos artefatos de integracao e qualidade ---
    fail_df.to_csv(d_int / "failures.csv", index=False)
    merged.to_csv(d_int / "features_primary_window.csv", index=False)
    for k, v in feat_by_window.items():
        v.to_csv(d_sens / ("features_%s.csv" % k), index=False)
    qc_df.to_csv(d_qual / "subject_qc.csv", index=False)
    qc_dist.to_csv(d_qual / "qc_distribution.csv", index=False)           # B5
    fs_df.to_csv(d_qual / "sampling_rate_verification.csv", index=False)  # B2
    settle_df.to_csv(d_qual / "temporal_settling.csv", index=False)
    blocks_df.to_csv(d_qual / "temporal_blocks.csv", index=False)

    fs_summary = {"fs_assumed_hz": float(cfg.acquisition.fs),
                  "fs_is_assumption": bool(cfg.acquisition.fs_is_assumption),
                  "n_subjects_supporting_fs": int(fs_df.get("fs_supported",
                                                            pd.Series(dtype=bool)).sum()),
                  "n_subjects_evaluated": int(len(fs_df)),
                  "status_counts": (fs_df["fs_verification"].value_counts().to_dict()
                                    if "fs_verification" in fs_df else {}),
                  "note": ("A ausencia de pico de rede NAO refuta fs; deixa a suposicao sem "
                           "corroboracao. Se fs estiver errada, todas as bandas e todo o eixo "
                           "temporal escalam por um fator constante.")}

    summary: Dict[str, Any] = {
        "version": __version__, "mode": policy.mode,
        "selection_criterion": selection_criterion,
        "config_fingerprint": cfg.fingerprint(),
        "sampling_rate_verification": fs_summary,
        "qc_thresholds_frozen": asdict(cfg.qc_freeze),
        "qc_threshold_impact": qc_dist.to_dict("records"),
        "temporal_protocol": protocol, "consort": consort,
        "open_limitations": list(OPEN_LIMITATIONS)}
    if not protocol.get("protocol_valid", True):
        summary["protocol_warning"] = (
            "B4: o inicio derivado foi truncado por max_skip_seconds; a janela pode "
            "conter transitorio de acomodacao. Resultados temporais nao sao validos "
            "sob a regra declarada.")

    # --- Etapa 6: analise (somente em modo research e com coorte suficiente) ---
    ok_cohort = (len(analytic) >= policy.min_subjects_for_ml
                 and analytic["label"].nunique() == 2
                 and analytic["label"].value_counts().min() >= policy.min_per_class_for_ml)
    if policy.mode == "integration" or not ok_cohort:
        summary["blocked_analyses"] = list(BLOCKED_ANALYSES)
        summary["block_reason"] = ("integration mode" if policy.mode == "integration"
                                   else "coorte insuficiente")
    else:
        acfg = cfg.analysis
        conf = ut.confirmatory_group_tests(analytic, acfg)
        conf.to_csv(d_ana / "confirmatory_tests.csv", index=False)
        adj = ut.confound_adjusted_models(analytic, acfg)
        adj.to_csv(d_ana / "confound_adjusted_models.csv", index=False)
        batch = ut.batch_confound_report(analytic, acfg)
        cols = ut.primary_columns(analytic)
        X = analytic[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        y = analytic["label"].to_numpy(int)
        res = ut.nested_cv_scores(X, y, acfg)
        met = ut.metrics_from_probs(y, res["prob_mean"], acfg,
                                 oof_per_repeat=res["oof_prob_per_repeat"])  # B10
        covs = [c for c in acfg.covariates if c in analytic.columns]
        base = {"majority_class": {"balanced_accuracy": 0.5, "roc_auc": 0.5}}
        if covs:
            Xd = analytic[covs].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            if np.isfinite(Xd).any():
                rd = ut.nested_cv_scores(Xd, y, acfg)
                base["demographic"] = ut.metrics_from_probs(y, rd["prob_mean"], acfg,
                                                         rd["oof_prob_per_repeat"])
                rc = ut.nested_cv_scores(np.hstack([X, Xd]), y, acfg)
                base["eeg_plus_demographic"] = ut.metrics_from_probs(
                    y, rc["prob_mean"], acfg, rc["oof_prob_per_repeat"])
            else:
                base["demographic"] = {"error": "covariaveis sem valores finitos (ver B9)"}
        # B1: observado e nulo sob o MESMO estimador
        obs_matched = ut.matched_observed_auc(X, y, acfg)
        perm = ut.permutation_auc_test(X, y, obs_matched, acfg)

        # B8: TODAS as janelas com o mesmo numero de repeticoes da primaria
        reps_sens = int(acfg.sensitivity_outer_repeats or acfg.outer_repeats)
        sens = {}
        for k, dfw in feat_by_window.items():
            m = dfw.merge(metadata, on="subject_id", how="inner")
            m = m.loc[m["subject_id"].isin(analytic["subject_id"])]
            if len(m) < policy.min_subjects_for_ml or m["label"].nunique() < 2:
                continue
            ck = ut.confirmatory_group_tests(m, acfg)
            Xk = m[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            yk = m["label"].to_numpy(int)
            rk = ut.nested_cv_scores(Xk, yk, acfg, repeats=reps_sens)
            sens[k] = {"n": int(len(m)), "cv_repeats": reps_sens,
                       "roc_auc": float(roc_auc_score(yk, rk["prob_mean"])),
                       "effects": ck.set_index("feature")["hedges_g"].to_dict(),
                       "n_significant_fdr": int(ck["significant_fdr"].sum())}
            ck.to_csv(d_sens / ("confirmatory_%s.csv" % k), index=False)

        # B7: concordancia com IC e p binomial, sem veredito qualitativo
        gref = conf.set_index("feature")["hedges_g"]
        concord = {}
        for k, v in sens.items():
            if k == key_primary:
                continue
            c = ut.sign_agreement_with_ci(gref, pd.Series(v["effects"]),
                                       n_boot=acfg.sign_agreement_n_boot,
                                       seed=acfg.random_seed)
            c["roc_auc"] = v["roc_auc"]
            concord[k] = c

        summary["analysis"] = {
            "primary_window": {"start_s": primary_start, "seconds": primary_len},
            "cv_repeats_all_windows": reps_sens,
            "classification_eeg_primary": met, "baselines": base,
            "permutation_test": perm,
            "n_significant_fdr": int(conf["significant_fdr"].sum()),
            "confound_adjusted_available": bool(len(adj)),
            "batch_confound": batch, "window_sensitivity": sens,
            "window_concordance": concord,
            "interpretation_guard": (
                "Metricas sao de validacao interna em um unico conjunto, com poder para "
                "detectar apenas efeitos grandes (g >= 0,78 com 26 vs 28). "
                + ("Lote e diagnostico sao colineares (V = 1,00): NAO interprete como "
                   "efeito clinico." if batch.get("batch_confounded") else ""))}
        json.dump(summary["analysis"], open(d_ana / "analysis_summary.json", "w"),
                  indent=2, ensure_ascii=False, default=ut.json_default)

    manifest = {"version": __version__, "changelog": list(CHANGELOG_V24),
                "environment": ut.environment_manifest(), "config": asdict(cfg),
                "policy": asdict(policy), "selection_criterion": selection_criterion,
                "input_files": [{"name": p.name, "sha256": ut.sha256_file(p)} for p in files]}
    json.dump(summary, open(out / "run_summary.json", "w"), indent=2,
              ensure_ascii=False, default=ut.json_default)
    json.dump(manifest, open(out / "manifest.json", "w"), indent=2,
              ensure_ascii=False, default=ut.json_default)
    return {"summary": summary, "manifest": manifest, "features": merged,
            "analytic": analytic, "qc": qc_df, "qc_distribution": qc_dist,
            "fs_verification": fs_df, "settling": settle_df,
            "features_by_window": feat_by_window, "failures": fail_df}