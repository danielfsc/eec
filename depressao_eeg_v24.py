#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""depressao_eeg_v24 - Pipeline MODMA 3-canais (Fp1, Fpz, Fp2).

NOTA DE PROCEDENCIA (v24.0.0)
-----------------------------
Deriva de depressao_eeg_v24.0.0. Implementa as sete correcoes B1-B9 apontadas na
auditoria metodologica da v23, todas verificadas contra o inventario factual dos
55 registros. Nenhuma alegacao empirica nova e introduzida sem verificacao.

B1 TESTE DE PERMUTACAO INVALIDO. Na v23 a AUC observada vinha de 20 repeticoes de
   CV e cada AUC nula de 1 repeticao: estimadores diferentes, permutabilidade
   quebrada, p-valor sem interpretacao e viesado na direcao anticonservadora.
   v24: observado e nulo usam O MESMO estimador (mesmo numero de repeticoes).
B2 fs=250 Hz DECLARADO COMO FATO. O descritor que sustenta fs tambem especifica
   repouso de ~5 min, contradito por registros de ~1200 s. v24: fs volta a ser
   SUPOSICAO (fs_is_assumption=True) e e verificada empiricamente pela posicao do
   pico de rede eletrica (verify_sampling_rate), com relatorio por sujeito.
B3 ESTIMADOR DE ACOMODACAO AUTO-REFERENCIAL. A referencia de estabilidade era a
   mediana da segunda metade do proprio registro; sob deriva monotona (31/55
   sujeitos com Spearman < -0,5) o alvo tambem deriva. v24: referencia ESPECTRAL
   EXTERNA (mediana da coorte nos blocos terminais) e sinalizador explicito de
   deriva que atravessa o registro inteiro (drift_spans_record).
B4 CAP SILENCIOSO DE max_skip_seconds. v24: skip_capped/window_capped expostos no
   protocolo, com aviso; protocolo marcado como nao valido quando o cap atua.
B5 LIMIARES DE QC NUNCA CONFRONTADOS COM A DISTRIBUICAO REAL. v24: relatorio de
   distribuicao empirica por metrica com contagem de exclusoes que cada limiar
   produziria, e TRAVA: modo research exige limiares congelados com data.
B6 CONTAGEM CONSORT AMBIGUA (falhas de leitura somadas as de janela). v24: falhas
   segregadas por estagio, com sujeitos unicos por estagio.
B7 VEREDITO DE ROBUSTEZ ARBITRARIO (sign_agreement >= 0,7 sobre 9 features ocorre
   por acaso em ~9% das vezes). v24: concordancia com IC bootstrap e p binomial
   contra o acaso; rotulo qualitativo removido.
B8 JANELAS NAO COMPARAVEIS (20 repeticoes na primaria, 5 nas de sensibilidade).
   v24: mesmo numero de repeticoes em todas as janelas, registrado no artefato.
B9 PARSING DE SEXO FRAGIL (codificacao numerica 1/2 virava NaN silencioso e
   apagava os modelos ajustados). v24: parser explicito com relatorio e aviso.
Correcoes menores: decision_threshold reintroduzido e efetivamente usado;
expected_subject_ids preenchido com a coorte canonica de 55 e auditado por
contagens tipo CONSORT; IC da AUC passa a incorporar a variabilidade da CV;
sm.tools.tools.np substituido por numpy.

CONTRATO DE EVIDENCIA
---------------------
mode="integration": nenhuma inferencia, nenhum p-valor, nenhuma metrica de
classificacao. mode="research": exige limiares de QC congelados e produz apenas
estimativas de VALIDACAO INTERNA em um unico conjunto. Nada aqui e evidencia
clinica, biomarcador ou diagnostico.

CONFUNDIMENTO ESTRUTURAL (nao corrigivel por codigo)
----------------------------------------------------
Prefixo de ID e diagnostico sao perfeitamente colineares (V de Cramer = 1,00):
0201 -> 26/26 MDD; 0202 e 0203 -> 29/29 HC. Lote de aquisicao e grupo clinico sao
inseparaveis. Poder estatistico com 26 vs 28: apenas efeitos g >= 0,78 sao
detectaveis com 80% de poder (g >= 1,02 apos correcao para 9 testes).

UNIDADES
--------
Os TXT contem contagens inteiras de A/D. O fator contas->uV e INFERIDO por ancora
de coorte apenas para relatorio descritivo. Nenhuma feature primaria depende dele.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
import warnings
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import butter, filtfilt, iirnotch, welch as _welch

try:
    from sklearn.dummy import DummyClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                                 brier_score_loss, confusion_matrix, f1_score,
                                 matthews_corrcoef, roc_auc_score)
    from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    _SKLEARN_OK, _SKLEARN_ERR = True, ""
except Exception as _exc:  # pragma: no cover
    _SKLEARN_OK, _SKLEARN_ERR = False, "%s: %s" % (type(_exc).__name__, _exc)

__version__ = "24.0.0"

CHANGELOG_V24 = (
    "B1 teste de permutacao com estimador identico entre observado e nulo.",
    "B2 fs volta a ser suposicao declarada e e verificada pelo pico de rede eletrica.",
    "B3 referencia de acomodacao espectral externa (coorte) + flag de deriva global.",
    "B4 skip_capped/window_capped expostos; protocolo invalidado quando o cap atua.",
    "B5 distribuicao empirica de QC publicada; modo research exige limiares congelados.",
    "B6 falhas segregadas por estagio nas contagens CONSORT.",
    "B7 concordancia entre janelas com IC bootstrap e p binomial; veredito removido.",
    "B8 mesmo numero de repeticoes de CV em todas as janelas.",
    "B9 parser de sexo robusto a codificacao numerica, com relatorio.",
    "B11 decision_threshold usado; coorte canonica de 55 IDs auditada; numpy direto.",
)

OPEN_LIMITATIONS = (
    "CONFUNDIMENTO ESTRUTURAL: prefixo de ID perfeitamente colinear com diagnostico "
    "(V de Cramer = 1,00; 0201=26/26 MDD, 0202+0203=29/29 HC). Efeito de lote e efeito "
    "clinico sao INSEPARAVEIS. Nenhum ajuste estatistico corrige isto; nenhum resultado "
    "deste pipeline pode ser atribuido a fisiologia sem replicacao com lotes cruzados.",
    "PODER INSUFICIENTE: com 26 MDD vs 28 HC, efeitos g < 0,78 nao sao detectaveis com "
    "80% de poder (g < 1,02 apos correcao para 9 testes). Efeitos tipicos de EEG de "
    "repouso na depressao (g ~ 0,3-0,5) sao indetectaveis; achados significativos aqui "
    "terao magnitude inflada (erro tipo M).",
    "REDUNDANCIA ESPACIAL: correlacao mediana entre canais = 0,983 (38/55 acima de 0,95). "
    "Fp1/Fpz/Fp2 carregam essencialmente um unico sinal; features por canal sao "
    "quase-redundantes e assimetria frontal nao e sustentavel neste dataset.",
    "Sem conjunto externo: todas as metricas sao de validacao interna.",
    "fs=250 Hz e SUPOSICAO verificada apenas indiretamente pelo pico de rede eletrica; "
    "o descritor que a sustenta e contradito pela duracao observada dos registros.",
    "Com 3 eletrodos frontais nao ha ICA identificavel; artefatos oculares sao "
    "quantificados e usados para exclusao, nao removidos por decomposicao.",
    "O fator contas->uV e INFERIDO por ancora de coorte, nao calibrado.",
    "A derivacao do protocolo temporal e a exclusao por outlier de coorte usam "
    "estatisticas de todos os sujeitos: sao etapas de COORTE, cegas ao desfecho, mas "
    "nao independentes entre sujeitos, e devem ser descritas como tal.",
    "A referencia de acomodacao para RMS permanece intra-sujeito (blocos terminais), "
    "porque a escala em contas nao e comparavel entre sujeitos: deriva que atravessa "
    "todo o registro e sinalizada, nao corrigida.",
)

EXPECTED_SUBJECT_IDS_55: Tuple[str, ...] = (
    "02010001", "02010002", "02010003", "02010005", "02010006", "02010007",
    "02010008", "02010009", "02010010", "02010011", "02010012", "02010013",
    "02010014", "02010015", "02010016", "02010018", "02010019", "02010020",
    "02010021", "02010022", "02010024", "02010025", "02010026", "02010030",
    "02010035", "02010038", "02020004", "02020007", "02020008", "02020010",
    "02020011", "02020014", "02020015", "02020016", "02020018", "02020019",
    "02020021", "02020022", "02020023", "02020025", "02020026", "02020027",
    "02030001", "02030002", "02030003", "02030004", "02030005", "02030006",
    "02030007", "02030008", "02030009", "02030010", "02030016", "02030020",
    "02030021",
)


QC_THRESHOLD_MAP: Dict[str, Tuple[str, str]] = {
    "line_noise_ratio_post": ("max_line_noise_ratio", "greater"),
    "ocular_index": ("max_ocular_index", "greater"),
    "muscle_ratio": ("max_muscle_ratio", "greater"),
    "zero_diff_fraction": ("max_zero_diff_fraction", "greater"),
    "saturation_fraction": ("max_saturation_fraction", "greater"),
    "nonstationarity_cv": ("max_nonstationarity_cv", "greater"),
    "min_channel_corr": ("min_channel_corr", "less"),
}


@dataclass(frozen=True)
class QCFreezeConfig:
    """B5: contrato de congelamento dos limiares de controle de qualidade.

    Limiares de QC nunca confrontados com a distribuicao empirica podem excluir
    zero ou quase toda a coorte sem que ninguem perceba antes de rodar. O
    procedimento obrigatorio e: (i) executar em modo integration; (ii) ler
    qc_distribution.csv, que informa quantos sujeitos CADA limiar excluiria;
    (iii) fixar os limiares; (iv) registrar aqui a data e a justificativa.
    Enquanto ``thresholds_frozen`` for False, o modo research fica bloqueado.
    """
    thresholds_frozen: bool = False
    freeze_date: Optional[str] = None
    freeze_source: str = ""
    require_freeze_for_research: bool = True

    def __post_init__(self) -> None:
        if self.thresholds_frozen and not self.freeze_date:
            raise ValueError("Limiares declarados congelados exigem freeze_date.")

# ==================================================================================
# BLOCO 1/11 - CONFIGURACAO IMUTAVEL
# ==================================================================================

@dataclass(frozen=True)
class AcquisitionConfig:
    """Parametros de aquisicao do experimento pervasivo de 3 eletrodos."""
    fs: float = 250.0
    fs_is_assumption: bool = True   # B2: revertido para suposicao declarada
    fs_verification_line_hz: Tuple[float, ...] = (50.0, 60.0)
    fs_verification_min_prominence_log10: float = 0.30
    fs_verification_tol_hz: float = 0.50
    fs_source: str = ("Descritor MODMA: experimento pervasivo de 3 eletrodos, A/D de "
                      "24 bits a 250 Hz; lote unico verificado por esquema.")
    protocol_id: str = "MODMA_3ch_pervasive_resting"
    channel_names: Tuple[str, ...] = ("Fp1", "Fpz", "Fp2")
    channel_order_source: str = "Ordem posicional das colunas do TXT (documentacao do dispositivo)."
    container_bits: int = 32
    adc_bits: int = 24
    lsb_to_uv: Optional[float] = None

    def __post_init__(self) -> None:
        if self.fs <= 0:
            raise ValueError("fs deve ser positivo.")
        if len(self.channel_names) != 3:
            raise ValueError("Esperados exatamente 3 canais.")
        if self.adc_bits > self.container_bits:
            raise ValueError("adc_bits nao pode exceder container_bits.")


@dataclass(frozen=True)
class SchemaConfig:
    """Contrato de esquema verificado em TODOS os arquivos antes de qualquer feature."""
    expected_n_channels: int = 3
    require_finite: bool = True
    require_integer: bool = True
    max_duration_s: float = 7200.0
    allowed_tasks: Tuple[str, ...] = ("still", "unspecified")
    enforce_allowed_tasks: bool = True


@dataclass(frozen=True)
class TemporalProtocolConfig:
    """A5: definicao do protocolo temporal DERIVADA DOS DADOS, cega ao rotulo.

    Motivacao empirica: o inventario dos 55 registros mostrou queda sistematica
    de amplitude nos primeiros minutos (acomodacao de eletrodo/impedancia e
    relaxamento do participante). Fixar a janela em 30-270 s, como nas v21/v22,
    amostra justamente esse transitorio.

    Procedimento (todas as etapas usam apenas sinais e tempos, nunca rotulos):
      1. perfil por blocos de ``block_seconds`` do registro filtrado;
      2. referencia de estado estavel = mediana dos blocos da metade final;
      3. tempo de acomodacao do sujeito = primeiro instante a partir do qual TODOS
         os blocos seguintes ficam dentro de ``settle_tol_log2`` em log2(RMS) e de
         ``settle_tol_js`` em divergencia espectral de Jensen-Shannon;
      4. acomodacao da coorte = quantil ``settle_quantile`` dos tempos individuais,
         arredondado para cima na grade de blocos, limitado por ``max_skip_seconds``;
      5. duracao da janela = maior valor comum a todos os sujeitos elegiveis,
         limitado por ``target_window_seconds`` e por ``max_window_seconds``;
      6. sensibilidade: ``n_sensitivity_windows`` janelas NAO SOBREPOSTAS de mesma
         duracao dentro da regiao comum.

    ``min_duration_seconds`` e criterio de elegibilidade PRE-DECLARADO. No lote de
    55 arquivos exatamente um registro (72,8 s) fica abaixo de qualquer limiar
    razoavel; sua exclusao e contada e testada quanto a diferenca entre grupos.
    """
    mode: str = "derived"                 # derived | fixed
    block_seconds: float = 20.0
    min_duration_seconds: float = 600.0
    settle_tol_log2: float = 0.35
    settle_tol_js: float = 0.15
    settle_quantile: float = 0.80
    min_skip_seconds: float = 30.0
    max_skip_seconds: float = 420.0
    target_window_seconds: float = 240.0
    max_window_seconds: float = 480.0
    n_sensitivity_windows: int = 3
    require_sensitivity: bool = True
    max_within_window_drift_log2: float = 1.50
    max_within_window_js: float = 0.35
    settling_ref_mode: str = "cohort"     # B3: cohort | self
    n_terminal_blocks: int = 5
    tail_trend_warn: float = 0.50
    fixed_skip_seconds: float = 30.0
    fixed_window_seconds: float = 240.0

    def __post_init__(self) -> None:
        if self.mode not in {"derived", "fixed"}:
            raise ValueError("mode temporal deve ser derived|fixed.")
        if self.block_seconds <= 0 or self.target_window_seconds <= 0:
            raise ValueError("Parametros temporais devem ser positivos.")
        if not (0.0 < self.settle_quantile <= 1.0):
            raise ValueError("settle_quantile fora de (0,1].")
        if self.min_skip_seconds > self.max_skip_seconds:
            raise ValueError("min_skip_seconds > max_skip_seconds.")
        if self.n_sensitivity_windows < 1:
            raise ValueError("n_sensitivity_windows >= 1.")


@dataclass(frozen=True)
class PreprocConfig:
    """Pre-processamento no sinal CONTINUO, antes de epocar."""
    notch_hz: Optional[float] = 50.0
    notch_q: float = 30.0
    notch_mode: str = "auto"
    notch_line_ratio_threshold: float = 0.02
    bandpass_hz: Tuple[float, float] = (1.0, 40.0)
    filter_order: int = 4
    edge_trim_seconds: float = 2.0
    epoch_seconds: float = 4.0
    epoch_overlap: float = 0.0

    def __post_init__(self) -> None:
        lo, hi = self.bandpass_hz
        if not (0 < lo < hi):
            raise ValueError("bandpass invalida.")
        if self.notch_mode not in {"auto", "always", "never"}:
            raise ValueError("notch_mode deve ser auto|always|never.")
        if not (0.0 <= self.epoch_overlap < 1.0):
            raise ValueError("epoch_overlap fora de [0,1).")


@dataclass(frozen=True)
class QCConfig:
    """QC em tres niveis: epoca, participante e coorte. Sempre cego ao rotulo."""
    p2p_ratio_max: float = 5.0
    p2p_ratio_min: float = 0.2
    jump_ratio_max: float = 8.0
    lf_share_ratio_max: float = 2.5
    max_rejected_fraction: float = 0.5
    max_line_noise_ratio: float = 0.25
    max_ocular_index: float = 0.80
    max_muscle_ratio: float = 0.60
    max_zero_diff_fraction: float = 0.20
    max_saturation_fraction: float = 1e-4
    max_nonstationarity_cv: float = 1.00
    min_channel_corr: float = -0.50
    cohort_mad_z_max: float = 4.0
    cohort_outlier_metrics: Tuple[str, ...] = (
        "line_noise_ratio_post", "ocular_index", "muscle_ratio",
        "nonstationarity_cv", "log_rms_counts")
    wrap_fraction_warn: float = 0.60


@dataclass(frozen=True)
class ScaleInferenceConfig:
    """Inferencia empirica de referencia e de escala (descritiva)."""
    enabled: bool = True
    anchor_rms_uv: float = 15.0
    anchor_source: str = "Valor central declarado a priori para RMS 1-40 Hz frontal de repouso."
    average_reference_tol: float = 1e-6
    common_mode_high: float = 0.50
    gcd_max_samples: int = 200000


@dataclass(frozen=True)
class FeatureConfig:
    """Parametros PRE-ESPECIFICADOS de extracao (a janela vem de TemporalProtocolConfig)."""
    welch_nperseg_seconds: float = 2.0
    welch_overlap: float = 0.5
    bands: Tuple[Tuple[str, float, float], ...] = (
        ("delta", 1.0, 4.0), ("theta", 4.0, 8.0), ("alpha", 8.0, 13.0), ("beta", 13.0, 30.0))
    total_band: Tuple[float, float] = (1.0, 40.0)
    alpha_peak_band: Tuple[float, float] = (7.0, 13.0)
    alpha_min_peak_log10: float = 0.05
    aperiodic_fit_band: Tuple[float, float] = (3.0, 35.0)
    aperiodic_peak_sd: float = 1.5
    aperiodic_max_iter: int = 5
    ocular_band: Tuple[float, float] = (0.5, 3.0)
    muscle_band: Tuple[float, float] = (20.0, 40.0)
    perm_entropy_order: int = 3
    perm_entropy_delay: int = 1
    higuchi_kmax: int = 10
    aggregate: str = "median"
    min_good_epochs: int = 30

    def __post_init__(self) -> None:
        if self.aggregate not in ("median", "mean"):
            raise ValueError("aggregate deve ser median|mean.")


@dataclass(frozen=True)
class CohortConfig:
    """Definicao operacional da coorte-alvo e contagens tipo CONSORT."""
    expected_n: Optional[int] = 55
    expected_subject_ids: Optional[Tuple[str, ...]] = EXPECTED_SUBJECT_IDS_55
    duplicate_policy: str = "prefer_task"
    preferred_task: str = "still"
    require_one_file_per_subject: bool = True


@dataclass(frozen=True)
class AnalysisConfig:
    """Plano analitico congelado ANTES de ver resultados."""
    covariates: Tuple[str, ...] = ("age", "sex", "education_years")
    fdr_alpha: float = 0.05
    n_boot_effect: int = 2000
    outer_folds: int = 5
    outer_repeats: int = 20
    inner_folds: int = 5
    logistic_c_grid: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0)
    inner_scoring: str = "roc_auc"
    n_boot_metrics: int = 2000
    n_permutations: int = 200
    permutation_outer_repeats: int = 1
    batch_collinearity_warn: float = 0.60
    decision_threshold: float = 0.50            # B11
    permutation_matched_repeats: int = 5        # B1: identico para observado e nulo
    sensitivity_outer_repeats: Optional[int] = None  # B8: None => outer_repeats
    sign_agreement_n_boot: int = 2000
    random_seed: int = 20240517


@dataclass(frozen=True)
class EvidencePolicy:
    """Contrato de evidencia. 'integration' nunca executa inferencia/ML."""
    mode: str = "integration"
    min_subjects_for_ml: int = 20
    min_per_class_for_ml: int = 20
    require_both_classes: bool = True
    forbid_clinical_scales: bool = True

    def __post_init__(self) -> None:
        if self.mode not in {"integration", "research"}:
            raise ValueError("mode deve ser 'integration' ou 'research'.")


@dataclass(frozen=True)
class RunConfig:
    """Configuracao raiz imutavel com impressao digital reproduzivel."""
    acquisition: AcquisitionConfig = field(default_factory=AcquisitionConfig)
    schema: SchemaConfig = field(default_factory=SchemaConfig)
    temporal: TemporalProtocolConfig = field(default_factory=TemporalProtocolConfig)
    preproc: PreprocConfig = field(default_factory=PreprocConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    qc_freeze: QCFreezeConfig = field(default_factory=QCFreezeConfig)
    scale: ScaleInferenceConfig = field(default_factory=ScaleInferenceConfig)
    cohort: CohortConfig = field(default_factory=CohortConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    random_seed: int = 20240517

    def __post_init__(self) -> None:
        nyq = self.acquisition.fs / 2.0
        lo, hi = self.preproc.bandpass_hz
        if hi >= nyq:
            raise ValueError("Banda alta >= Nyquist.")
        if self.preproc.notch_hz is not None and self.preproc.notch_hz >= nyq:
            raise ValueError("notch acima de Nyquist.")

    def fingerprint(self) -> str:
        return hashlib.sha256(
            json.dumps(asdict(self), sort_keys=True, default=str).encode()).hexdigest()[:16]


BLOCKED_ANALYSES = ["ROC/AUC", "cross-validation", "classification metrics",
                    "accuracy/sensitivity/specificity/F1", "p-values",
                    "between-group hypothesis tests", "permutation tests",
                    "clinical inference", "biomarker claims", "MDD detection claims"]

PRIMARY_FEATURES = ("rel_delta", "rel_theta", "rel_alpha", "rel_beta", "spec_entropy",
                    "alpha_peak_hz", "aperiodic_slope", "hjorth_mobility", "hjorth_complexity")

CLINICAL_SCALE_TOKENS = ("phq", "gad", "psqi", "ctq", "les", "ssrs")

# ==================================================================================
# BLOCO 2/11 - LEITURA, ESQUEMA E METADADOS
# ==================================================================================

_NAME_RE = re.compile(r"^(?P<sid>\d{8})(?:_(?P<task>[A-Za-z]+))?\.txt$")
TASK_UNSPECIFIED = "unspecified"


@dataclass(frozen=True)
class SubjectRecording:
    subject_id: str
    task: str
    data_counts: np.ndarray
    fs: float
    source_name: str
    parse_report: Dict[str, Any]


def parse_modma_filename(name: str) -> Tuple[str, str]:
    """Extrai (subject_id, task) do nome do arquivo; task ausente -> 'unspecified'."""
    m = _NAME_RE.match(name.strip())
    if not m:
        raise ValueError("Nome fora do padrao MODMA: %r" % name)
    return m.group("sid"), (m.group("task") or TASK_UNSPECIFIED)


def fix_integer_wraparound(arr: np.ndarray, container_bits: int) -> Tuple[np.ndarray, int]:
    """Reinterpreta inteiros sem sinal como complemento de dois (v > 2**31 -> v - 2**32)."""
    half, full = 2 ** (container_bits - 1), 2 ** container_bits
    mask = arr > half
    n = int(mask.sum())
    if n == 0:
        return arr, 0
    out = arr.copy()
    out[mask] = out[mask] - full
    return out, n


def load_modma_txt(path: str | Path, acq: AcquisitionConfig, source_name: str,
                   schema: SchemaConfig) -> SubjectRecording:
    """Le um TXT MODMA de 3 canais e valida o esquema antes de devolver o registro."""
    sid, task = parse_modma_filename(source_name)
    arr = np.loadtxt(str(path), dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != schema.expected_n_channels:
        raise ValueError("%s: esperado (N,%d), obtido %s."
                         % (source_name, schema.expected_n_channels, arr.shape))
    if schema.require_finite and not np.isfinite(arr).all():
        raise ValueError("%s: valores nao finitos." % source_name)
    if schema.require_integer and not np.allclose(arr, np.round(arr)):
        raise ValueError("%s: valores nao inteiros." % source_name)
    if schema.enforce_allowed_tasks and task not in schema.allowed_tasks:
        raise ValueError("%s: tarefa '%s' fora do protocolo." % (source_name, task))
    arr, n_wrap = fix_integer_wraparound(arr, acq.container_bits)
    x = np.ascontiguousarray(arr.T)
    dur = x.shape[1] / acq.fs
    if dur > schema.max_duration_s:
        raise ValueError("%s: duracao %.1f s acima do maximo." % (source_name, dur))
    rep = {"n_samples": int(x.shape[1]), "duration_s": float(dur),
           "n_wraparound_fixed": int(n_wrap),
           "wraparound_fraction": float(n_wrap / x.size),
           "dc_offset_counts": [float(v) for v in x.mean(axis=1)],
           "task_parsed": task}
    return SubjectRecording(sid, task, x, acq.fs, source_name, rep)


def canonical_subject_id(values) -> pd.Series:
    """Normaliza identificadores para 8 digitos com zeros a esquerda."""
    s = pd.Series(values).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(8)


def discover_eeg_files(eeg_dir: str | Path) -> Tuple[List[Path], List[Dict[str, str]]]:
    """Descobre TXT validos; nomes fora do padrao vao para quarentena (nao abortam)."""
    root = Path(eeg_dir)
    if not root.is_dir():
        raise FileNotFoundError("Diretorio EEG inexistente: %s" % root)
    valid, quarantined = [], []
    for p in sorted(root.glob("*.txt")):
        try:
            parse_modma_filename(p.name)
            valid.append(p)
        except ValueError as exc:
            quarantined.append({"file": p.name, "stage": "filename", "error": str(exc)[:200]})
    if not valid:
        raise FileNotFoundError("Nenhum TXT MODMA valido encontrado.")
    return valid, quarantined


def parse_sex_column(values: pd.Series) -> Tuple[pd.Series, Dict[str, Any]]:
    """B9: parser de sexo robusto a codificacao numerica, com relatorio explicito.

    ERRO CORRIGIDO. A v23 fazia ``str[0]`` e mapeava apenas 'M'/'F'. Se a planilha
    codificasse sexo como 1/2 - convencao comum, inclusive no MODMA - o resultado
    era NaN para TODOS os participantes. Em ``confound_adjusted_models`` o
    ``dropna()`` esvaziava o dataframe, ``len(sub) < 20`` disparava e cada feature
    era pulada SEM erro: os modelos ajustados por confundidores simplesmente nao
    rodavam, e nada no relatorio deixava isso evidente.

    Convencao de saida: 1.0 = masculino, 0.0 = feminino, NaN = nao mapeado.
    A codificacao numerica assumida (1=M, 2=F) e declarada no relatorio para que
    possa ser contestada; se estiver invertida, o sinal do coeficiente de sexo
    inverte, sem afetar as demais covariaveis.
    """
    raw = pd.Series(values)
    s = raw.astype(str).str.strip().str.upper()
    out = pd.Series(np.nan, index=raw.index, dtype=float)
    out[s.str.startswith("M") | s.isin({"MALE", "MASCULINO", "H"})] = 1.0
    out[s.str.startswith("F") | s.isin({"FEMALE", "FEMININO", "W"})] = 0.0
    num = pd.to_numeric(raw, errors="coerce")
    coding = "text"
    unresolved = out.isna() & num.notna()
    if unresolved.any():
        vals = set(np.unique(num[unresolved].to_numpy()))
        if vals <= {1.0, 2.0}:
            out[unresolved] = np.where(num[unresolved] == 1.0, 1.0, 0.0); coding = "numeric_1M_2F"
        elif vals <= {0.0, 1.0}:
            out[unresolved] = num[unresolved].astype(float); coding = "numeric_1M_0F"
    report = {"n_total": int(len(raw)), "n_mapped": int(out.notna().sum()),
              "n_unmapped": int(out.isna().sum()), "coding_detected": coding,
              "distinct_input_values": sorted(map(str, pd.unique(raw.astype(str))))[:12]}
    if report["n_mapped"] == 0:
        warnings.warn("B9: nenhuma linha de sexo pode ser mapeada (valores: %s)."
                      % report["distinct_input_values"], RuntimeWarning)
    return out, report


def audit_expected_cohort(observed_ids: Sequence[str],
                          ccfg: CohortConfig) -> Dict[str, Any]:
    """B11: confere a coorte observada contra a lista canonica de IDs esperados.

    ``expected_n`` sozinho responde QUANTOS faltam; a lista de IDs responde QUAIS.
    Sem identidade declarada, um arquivo ausente e um arquivo nunca esperado sao
    indistinguiveis, e as contagens tipo CONSORT ficam sem denominador.
    """
    obs = sorted(set(map(str, observed_ids)))
    if not ccfg.expected_subject_ids:
        return {"expected_declared": False, "n_observed": len(obs)}
    exp = sorted(set(ccfg.expected_subject_ids))
    return {"expected_declared": True, "n_expected": len(exp), "n_observed": len(obs),
            "missing_ids": sorted(set(exp) - set(obs)),
            "unexpected_ids": sorted(set(obs) - set(exp)),
            "coverage": float(len(set(exp) & set(obs)) / len(exp))}


def load_modma_metadata(path: str | Path) -> pd.DataFrame:
    """Le a planilha, canoniza IDs, mapeia MDD/HC e extrai covariaveis demograficas."""
    p = Path(path)
    raw = pd.read_csv(p) if p.suffix.lower() == ".csv" else pd.read_excel(p)
    if raw.empty:
        raise ValueError("Planilha de metadados vazia.")
    raw = raw.loc[:, [c for c in raw.columns if not str(c).lower().startswith("unnamed")]]
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]
    sidcols = [c for c in raw.columns
               if c in {"subject_id", "subject", "id", "participant_id"} or c.startswith("subject")]
    groupcols = [c for c in raw.columns
                 if c in {"type", "group", "class", "diagnosis", "diagnostic_group"}]
    if len(sidcols) != 1 or len(groupcols) != 1:
        raise ValueError("Esperada exatamente uma coluna de ID e uma de grupo.")
    out = raw.rename(columns={sidcols[0]: "subject_id", groupcols[0]: "label_raw"}).copy()
    out["subject_id"] = canonical_subject_id(out["subject_id"])
    if out["subject_id"].duplicated().any():
        raise ValueError("IDs duplicados nos metadados.")
    lab = out["label_raw"].astype(str).str.upper().str.strip()
    out["label"] = np.where(lab.str.startswith("MDD"), 1, np.where(lab.str.startswith("HC"), 0, -1))
    if (out["label"] < 0).any():
        raise ValueError("Rotulo nao mapeavel para MDD/HC.")
    for c in list(out.columns):
        if c.startswith("education"):
            out = out.rename(columns={c: "education_years"})
        elif c.startswith("gender") or c.startswith("sex"):
            out["sex"], sex_report = parse_sex_column(out[c])
    if "sex" in out.columns and not out["sex"].notna().any():
        warnings.warn("B9: coluna de sexo presente mas totalmente nao mapeada; os modelos "
                      "ajustados por confundidores ficariam silenciosamente vazios.",
                      RuntimeWarning)
    out["id_prefix"] = out["subject_id"].str[:4]
    keep = [c for c in ("subject_id", "label", "label_raw", "age", "sex",
                        "education_years", "id_prefix") if c in out.columns]
    return out[keep].drop_duplicates("subject_id").reset_index(drop=True)

# ==================================================================================
# BLOCO 3/11 - ESPECTRO, FILTRAGEM E DIAGNOSTICOS
# ==================================================================================

def psd_continuous(x: np.ndarray, fs: float, nperseg_s: float = 4.0):
    """PSD de Welch do sinal continuo (janelas de nperseg_s segundos)."""
    nper = min(int(round(nperseg_s * fs)), x.shape[-1])
    return _welch(x, fs=fs, nperseg=nper, noverlap=nper // 2, detrend="constant", axis=-1)


def band_share(f: np.ndarray, psd: np.ndarray, lo: float, hi: float,
               ref: Tuple[float, float]) -> np.ndarray:
    """Fracao adimensional da potencia de ``ref`` contida em [lo, hi)."""
    num = psd[..., (f >= lo) & (f < hi)].sum(axis=-1)
    den = psd[..., (f >= ref[0]) & (f < ref[1])].sum(axis=-1) + 1e-30
    return num / den


def line_noise_ratio(x: np.ndarray, fs: float, hz: float = 50.0, bw: float = 1.0,
                     ref: Tuple[float, float] = (1.0, 45.0)) -> float:
    """Evidencia MEDIDA de interferencia de rede eletrica; base da decisao de notch."""
    f, p = psd_continuous(x, fs)
    return float(np.median(np.atleast_1d(band_share(f, p, hz - bw, hz + bw, ref))))


def verify_sampling_rate(x: np.ndarray, acq: AcquisitionConfig) -> Dict[str, Any]:
    """B2: verificacao EMPIRICA da taxa de amostragem pelo pico de rede eletrica.

    ``fs`` nao e lida do arquivo: os TXT do MODMA sao matrizes de contas sem
    cabecalho. O descritor que a documenta especifica repouso de ~5 min, mas os
    registros tem ~1200 s - a mesma fonte e contraditada pelos dados. Portanto
    ``fs`` e SUPOSICAO, e esta funcao busca corroboracao independente.

    Metodo: sob a suposicao ``fs``, a interferencia de rede aparece em 50 Hz (China,
    MODMA) ou 60 Hz. Se o pico observado no espectro nao coincidir com nenhuma das
    frequencias candidatas dentro de ``tol``, a suposicao de ``fs`` fica SEM
    corroboracao e todo o eixo de frequencia escala por um fator desconhecido.

    LIMITACAO: a ausencia de pico de rede NAO refuta ``fs`` (o aparelho pode ter
    filtro de rede interno); apenas deixa a suposicao sem verificacao. A presenca
    do pico na posicao esperada e evidencia consistente, nao prova.
    """
    y = np.asarray(x, float)
    y = y - y.mean(axis=-1, keepdims=True)
    f, p = psd_continuous(y, acq.fs, nperseg_s=8.0)
    p = np.median(np.atleast_2d(p), axis=0)
    band = (f >= 35.0) & (f <= min(70.0, acq.fs / 2.0 - 1.0))
    if band.sum() < 8:
        return {"fs_assumed_hz": float(acq.fs), "fs_verification": "insufficient_band",
                "fs_supported": False}
    fb, pb = f[band], np.log10(p[band] + 1e-30)
    i = int(np.argmax(pb))
    peak_hz = float(fb[i])
    prominence = float(pb[i] - np.median(pb))
    cand = [c for c in acq.fs_verification_line_hz
            if abs(peak_hz - c) <= acq.fs_verification_tol_hz]
    has_peak = prominence >= acq.fs_verification_min_prominence_log10
    if has_peak and cand:
        status, supported = "line_peak_at_expected_frequency", True
    elif has_peak:
        status, supported = "line_peak_at_unexpected_frequency", False
    else:
        status, supported = "no_line_peak_detected_inconclusive", False
    implied = [float(acq.fs * c / peak_hz) for c in acq.fs_verification_line_hz] if has_peak else []
    return {"fs_assumed_hz": float(acq.fs), "fs_is_assumption": bool(acq.fs_is_assumption),
            "line_peak_hz_under_assumed_fs": peak_hz,
            "line_peak_prominence_log10": prominence,
            "fs_verification": status, "fs_supported": bool(supported),
            "implied_fs_if_peak_is_line_hz": implied}


def preprocess_continuous(x_counts: np.ndarray, acq: AcquisitionConfig,
                          pre: PreprocConfig) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Notch condicional + passa-banda de fase zero + corte de borda, no sinal continuo."""
    y = np.asarray(x_counts, float)
    nyq = acq.fs / 2.0
    pre_ratio = line_noise_ratio(y - y.mean(axis=-1, keepdims=True), acq.fs, pre.notch_hz or 50.0)
    if pre.notch_mode == "always":
        do_notch = pre.notch_hz is not None
    elif pre.notch_mode == "never":
        do_notch = False
    else:
        do_notch = (pre.notch_hz is not None) and (pre_ratio > pre.notch_line_ratio_threshold)
    if do_notch:
        bn, an = iirnotch(pre.notch_hz, Q=pre.notch_q, fs=acq.fs)
        y = filtfilt(bn, an, y, axis=-1)
    lo, hi = pre.bandpass_hz
    b, a = butter(pre.filter_order, [lo / nyq, hi / nyq], btype="band")
    y = filtfilt(b, a, y, axis=-1)
    k = int(round(pre.edge_trim_seconds * acq.fs))
    if k > 0:
        if y.shape[1] <= 2 * k:
            raise ValueError("Sinal curto demais para o corte de borda.")
        y = y[:, k:-k]
    y = np.ascontiguousarray(y - y.mean(axis=1, keepdims=True))
    diag = {"notch_mode": pre.notch_mode, "notch_applied": bool(do_notch),
            "line_ratio_pre": float(pre_ratio),
            "line_ratio_post": float(line_noise_ratio(y, acq.fs, pre.notch_hz or 50.0)),
            "edge_trim_s": float(pre.edge_trim_seconds),
            "n_samples_after_trim": int(y.shape[1])}
    return y, diag

# ==================================================================================
# BLOCO 4/11 - A5: PROTOCOLO TEMPORAL DERIVADO DOS DADOS
# ==================================================================================

_PROFILE_BANDS = (("delta", 1.0, 4.0), ("theta", 4.0, 8.0), ("alpha", 8.0, 13.0),
                  ("beta", 13.0, 30.0), ("gamma_low", 30.0, 40.0))


def _spectral_profile(seg: np.ndarray, fs: float) -> np.ndarray:
    """Vetor de participacao por banda (soma 1), medio entre canais: assinatura de estado."""
    f, p = psd_continuous(seg, fs, nperseg_s=min(4.0, seg.shape[-1] / fs))
    v = np.array([np.mean(band_share(f, p, lo, hi, (1.0, 40.0))) for _, lo, hi in _PROFILE_BANDS])
    s = v.sum()
    return v / s if s > 0 else np.full(v.shape, 1.0 / v.size)


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Divergencia de Jensen-Shannon (base 2) entre dois perfis espectrais em [0,1]."""
    p = np.clip(np.asarray(p, float), 1e-12, None); p = p / p.sum()
    q = np.clip(np.asarray(q, float), 1e-12, None); q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: float(np.sum(a * np.log2(a / b)))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def block_profile(filtered: np.ndarray, fs: float, block_seconds: float) -> pd.DataFrame:
    """A5: descreve o registro em blocos consecutivos (RMS, ocular, bandas).

    Cada linha e um bloco; as colunas descrevem o ESTADO do sinal naquele bloco.
    E a materia-prima para decidir quando o registro se estabiliza.
    """
    n = int(round(block_seconds * fs))
    nb = filtered.shape[1] // n
    rows = []
    for i in range(nb):
        seg = filtered[:, i * n:(i + 1) * n]
        f, p = psd_continuous(seg, fs, nperseg_s=min(4.0, block_seconds))
        prof = _spectral_profile(seg, fs)
        rows.append({"block_index": i, "t_start_s": i * block_seconds,
                     "rms_counts": float(np.median(np.sqrt((seg ** 2).mean(axis=1)))),
                     "ocular_index": float(np.median(band_share(f, p, 0.5, 3.0, (1.0, 40.0)))),
                     "muscle_ratio": float(np.median(band_share(f, p, 20.0, 40.0, (1.0, 40.0)))),
                     **{"prof_%s" % nm: float(prof[j])
                        for j, (nm, _, _) in enumerate(_PROFILE_BANDS)}})
    return pd.DataFrame(rows)


def _profile_matrix(df: pd.DataFrame) -> np.ndarray:
    return df[["prof_%s" % nm for nm, _, _ in _PROFILE_BANDS]].to_numpy(float)


def subject_settling_time(prof: pd.DataFrame, tcfg: TemporalProtocolConfig,
                          cohort_ref_profile: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """B3: instante a partir do qual o registro do sujeito fica ESTAVEL.

    PROBLEMA CORRIGIDO. Na v23 a referencia de estabilidade era a mediana da
    SEGUNDA METADE DO PROPRIO REGISTRO. Sob deriva monotona - e 31 dos 55
    registros tem Spearman(RMS, tempo) < -0,5 - a referencia tambem deriva: o
    metodo compara o sinal com um alvo movel e nao pode, por construcao, detectar
    deriva que atravessa o registro inteiro. Isso tende a declarar estabilidade
    cedo demais.

    CORRECAO. A referencia ESPECTRAL passa a ser EXTERNA ao sujeito: mediana da
    coorte nos blocos terminais (``cohort_ref_profile``), calculada uma unica vez
    e cega ao rotulo. O perfil espectral e adimensional (fracoes de banda), logo
    comparavel entre sujeitos.

    LIMITACAO DECLARADA. A referencia de RMS permanece intra-sujeito (mediana dos
    ``n_terminal_blocks`` finais), porque a amplitude em contas depende do ganho
    de contato e NAO e comparavel entre sujeitos. Deriva que abrange todo o
    registro nao e corrigida: e SINALIZADA por ``drift_spans_record``, derivado da
    tendencia de RMS nos blocos terminais. Quando esse sinalizador esta ativo, o
    tempo de acomodacao do sujeito deve ser lido como limite inferior.
    """
    n = len(prof)
    if n < 4:
        return {"settling_time_s": float("nan"), "n_blocks": int(n),
                "settling_status": "insufficient_blocks", "drift_spans_record": False,
                "ref_mode": tcfg.settling_ref_mode}
    k = max(2, min(int(tcfg.n_terminal_blocks), n // 2))
    tail = prof.iloc[-k:]
    rms_ref = float(np.median(tail["rms_counts"])) + 1e-30
    M = _profile_matrix(prof)
    if tcfg.settling_ref_mode == "cohort" and cohort_ref_profile is not None:
        prof_ref = np.asarray(cohort_ref_profile, float)
        ref_used = "cohort_terminal_blocks"
    else:
        prof_ref = np.median(_profile_matrix(tail), axis=0)
        ref_used = "self_terminal_blocks"
    dev_log2 = np.abs(np.log2((prof["rms_counts"].to_numpy(float) + 1e-30) / rms_ref))
    dev_js = np.array([jensen_shannon(M[i], prof_ref) for i in range(n)])
    stable = (dev_log2 <= tcfg.settle_tol_log2) & (dev_js <= tcfg.settle_tol_js)
    idx = len(stable)
    for i in range(len(stable) - 1, -1, -1):
        if stable[i]:
            idx = i
        else:
            break
    t = float(prof["t_start_s"].iloc[idx]) if idx < len(prof) else float("nan")
    tail_rho = float(stats.spearmanr(tail["t_start_s"], tail["rms_counts"]).statistic) \
        if k >= 3 else float("nan")
    spans = bool(np.isfinite(tail_rho) and abs(tail_rho) >= tcfg.tail_trend_warn)
    return {"settling_time_s": t, "n_blocks": int(n),
            "frac_blocks_stable": float(np.mean(stable)),
            "rms_trend_spearman": float(stats.spearmanr(prof["t_start_s"],
                                                        prof["rms_counts"]).statistic),
            "tail_rms_trend_spearman": tail_rho,
            "drift_spans_record": spans,
            "ref_mode": tcfg.settling_ref_mode, "reference_used": ref_used,
            "n_terminal_blocks": int(k),
            "dev_log2_first_block": float(dev_log2[0]), "dev_js_first_block": float(dev_js[0]),
            "settling_status": ("ok" if np.isfinite(t) else "never_stable")}


def cohort_terminal_profile(block_frames: Sequence[pd.DataFrame],
                            tcfg: TemporalProtocolConfig) -> Optional[np.ndarray]:
    """B3: perfil espectral de referencia EXTERNO, mediana da coorte nos blocos finais.

    Calculado uma unica vez, sobre todos os sujeitos, sem qualquer acesso ao rotulo.
    E etapa de COORTE (nao pre-processamento individual) e deve ser descrita como tal.
    """
    rows = []
    for prof in block_frames:
        if len(prof) < 4:
            continue
        k = max(2, min(int(tcfg.n_terminal_blocks), len(prof) // 2))
        rows.append(np.median(_profile_matrix(prof.iloc[-k:]), axis=0))
    if not rows:
        return None
    ref = np.median(np.vstack(rows), axis=0)
    total = ref.sum()
    return ref / total if total > 0 else ref


def derive_temporal_protocol(settle_df: pd.DataFrame, durations: Dict[str, float],
                             tcfg: TemporalProtocolConfig) -> Dict[str, Any]:
    """A5/B4: converte os tempos de acomodacao individuais em UMA regra de coorte.

    Passos: (i) elegibilidade por duracao minima pre-declarada; (ii) inicio comum =
    quantil dos tempos de acomodacao, na grade de blocos, limitado por
    [min_skip_seconds, max_skip_seconds]; (iii) duracao comum = menor sobra entre os
    elegiveis, limitada por target/max; (iv) janelas de sensibilidade nao sobrepostas.
    Nenhum rotulo participa de qualquer etapa.

    B4 CORRIGIDO. Na v23, se o quantil dos tempos de acomodacao excedesse
    ``max_skip_seconds``, o inicio era truncado EM SILENCIO e a janela podia voltar
    a cair no transitorio - exatamente o problema que a derivacao existe para
    evitar. Agora o truncamento e explicito (``skip_capped``, ``window_capped``),
    emite ``warnings.warn`` e marca ``protocol_valid=False``, para que o resultado
    nao seja lido como se a regra tivesse sido respeitada.
    """
    if tcfg.mode == "fixed":
        skip, win = tcfg.fixed_skip_seconds, tcfg.fixed_window_seconds
        eligible = [s for s, d in durations.items() if d >= skip + win]
        return {"rule": "fixed", "skip_seconds": skip, "window_seconds": win,
                "eligible_subjects": sorted(eligible),
                "excluded_short": sorted(set(durations) - set(eligible)),
                "sensitivity_windows": [(skip, win)], "settling_quantile_used": None,
                "skip_capped": False, "window_capped": False, "protocol_valid": True,
                "n_drift_spans_record": 0}
    eligible = [s for s, d in durations.items() if d >= tcfg.min_duration_seconds]
    excluded = sorted(set(durations) - set(eligible))
    if not eligible:
        raise ValueError("Nenhum sujeito atinge min_duration_seconds=%.0f s."
                         % tcfg.min_duration_seconds)
    sdf = settle_df.set_index("subject_id").reindex(eligible)
    st = pd.to_numeric(sdf["settling_time_s"], errors="coerce").to_numpy(float)
    st = st[np.isfinite(st)]
    n_spans = int(pd.Series(sdf.get("drift_spans_record", pd.Series(dtype=bool))).fillna(False).sum())
    q = float(np.quantile(st, tcfg.settle_quantile)) if st.size else tcfg.min_skip_seconds
    skip_raw = float(np.ceil(q / tcfg.block_seconds) * tcfg.block_seconds)
    skip = float(min(max(skip_raw, tcfg.min_skip_seconds), tcfg.max_skip_seconds))
    skip_capped = bool(skip_raw > tcfg.max_skip_seconds)
    if skip_capped:
        warnings.warn("B4: inicio derivado (%.0f s) excede max_skip_seconds (%.0f s) e foi "
                      "truncado; a janela pode conter transitorio de acomodacao."
                      % (skip_raw, tcfg.max_skip_seconds), RuntimeWarning)
    avail = min(durations[s] for s in eligible) - skip
    win_raw = min(avail, tcfg.target_window_seconds)
    win = float(np.floor(min(win_raw, tcfg.max_window_seconds)))
    window_capped = bool(avail < tcfg.target_window_seconds)
    if window_capped:
        warnings.warn("B4: duracao comum limitada pela sobra disponivel (%.0f s) e nao pelo "
                      "alvo (%.0f s)." % (avail, tcfg.target_window_seconds), RuntimeWarning)
    if win <= 0:
        raise ValueError("Janela comum nao positiva; revise min_duration_seconds.")
    k = max(1, int(np.floor(avail / win)))
    k = min(k, tcfg.n_sensitivity_windows)
    if k > 1:
        span = avail - win
        starts = [float(np.floor(skip + span * i / (k - 1))) for i in range(k)]
        for i in range(1, k):  # garante nao sobreposicao
            starts[i] = max(starts[i], starts[i - 1] + win)
        starts = [s for s in starts if s + win <= skip + avail + 1e-9]
    else:
        starts = [skip]
    return {"rule": "derived_from_data", "skip_seconds": skip, "window_seconds": win,
            "skip_seconds_before_cap": skip_raw,
            "skip_capped": skip_capped, "window_capped": window_capped,
            "protocol_valid": bool(not skip_capped),
            "settling_quantile_used": tcfg.settle_quantile,
            "settling_time_quantile_s": q,
            "cohort_settling_median_s": float(np.median(st)) if st.size else None,
            "cohort_settling_max_s": float(np.max(st)) if st.size else None,
            "available_after_skip_s": float(avail),
            "n_drift_spans_record": n_spans,
            "settling_reference_mode": tcfg.settling_ref_mode,
            "eligible_subjects": sorted(eligible), "excluded_short": excluded,
            "sensitivity_windows": [(float(s), float(win)) for s in starts]}


def extract_window(filtered: np.ndarray, fs: float, start_s: float,
                   window_s: float) -> np.ndarray:
    """Recorta a janela [start_s, start_s + window_s) do sinal ja filtrado."""
    a, b = int(round(start_s * fs)), int(round((start_s + window_s) * fs))
    if filtered.shape[1] < b:
        raise ValueError("Sinal insuficiente: %d < %d amostras." % (filtered.shape[1], b))
    return np.ascontiguousarray(filtered[:, a:b])


def within_window_drift(window: np.ndarray, fs: float,
                        block_seconds: float) -> Dict[str, float]:
    """A5: mede se o ESTADO muda dentro da janela escolhida (deriva residual)."""
    prof = block_profile(window, fs, min(block_seconds, window.shape[1] / fs / 3.0))
    if len(prof) < 2:
        return {"drift_log2_rms": float("nan"), "drift_js": float("nan"),
                "drift_rms_spearman": float("nan")}
    M = _profile_matrix(prof)
    med = np.median(M, axis=0)
    r = prof["rms_counts"].to_numpy(float)
    return {"drift_log2_rms": float(np.log2((r.max() + 1e-30) / (r.min() + 1e-30))),
            "drift_js": float(max(jensen_shannon(M[i], med) for i in range(len(prof)))),
            "drift_rms_spearman": float(stats.spearmanr(prof["t_start_s"], r).statistic)}

# ==================================================================================
# BLOCO 5/11 - EPOCAS E CONTROLE DE QUALIDADE
# ==================================================================================

def epoch_signal(x: np.ndarray, fs: float, epoch_seconds: float,
                 overlap: float = 0.0) -> np.ndarray:
    """Segmenta a janela em epocas -> (n_epochs, n_channels, n_times)."""
    n = int(round(epoch_seconds * fs))
    step = max(1, int(round(n * (1.0 - overlap))))
    starts = list(range(0, x.shape[1] - n + 1, step))
    if not starts:
        raise ValueError("Nenhuma epoca gerada.")
    return np.stack([x[:, s:s + n] for s in starts], axis=0)


def compute_psd(epochs: np.ndarray, fs: float, fcfg: FeatureConfig):
    """PSD de Welch por epoca e canal."""
    nper = int(round(fcfg.welch_nperseg_seconds * fs))
    return _welch(epochs, fs=fs, nperseg=nper,
                  noverlap=int(round(nper * fcfg.welch_overlap)),
                  detrend="constant", axis=-1)


def epoch_rejection_mask(epochs: np.ndarray, f: np.ndarray, psd: np.ndarray,
                         qc: QCConfig, fcfg: FeatureConfig) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Rejeicao de epoca relativa ao proprio sujeito, adimensional e cega ao rotulo."""
    finite = np.isfinite(epochs).all(axis=(1, 2))
    p2p = np.ptp(epochs, axis=2)
    jump = np.abs(np.diff(epochs, axis=2)).max(axis=2)
    nonflat = p2p.min(axis=1) > 0
    r_p2p = p2p / (np.median(p2p, axis=0) + 1e-30)
    r_jump = jump / (np.median(jump, axis=0) + 1e-30)
    lf = band_share(f, psd, fcfg.ocular_band[0], fcfg.ocular_band[1], fcfg.total_band)
    r_lf = lf / (np.median(lf, axis=0) + 1e-30)
    ok_hi = (r_p2p <= qc.p2p_ratio_max).all(axis=1)
    ok_lo = (r_p2p >= qc.p2p_ratio_min).all(axis=1)
    ok_jump = (r_jump <= qc.jump_ratio_max).all(axis=1)
    ok_lf = (r_lf <= qc.lf_share_ratio_max).all(axis=1)
    good = finite & nonflat & ok_hi & ok_lo & ok_jump & ok_lf
    diag = {"n_epochs": int(epochs.shape[0]), "n_good": int(good.sum()),
            "n_rejected": int((~good).sum()),
            "frac_rejected": float(1.0 - good.mean()) if good.size else float("nan"),
            "n_rej_nonfinite": int((~finite).sum()), "n_rej_flat": int((~nonflat).sum()),
            "n_rej_p2p_high": int((~ok_hi).sum()), "n_rej_p2p_low": int((~ok_lo).sum()),
            "n_rej_jump": int((~ok_jump).sum()), "n_rej_ocular_lf": int((~ok_lf).sum()),
            "rule": "relative_intrasubject_p2p_jump_lowfreq"}
    return good, diag


def nonstationarity_cv(x: np.ndarray, fs: float, block_s: float = 10.0) -> float:
    """Coeficiente de variacao da potencia em blocos: instabilidade do registro."""
    n = int(round(block_s * fs))
    if x.shape[1] < 3 * n:
        return float("nan")
    nb = x.shape[1] // n
    v = np.array([x[:, i * n:(i + 1) * n].var(axis=1).mean() for i in range(nb)])
    m = float(np.mean(v))
    return float(np.std(v) / m) if m > 0 else float("nan")


def subject_qc_metrics(rec: SubjectRecording, filtered: np.ndarray, acq: AcquisitionConfig,
                       fcfg: FeatureConfig, filt_diag: Dict[str, Any]) -> Dict[str, Any]:
    """QC do REGISTRO INTEIRO: tecnica, adimensional e cega ao rotulo."""
    raw = np.asarray(rec.data_counts, float)
    f_raw, p_raw = psd_continuous(raw - raw.mean(axis=1, keepdims=True), acq.fs)
    f_flt, p_flt = psd_continuous(filtered, acq.fs)
    lo24, hi24 = -(2 ** (acq.adc_bits - 1)), 2 ** (acq.adc_bits - 1) - 1
    diff = np.diff(raw, axis=1)
    corr = np.corrcoef(filtered)
    off = corr[np.triu_indices(corr.shape[0], k=1)]
    rms = float(np.median(np.sqrt((filtered ** 2).mean(axis=1))))
    return {"subject_id": rec.subject_id, "source_name": rec.source_name, "task": rec.task,
            "n_samples": int(raw.shape[1]), "duration_s": float(raw.shape[1] / acq.fs),
            "finite": bool(np.isfinite(raw).all()),
            "flat_channels": int(np.sum(np.ptp(raw, axis=1) == 0)),
            "n_wraparound_fixed": int(rec.parse_report.get("n_wraparound_fixed", 0)),
            "saturation_fraction": float((np.sum(raw <= lo24) + np.sum(raw >= hi24)) / raw.size),
            "zero_diff_fraction": float(np.mean(diff == 0)) if diff.size else float("nan"),
            "line_noise_ratio_raw": float(np.median(band_share(f_raw, p_raw, 49.0, 51.0, (1.0, 45.0)))),
            "line_noise_ratio_post": float(filt_diag.get("line_ratio_post", np.nan)),
            "notch_applied": bool(filt_diag.get("notch_applied", False)),
            "ocular_index": float(np.median(band_share(f_flt, p_flt, *fcfg.ocular_band,
                                                       ref=fcfg.total_band))),
            "muscle_ratio": float(np.median(band_share(f_flt, p_flt, *fcfg.muscle_band,
                                                       ref=fcfg.total_band))),
            "nonstationarity_cv": nonstationarity_cv(filtered, acq.fs),
            "min_channel_corr": float(np.min(off)) if off.size else float("nan"),
            "rms_counts": rms, "log_rms_counts": float(np.log10(rms + 1e-30))}


def qc_distribution_report(qc_df: pd.DataFrame, qc: QCConfig) -> pd.DataFrame:
    """B5: confronta CADA limiar de QC com a distribuicao empirica da coorte.

    ERRO CORRIGIDO. Na v23 os limiares (max_ocular_index=0,80, max_line_noise_ratio
    =0,25, max_muscle_ratio=0,60, max_nonstationarity_cv=1,00) nunca foram
    comparados com os valores realmente observados. O numero de exclusoes que eles
    produzem era desconhecido a priori e podia variar de zero a quase toda a coorte
    - e ninguem saberia antes de rodar.

    Esta tabela informa, por metrica: n valido, minimo, quartis, maximo, o limiar
    vigente e QUANTOS sujeitos ele excluiria (contagem e fracao). E o insumo
    obrigatorio para congelar os limiares antes de habilitar o modo research.
    Nao usa rotulo em nenhum momento.
    """
    rows = []
    for col, (attr, direction) in QC_THRESHOLD_MAP.items():
        if col not in qc_df.columns:
            continue
        v = pd.to_numeric(qc_df[col], errors="coerce").to_numpy(float)
        fin = v[np.isfinite(v)]
        thr = float(getattr(qc, attr))
        n_ex = int(np.sum(fin > thr) if direction == "greater" else np.sum(fin < thr))
        rows.append({"metric": col, "threshold_param": attr, "threshold_value": thr,
                     "direction": direction, "n_valid": int(fin.size),
                     "min": float(fin.min()) if fin.size else np.nan,
                     "p05": float(np.percentile(fin, 5)) if fin.size else np.nan,
                     "q25": float(np.percentile(fin, 25)) if fin.size else np.nan,
                     "median": float(np.median(fin)) if fin.size else np.nan,
                     "q75": float(np.percentile(fin, 75)) if fin.size else np.nan,
                     "p95": float(np.percentile(fin, 95)) if fin.size else np.nan,
                     "max": float(fin.max()) if fin.size else np.nan,
                     "n_excluded_by_threshold": n_ex,
                     "frac_excluded_by_threshold": float(n_ex / fin.size) if fin.size else np.nan})
    for m in qc.cohort_outlier_metrics:
        col = "z__" + m
        if col not in qc_df.columns:
            continue
        z = pd.to_numeric(qc_df[col], errors="coerce").to_numpy(float)
        fin = z[np.isfinite(z)]
        n_ex = int(np.sum(np.abs(fin) > qc.cohort_mad_z_max))
        rows.append({"metric": col, "threshold_param": "cohort_mad_z_max",
                     "threshold_value": float(qc.cohort_mad_z_max), "direction": "abs_greater",
                     "n_valid": int(fin.size),
                     "min": float(fin.min()) if fin.size else np.nan,
                     "p05": float(np.percentile(fin, 5)) if fin.size else np.nan,
                     "q25": float(np.percentile(fin, 25)) if fin.size else np.nan,
                     "median": float(np.median(fin)) if fin.size else np.nan,
                     "q75": float(np.percentile(fin, 75)) if fin.size else np.nan,
                     "p95": float(np.percentile(fin, 95)) if fin.size else np.nan,
                     "max": float(fin.max()) if fin.size else np.nan,
                     "n_excluded_by_threshold": n_ex,
                     "frac_excluded_by_threshold": float(n_ex / fin.size) if fin.size else np.nan})
    return pd.DataFrame(rows)


def assert_qc_thresholds_frozen(freeze: QCFreezeConfig, policy: EvidencePolicy) -> None:
    """B5: trava dura - modo research exige limiares de QC congelados e datados."""
    if policy.mode != "research" or not freeze.require_freeze_for_research:
        return
    if not freeze.thresholds_frozen:
        raise RuntimeError(
            "B5: modo research bloqueado. Os limiares de QC ainda nao foram congelados. "
            "Execute em modo integration, inspecione 02_quality/qc_distribution.csv (que "
            "informa quantos sujeitos cada limiar excluiria), fixe os valores e declare "
            "QCFreezeConfig(thresholds_frozen=True, freeze_date=..., freeze_source=...).")


def cohort_qc_decision(qc_df: pd.DataFrame, qc: QCConfig) -> pd.DataFrame:
    """Decisao de exclusao PRE-ESPECIFICADA: limiares absolutos + outlier de coorte (MAD)."""
    out = qc_df.copy()
    reasons: List[List[str]] = [[] for _ in range(len(out))]

    def _abs_rule(col: str, thr: float, greater: bool, tag: str) -> None:
        if col not in out.columns:
            return
        v = pd.to_numeric(out[col], errors="coerce").to_numpy(float)
        bad = (v > thr) if greater else (v < thr)
        for i in np.where(np.isfinite(v) & bad)[0]:
            reasons[i].append("%s=%.4g" % (tag, v[i]))

    _abs_rule("line_noise_ratio_post", qc.max_line_noise_ratio, True, "line_noise")
    _abs_rule("ocular_index", qc.max_ocular_index, True, "ocular")
    _abs_rule("muscle_ratio", qc.max_muscle_ratio, True, "muscle")
    _abs_rule("zero_diff_fraction", qc.max_zero_diff_fraction, True, "frozen_samples")
    _abs_rule("saturation_fraction", qc.max_saturation_fraction, True, "saturation")
    _abs_rule("nonstationarity_cv", qc.max_nonstationarity_cv, True, "nonstationarity")
    _abs_rule("min_channel_corr", qc.min_channel_corr, False, "channel_corr")
    if "flat_channels" in out.columns:
        for i in np.where(pd.to_numeric(out["flat_channels"], errors="coerce").fillna(0) > 0)[0]:
            reasons[i].append("flat_channel")
    for m in qc.cohort_outlier_metrics:
        if m not in out.columns:
            continue
        v = pd.to_numeric(out[m], errors="coerce").to_numpy(float)
        med = np.nanmedian(v)
        mad = np.nanmedian(np.abs(v - med))
        scale = 1.4826 * mad if mad > 0 else np.nanstd(v)
        z = (v - med) / (scale + 1e-30)
        out["z__" + m] = z
        for i in np.where(np.isfinite(z) & (np.abs(z) > qc.cohort_mad_z_max))[0]:
            reasons[i].append("cohort_outlier_%s(z=%.1f)" % (m, z[i]))
    out["qc_reasons"] = ["; ".join(r) for r in reasons]
    out["qc_pass"] = [len(r) == 0 for r in reasons]
    return out

# ==================================================================================
# BLOCO 6/11 - FEATURES
# ==================================================================================

def _mask(f: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return (f >= lo) & (f < hi)


def relative_band_powers(f, psd, fcfg: FeatureConfig) -> Dict[str, np.ndarray]:
    """Potencia RELATIVA por banda: adimensional, imune ao ganho de contato."""
    tot = psd[..., _mask(f, *fcfg.total_band)].sum(axis=-1) + 1e-30
    return {n: psd[..., _mask(f, lo, hi)].sum(axis=-1) / tot for n, lo, hi in fcfg.bands}


def spectral_entropy(f, psd, fcfg: FeatureConfig) -> np.ndarray:
    """Entropia de Shannon do espectro normalizado, em [0,1]: 1 = espectro plano."""
    sel = psd[..., _mask(f, *fcfg.total_band)]
    p = sel / (sel.sum(axis=-1, keepdims=True) + 1e-30)
    return -(p * np.log(p + 1e-30)).sum(axis=-1) / np.log(p.shape[-1])


def robust_aperiodic_fit(f, psd, fcfg: FeatureConfig):
    """Ajuste log-log robusto do fundo 1/f, removendo picos iterativamente.

    Devolve (inclinacao, R2, espectro achatado). A inclinacao e o expoente
    aperiodico com sinal invertido: quanto maior, mais "inclinado" o espectro.
    """
    m = _mask(f, *fcfg.aperiodic_fit_band) & (f > 0)
    lf = np.log10(f[m])
    ly = np.log10(psd[..., m] + 1e-30)
    shape = ly.shape[:-1]
    lyf = ly.reshape(-1, ly.shape[-1])
    slopes = np.empty(lyf.shape[0]); r2s = np.empty(lyf.shape[0])
    flat = np.empty_like(lyf)
    X = np.vstack([lf, np.ones_like(lf)]).T
    for i in range(lyf.shape[0]):
        y = lyf[i]
        keep = np.ones_like(y, bool)
        beta = np.array([0.0, 0.0])
        for _ in range(fcfg.aperiodic_max_iter):
            beta, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            res = y - X @ beta
            sd = np.std(res[keep]) + 1e-30
            new = res <= fcfg.aperiodic_peak_sd * sd
            if new.sum() < 5 or np.array_equal(new, keep):
                keep = new if new.sum() >= 5 else keep
                break
            keep = new
        pred = X @ beta
        ss = np.sum((y[keep] - pred[keep]) ** 2)
        st = np.sum((y[keep] - y[keep].mean()) ** 2) + 1e-30
        slopes[i] = beta[0]; r2s[i] = 1.0 - ss / st
        flat[i] = y - pred
    return (slopes.reshape(shape), r2s.reshape(shape),
            flat.reshape(ly.shape), f[m])


def alpha_peak_from_flat(flat: np.ndarray, f_fit: np.ndarray,
                         fcfg: FeatureConfig) -> np.ndarray:
    """Pico alfa detectado no espectro ACHATADO; sem pico identificavel -> NaN."""
    m = _mask(f_fit, *fcfg.alpha_peak_band)
    if not m.any():
        return np.full(flat.shape[:-1], np.nan)
    sub = flat[..., m]
    idx = np.argmax(sub, axis=-1)
    height = np.take_along_axis(sub, idx[..., None], axis=-1)[..., 0]
    peak = f_fit[m][idx]
    return np.where(height >= fcfg.alpha_min_peak_log10, peak, np.nan)


def hjorth_params(epochs: np.ndarray):
    """Mobilidade e complexidade de Hjorth: razoes de desvios, invariantes a escala."""
    d1 = np.diff(epochs, axis=-1)
    d2 = np.diff(d1, axis=-1)
    s0 = epochs.std(axis=-1) + 1e-30
    s1 = d1.std(axis=-1) + 1e-30
    s2 = d2.std(axis=-1) + 1e-30
    mob = s1 / s0
    return mob, (s2 / s1) / mob


def permutation_entropy(x: np.ndarray, order: int = 3, delay: int = 1) -> float:
    """Entropia de permutacao normalizada: previsibilidade dos padroes ordinais."""
    n = len(x) - (order - 1) * delay
    if n <= 1:
        return float("nan")
    idx = np.arange(order) * delay
    emb = np.stack([x[i + idx] for i in range(n)])
    perm = np.argsort(emb, axis=1)
    _, counts = np.unique(perm, axis=0, return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum() / np.log(math.factorial(order)))


def higuchi_fd(x: np.ndarray, kmax: int = 10) -> float:
    """Dimensao fractal de Higuchi: rugosidade da serie no tempo."""
    n = len(x)
    lk = []
    ks = []
    for k in range(1, kmax + 1):
        lm = []
        for m in range(k):
            idx = np.arange(m, n, k)
            if idx.size < 2:
                continue
            lm.append(np.abs(np.diff(x[idx])).sum() * (n - 1) / (idx.size - 1) / k)
        if lm:
            lk.append(np.log(np.mean(lm) + 1e-30)); ks.append(np.log(1.0 / k))
    if len(ks) < 2:
        return float("nan")
    return float(np.polyfit(ks, lk, 1)[0])


def c0_complexity(x: np.ndarray) -> float:
    """Complexidade C0: fracao de energia nao explicada pelas componentes dominantes."""
    X = np.fft.fft(x)
    p = np.abs(X) ** 2
    mean_p = p.mean()
    Y = np.where(p > mean_p, X, 0.0)
    y = np.fft.ifft(Y)
    num = float(np.sum(np.abs(x - y) ** 2))
    den = float(np.sum(np.abs(x) ** 2)) + 1e-30
    return num / den


def epoch_feature_arrays(epochs: np.ndarray, fs: float, fcfg: FeatureConfig,
                         channels: Sequence[str]) -> Dict[str, np.ndarray]:
    """Calcula todas as features no nivel EPOCA x CANAL."""
    f, psd = compute_psd(epochs, fs, fcfg)
    out: Dict[str, np.ndarray] = {}
    for name, val in relative_band_powers(f, psd, fcfg).items():
        out["rel_" + name] = val
    out["spec_entropy"] = spectral_entropy(f, psd, fcfg)
    slope, r2, flat, f_fit = robust_aperiodic_fit(f, psd, fcfg)
    out["aperiodic_slope"] = slope
    out["aperiodic_r2"] = r2
    out["alpha_peak_hz"] = alpha_peak_from_flat(flat, f_fit, fcfg)
    out["hjorth_mobility"], out["hjorth_complexity"] = hjorth_params(epochs)
    ne, nc, _ = epochs.shape
    pe = np.empty((ne, nc)); hf = np.empty((ne, nc)); c0 = np.empty((ne, nc))
    for i in range(ne):
        for j in range(nc):
            seg = epochs[i, j]
            pe[i, j] = permutation_entropy(seg, fcfg.perm_entropy_order, fcfg.perm_entropy_delay)
            hf[i, j] = higuchi_fd(seg, fcfg.higuchi_kmax)
            c0[i, j] = c0_complexity(seg)
    out["perm_entropy"], out["higuchi_fd"], out["c0_complexity"] = pe, hf, c0
    a = psd[..., _mask(f, 8.0, 13.0)].sum(axis=-1)
    out["faa"] = np.log(a[:, 2] + 1e-30) - np.log(a[:, 0] + 1e-30)
    return out


def aggregate_features(arrays: Dict[str, np.ndarray], channels: Sequence[str],
                       fcfg: FeatureConfig) -> Dict[str, float]:
    """Agrega epoca->participante e canal->participante. p__ = confirmatorio, e__ = exploratorio."""
    agg = np.nanmedian if fcfg.aggregate == "median" else np.nanmean
    row: Dict[str, float] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in PRIMARY_FEATURES:
            v = arrays[k]
            row["p__" + k] = float(agg(np.nanmedian(v, axis=1)))
            for j, ch in enumerate(channels):
                row["e__%s_%s" % (k, ch)] = float(agg(v[:, j]))
            q = np.nanpercentile(np.nanmedian(v, axis=1), [25, 75])
            row["e__iqr_" + k] = float(q[1] - q[0])
        row["e__faa"] = float(agg(arrays["faa"]))
        row["e__aperiodic_r2"] = float(agg(np.nanmedian(arrays["aperiodic_r2"], axis=1)))
        row["e__alpha_detect_rate"] = float(np.mean(np.isfinite(arrays["alpha_peak_hz"])))
        for k in ("perm_entropy", "higuchi_fd", "c0_complexity"):
            v = arrays[k]
            row["e__" + k] = float(agg(np.nanmedian(v, axis=1)))
            for j, ch in enumerate(channels):
                row["e__%s_%s" % (k, ch)] = float(agg(v[:, j]))
    return row


def primary_columns(df: pd.DataFrame) -> List[str]:
    return sorted(c for c in df.columns if c.startswith("p__"))


def exploratory_columns(df: pd.DataFrame) -> List[str]:
    return sorted(c for c in df.columns if c.startswith("e__"))


def assert_primary_only(cols: Sequence[str]) -> None:
    """Barreira dura confirmatorio/exploratorio."""
    bad = [c for c in cols if not c.startswith("p__")]
    if bad:
        raise RuntimeError("Colunas nao primarias no modelo confirmatorio: %s" % bad[:8])


def assert_no_circular_features(cols: Sequence[str]) -> None:
    """Escalas clinicas e o proprio rotulo nunca podem ser preditores."""
    bad = [c for c in cols if any(t in c.lower() for t in CLINICAL_SCALE_TOKENS)
           or c.lower() in {"label", "y", "type", "label_raw"}]
    if bad:
        raise RuntimeError("Features circulares/clinicas proibidas: %s" % bad[:8])

# ==================================================================================
# BLOCO 7/11 - ESTATISTICA
# ==================================================================================

def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Tamanho de efeito padronizado com correcao para amostras pequenas."""
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    b = np.asarray(b, float); b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    na, nb = a.size, b.size
    sp = math.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp <= 0:
        return float("nan")
    d = (a.mean() - b.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (na + nb) - 9.0)
    return float(d * j)


def bootstrap_ci(func, *arrays, n_boot: int = 2000, alpha: float = 0.05,
                 seed: int = 0) -> Tuple[float, float]:
    """IC percentil por reamostragem de PARTICIPANTES (unidade independente)."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        res = [a[rng.integers(0, len(a), len(a))] for a in arrays]
        v = func(*res)
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 100 * alpha / 2)),
            float(np.percentile(vals, 100 * (1 - alpha / 2))))


def bh_fdr(p: Sequence[float]) -> np.ndarray:
    """Correcao de Benjamini-Hochberg para multiplas comparacoes."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def cramers_v(x: Sequence, y: Sequence) -> float:
    """V de Cramer entre duas variaveis categoricas (0 = independentes, 1 = colineares)."""
    tab = pd.crosstab(pd.Series(x), pd.Series(y)).to_numpy()
    if tab.size == 0 or tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    n = tab.sum()
    return float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))

# ==================================================================================
# BLOCO 8/11 - ANALISE CONFIRMATORIA E MODELAGEM
# ==================================================================================

def confirmatory_group_tests(df: pd.DataFrame, acfg: AnalysisConfig) -> pd.DataFrame:
    """Testes da familia confirmatoria (p__), com efeito, IC e FDR de Benjamini-Hochberg."""
    cols = primary_columns(df)
    assert_primary_only(cols)
    assert_no_circular_features(cols)
    y = df["label"].to_numpy(int)
    rows = []
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        a, b = v[y == 1], v[y == 0]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if a.size < 3 or b.size < 3:
            rows.append({"feature": c, "n_mdd": a.size, "n_hc": b.size, "p_raw": np.nan,
                         "hedges_g": np.nan, "g_ci_low": np.nan, "g_ci_high": np.nan})
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        g = hedges_g(a, b)
        lo, hi = bootstrap_ci(hedges_g, a, b, n_boot=acfg.n_boot_effect,
                              seed=acfg.random_seed)
        rows.append({"feature": c, "n_mdd": int(a.size), "n_hc": int(b.size),
                     "median_mdd": float(np.median(a)), "median_hc": float(np.median(b)),
                     "u_statistic": float(u.statistic), "p_raw": float(u.pvalue),
                     "hedges_g": g, "g_ci_low": lo, "g_ci_high": hi})
    out = pd.DataFrame(rows)
    ok = out["p_raw"].notna()
    out.loc[ok, "p_fdr"] = bh_fdr(out.loc[ok, "p_raw"].to_numpy(float))
    out["significant_fdr"] = out["p_fdr"] < acfg.fdr_alpha
    return out.sort_values("p_raw").reset_index(drop=True)


def confound_adjusted_models(df: pd.DataFrame, acfg: AnalysisConfig) -> pd.DataFrame:
    """Regressao logistica por feature ajustada por idade, sexo e escolaridade."""
    import statsmodels.api as sm  # dependencia opcional
    cols = primary_columns(df)
    covs = [c for c in acfg.covariates if c in df.columns]
    rows = []
    for c in cols:
        sub = df[[c, "label"] + covs].apply(pd.to_numeric, errors="coerce").dropna()
        if len(sub) < 20:
            continue
        X = sm.add_constant(np.column_stack(
            [stats.zscore(sub[c].to_numpy(float))] +
            [stats.zscore(sub[k].to_numpy(float)) for k in covs]))
        try:
            res = sm.Logit(sub["label"].to_numpy(int), X).fit(disp=0)
            rows.append({"feature": c, "n": len(sub), "beta_z": float(res.params[1]),
                         "or_per_sd": float(np.exp(res.params[1])),
                         "p_value": float(res.pvalues[1]),
                         "ci_low": float(np.exp(res.conf_int()[1][0])),
                         "ci_high": float(np.exp(res.conf_int()[1][1])),
                         "covariates": ",".join(covs)})
        except Exception as exc:
            rows.append({"feature": c, "n": len(sub), "error": str(exc)[:120]})
    out = pd.DataFrame(rows)
    if "p_value" in out.columns and out["p_value"].notna().any():
        ok = out["p_value"].notna()
        out.loc[ok, "p_fdr"] = bh_fdr(out.loc[ok, "p_value"].to_numpy(float))
    return out


def _build_pipeline(c_grid, inner_folds, scoring, seed):
    """Pipeline scikit-learn com TODA transformacao dentro do fold (evita vazamento)."""
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("sc", StandardScaler()),
                     ("clf", LogisticRegression(penalty="l2", solver="liblinear",
                                                max_iter=5000))])
    return GridSearchCV(pipe, {"clf__C": list(c_grid)},
                        cv=StratifiedKFold(inner_folds, shuffle=True, random_state=seed),
                        scoring=scoring, refit=True, n_jobs=1)


def nested_cv_scores(X: np.ndarray, y: np.ndarray, acfg: AnalysisConfig,
                     repeats: Optional[int] = None) -> Dict[str, Any]:
    """Validacao cruzada ANINHADA no nivel do participante.

    O fold externo estima desempenho; o fold interno escolhe o hiperparametro.
    Nada e ajustado fora do fold, o que evita otimismo por vazamento de dados.
    """
    if not _SKLEARN_OK:
        raise RuntimeError("scikit-learn indisponivel: %s" % _SKLEARN_ERR)
    reps = acfg.outer_repeats if repeats is None else repeats
    cv = RepeatedStratifiedKFold(n_splits=acfg.outer_folds, n_repeats=reps,
                                 random_state=acfg.random_seed)
    oof = np.full((reps, len(y)), np.nan)
    for i, (tr, te) in enumerate(cv.split(X, y)):
        gs = _build_pipeline(acfg.logistic_c_grid, acfg.inner_folds,
                             acfg.inner_scoring, acfg.random_seed)
        gs.fit(X[tr], y[tr])
        oof[i // acfg.outer_folds, te] = gs.predict_proba(X[te])[:, 1]
    return {"oof_prob_per_repeat": oof, "prob_mean": np.nanmean(oof, axis=0)}


def metrics_from_probs(y: np.ndarray, p: np.ndarray, acfg: AnalysisConfig,
                       oof_per_repeat: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Metricas de classificacao com IC bootstrap por participante.

    B11: o limiar de decisao volta a ser parametro explicito (``decision_threshold``),
    em vez do 0,5 fixo no corpo da funcao.

    B10: quando ``oof_per_repeat`` e fornecido, o IC da AUC passa a reamostrar
    conjuntamente PARTICIPANTES e REPETICOES da validacao cruzada. O IC da v23
    tratava as probabilidades medias como fixas e, por isso, subestimava a
    incerteza real, que inclui a variabilidade das particoes.
    """
    thr = float(acfg.decision_threshold)
    yhat = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yhat, labels=[0, 1]).ravel()
    auc = float(roc_auc_score(y, p))
    lo, hi = bootstrap_ci(lambda idx: float(roc_auc_score(y[idx], p[idx]))
                          if len(np.unique(y[idx])) == 2 else np.nan,
                          np.arange(len(y)), n_boot=acfg.n_boot_metrics,
                          seed=acfg.random_seed)
    out = {"n": int(len(y)), "n_mdd": int(y.sum()), "n_hc": int((1 - y).sum()),
           "decision_threshold": thr,
           "roc_auc": auc, "roc_auc_ci95_subjects_only": [lo, hi],
           "pr_auc": float(average_precision_score(y, p)),
           "balanced_accuracy": float(balanced_accuracy_score(y, yhat)),
           "sensitivity": float(tp / (tp + fn)) if (tp + fn) else float("nan"),
           "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
           "f1": float(f1_score(y, yhat)), "mcc": float(matthews_corrcoef(y, yhat)),
           "brier": float(brier_score_loss(y, p)),
           "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}}
    if oof_per_repeat is not None and np.ndim(oof_per_repeat) == 2:
        M = np.asarray(oof_per_repeat, float)
        n_rep, n_sub = M.shape
        rng = np.random.default_rng(acfg.random_seed)
        vals = []
        for _ in range(acfg.n_boot_metrics):
            ri = rng.integers(0, n_rep, n_rep)
            si = rng.integers(0, n_sub, n_sub)
            pb = np.nanmean(M[ri][:, si], axis=0)
            yb = y[si]
            if len(np.unique(yb)) == 2 and np.isfinite(pb).all():
                vals.append(roc_auc_score(yb, pb))
        if vals:
            out["roc_auc_ci95_subjects_and_cv"] = [float(np.percentile(vals, 2.5)),
                                                   float(np.percentile(vals, 97.5))]
            out["roc_auc_ci_note"] = ("IC que incorpora a variabilidade das particoes da CV; "
                                      "e o intervalo a reportar.")
    return out


def permutation_auc_test(X: np.ndarray, y: np.ndarray, observed_matched: float,
                         acfg: AnalysisConfig) -> Dict[str, Any]:
    """B1: teste de permutacao com ESTIMADOR IDENTICO entre observado e nulo.

    ERRO CORRIGIDO. Na v23 a AUC observada vinha de ``outer_repeats`` = 20
    repeticoes de CV (media de probabilidades OOF, portanto menos ruidosa) e cada
    AUC nula de UMA repeticao. Como a media sobre repeticoes altera a distribuicao
    amostral da estatistica, observado e nulo deixavam de ser permutaveis e o
    p-valor perdia interpretacao, com vies anticonservador.

    Agora ambos usam ``permutation_matched_repeats``. ``observed_matched`` DEVE ser
    calculado com o mesmo numero de repeticoes (ver ``matched_observed_auc``).
    """
    reps = int(acfg.permutation_matched_repeats)
    rng = np.random.default_rng(acfg.random_seed)
    null = []
    for _ in range(acfg.n_permutations):
        yp = rng.permutation(y)
        r = nested_cv_scores(X, yp, acfg, repeats=reps)
        null.append(float(roc_auc_score(yp, r["prob_mean"])))
    null = np.asarray(null, float)
    p = float((np.sum(null >= observed_matched) + 1) / (len(null) + 1))
    return {"n_permutations": int(len(null)),
            "estimator_repeats_observed_and_null": reps,
            "observed_auc_matched_estimator": float(observed_matched),
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
            "null_p95": float(np.percentile(null, 95)), "p_value": p,
            "note": ("Observado e nulo usam o mesmo numero de repeticoes de CV; "
                     "o p-valor NAO deve ser comparado com a AUC de %d repeticoes."
                     % acfg.outer_repeats)}


def matched_observed_auc(X: np.ndarray, y: np.ndarray, acfg: AnalysisConfig) -> float:
    """B1: AUC observada sob o MESMO estimador usado na distribuicao nula."""
    r = nested_cv_scores(X, y, acfg, repeats=int(acfg.permutation_matched_repeats))
    return float(roc_auc_score(y, r["prob_mean"]))


def sign_agreement_with_ci(g_ref: pd.Series, g_alt: pd.Series, n_boot: int = 2000,
                           seed: int = 0) -> Dict[str, Any]:
    """B7: concordancia de sinal entre janelas COM incerteza, sem veredito binario.

    ERRO CORRIGIDO. A v23 rotulava "consistente entre janelas" quando
    sign_agreement >= 0,7 sobre 9 features. Sob a hipotese nula de sinais
    aleatorios, P(>= 7 de 9 concordancias) = 0,090: o rotulo aparece por acaso em
    ~9% das vezes. Um limiar arbitrario sobre poucas features nao e verificacao de
    robustez. Agora reportamos a proporcao, o IC bootstrap por feature e o p-valor
    binomial exato contra o acaso (p = 0,5), sem rotulo qualitativo.
    """
    common = g_ref.index.intersection(g_alt.index)
    a = pd.to_numeric(g_ref[common], errors="coerce").to_numpy(float)
    b = pd.to_numeric(g_alt[common], errors="coerce").to_numpy(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = int(a.size)
    if n == 0:
        return {"n_features": 0, "sign_agreement": float("nan"),
                "ci95": [float("nan"), float("nan")], "p_binomial_vs_chance": float("nan")}
    agree = (np.sign(a) == np.sign(b)).astype(float)
    rng = np.random.default_rng(seed)
    boot = np.array([agree[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    k = int(agree.sum())
    p_bin = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)
    rho = stats.spearmanr(a, b).statistic if n >= 3 else float("nan")
    return {"n_features": n, "n_agreeing": k, "sign_agreement": float(agree.mean()),
            "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "p_binomial_vs_chance": p_bin,
            "spearman_effects": float(rho) if rho == rho else float("nan"),
            "note": ("Sem rotulo qualitativo: com poucas features a concordancia alta "
                     "e frequente por acaso (P(>=7/9) = 0,090).")}


def batch_confound_report(df: pd.DataFrame, acfg: AnalysisConfig) -> Dict[str, Any]:
    """Mede a colinearidade entre lote (prefixo de ID) e diagnostico, e o poder do lote.

    Se o prefixo separa perfeitamente os grupos, qualquer diferenca tecnica entre
    lotes e matematicamente indistinguivel do efeito clinico. O relatorio expoe
    esse fato e, quando V >= limiar, marca todos os resultados como confundidos.
    """
    if "id_prefix" not in df.columns:
        return {"available": False}
    v = cramers_v(df["id_prefix"], df["label"])
    tab = pd.crosstab(df["id_prefix"], df["label"]).to_dict()
    perfect = bool(df.groupby("id_prefix")["label"].nunique().max() == 1)
    return {"available": True, "cramers_v": v, "crosstab": tab,
            "prefix_perfectly_predicts_label": perfect,
            "batch_confounded": bool(perfect or v >= acfg.batch_collinearity_warn),
            "interpretation": ("Prefixo de ID determina o diagnostico: efeito de lote e "
                               "efeito clinico sao INSEPARAVEIS neste conjunto."
                               if perfect else
                               "Colinearidade parcial entre lote e diagnostico; ajustar e "
                               "reportar sensibilidade.")}

# ==================================================================================
# BLOCO 9/11 - ORQUESTRACAO
# ==================================================================================

def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def environment_manifest() -> Dict[str, Any]:
    """Versoes de software que afetam a reproducibilidade numerica."""
    import scipy
    env = {"python": sys.version.split()[0], "platform": platform.platform(),
           "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
           "utc_created": datetime.now(timezone.utc).isoformat()}
    env["scikit_learn"] = __import__("sklearn").__version__ if _SKLEARN_OK else "not_installed"
    return env


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def extract_features_for_window(filtered: np.ndarray, rec: SubjectRecording, cfg: RunConfig,
                                start_s: float, window_s: float) -> Dict[str, Any]:
    """Recorta a janela, epoca, aplica QC de epoca e agrega as features do participante."""
    acq, pre, qc, fcfg = cfg.acquisition, cfg.preproc, cfg.qc, cfg.features
    win = extract_window(filtered, acq.fs, start_s, window_s)
    drift = within_window_drift(win, acq.fs, cfg.temporal.block_seconds)
    epochs = epoch_signal(win, acq.fs, pre.epoch_seconds, pre.epoch_overlap)
    f, psd = compute_psd(epochs, acq.fs, fcfg)
    good, diag = epoch_rejection_mask(epochs, f, psd, qc, fcfg)
    if diag["frac_rejected"] > qc.max_rejected_fraction:
        raise ValueError("%s: fracao de epocas rejeitadas %.3f > %.2f."
                         % (rec.subject_id, diag["frac_rejected"], qc.max_rejected_fraction))
    if int(good.sum()) < fcfg.min_good_epochs:
        raise ValueError("%s: %d epocas validas < %d."
                         % (rec.subject_id, int(good.sum()), fcfg.min_good_epochs))
    arrays = epoch_feature_arrays(epochs[good], acq.fs, fcfg, acq.channel_names)
    row: Dict[str, Any] = {"subject_id": rec.subject_id, "source_name": rec.source_name,
                           "task": rec.task, "window_start_s": float(start_s),
                           "window_seconds": float(window_s),
                           "n_epochs_window": int(len(epochs)),
                           "n_epochs_used": int(good.sum()),
                           "frac_epochs_rejected": float(diag["frac_rejected"])}
    row.update({("qc_" + k): v for k, v in diag.items() if k.startswith("n_rej")})
    row.update({("drift_" + k.replace("drift_", "")): v for k, v in drift.items()})
    row.update(aggregate_features(arrays, acq.channel_names, cfg.features))
    return row


def run_v24(eeg_dir: str | Path, metadata_path: str | Path, output_dir: str | Path,
            selection_criterion: str, policy: EvidencePolicy = EvidencePolicy(),
            cfg: RunConfig = RunConfig()) -> Dict[str, Any]:
    """Executa o pipeline completo: integracao -> protocolo temporal -> QC -> features -> analise.

    Artefatos em quatro pastas: 01_integration, 02_quality, 03_analysis, 04_sensitivity.
    Falhas por arquivo NAO abortam o lote; vao para failures.csv, agora segregadas
    por estagio (B6). O modo research so e liberado com limiares de QC congelados (B5).
    """
    if not isinstance(selection_criterion, str) or len(selection_criterion.strip()) < 10:
        raise ValueError("selection_criterion obrigatorio (>=10 caracteres).")
    assert_qc_thresholds_frozen(cfg.qc_freeze, policy)  # B5
    out = Path(output_dir)
    d_int, d_qual = out / "01_integration", out / "02_quality"
    d_ana, d_sens = out / "03_analysis", out / "04_sensitivity"
    for d in (d_int, d_qual, d_ana, d_sens):
        d.mkdir(parents=True, exist_ok=True)

    files, failures = discover_eeg_files(eeg_dir)
    metadata = load_modma_metadata(metadata_path)

    # --- Etapa 1: leitura, esquema, filtragem, QC e perfil temporal (cego ao rotulo) ---
    recs: Dict[str, SubjectRecording] = {}
    filtered_map: Dict[str, np.ndarray] = {}
    qc_rows, block_rows, fs_rows = [], [], []
    for fp in files:
        try:
            rec = load_modma_txt(fp, cfg.acquisition, fp.name, cfg.schema)
            if rec.subject_id in recs:
                raise ValueError("Sujeito duplicado: %s" % rec.subject_id)
            fsv = verify_sampling_rate(rec.data_counts, cfg.acquisition)  # B2
            fsv["subject_id"] = rec.subject_id
            filt, fdiag = preprocess_continuous(rec.data_counts, cfg.acquisition, cfg.preproc)
            prof = block_profile(filt, cfg.acquisition.fs, cfg.temporal.block_seconds)
            prof["subject_id"] = rec.subject_id
            recs[rec.subject_id] = rec
            filtered_map[rec.subject_id] = filt
            qc_rows.append(subject_qc_metrics(rec, filt, cfg.acquisition, cfg.features, fdiag))
            block_rows.append(prof)
            fs_rows.append(fsv)
        except Exception as exc:
            failures.append({"file": fp.name, "subject_id": fp.stem.split("_")[0],
                             "stage": "read_or_profile",
                             "error": "%s: %s" % (type(exc).__name__, str(exc)[:200])})
    if not recs:
        raise RuntimeError("Nenhum registro utilizavel.")

    # B3: referencia espectral EXTERNA, mediana da coorte nos blocos terminais
    cohort_ref = cohort_terminal_profile(block_rows, cfg.temporal)
    settle_rows = []
    for prof in block_rows:
        sid = str(prof["subject_id"].iloc[0])
        st = subject_settling_time(prof, cfg.temporal, cohort_ref_profile=cohort_ref)
        st["subject_id"] = sid
        st["duration_s"] = recs[sid].parse_report["duration_s"]
        settle_rows.append(st)

    qc_df = cohort_qc_decision(pd.DataFrame(qc_rows), cfg.qc)
    qc_dist = qc_distribution_report(qc_df, cfg.qc)                      # B5
    settle_df = pd.DataFrame(settle_rows)
    blocks_df = pd.concat(block_rows, ignore_index=True)
    fs_df = pd.DataFrame(fs_rows)

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
                rows.append(extract_features_for_window(filtered_map[sid], rec, cfg, s0, L))
            except Exception as exc:
                failures.append({"file": rec.source_name, "subject_id": sid,
                                 "stage": "features_w%.0f" % s0,
                                 "error": "%s: %s" % (type(exc).__name__, str(exc)[:200])})
        feat_by_window["w_%.0f_%.0f" % (s0, L)] = pd.DataFrame(rows)

    key_primary = "w_%.0f_%.0f" % (primary_start, primary_len)
    features = feat_by_window[key_primary]

    # --- Etapa 4: pareamento e coorte analitica ---
    merged = features.merge(metadata, on="subject_id", how="left")
    unmatched = merged.loc[merged["label"].isna(), "subject_id"].tolist()
    merged = merged.loc[merged["label"].notna()].copy()
    merged["label"] = merged["label"].astype(int)
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
               "cohort_audit": audit_expected_cohort(list(recs), cfg.cohort),
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
        conf = confirmatory_group_tests(analytic, acfg)
        conf.to_csv(d_ana / "confirmatory_tests.csv", index=False)
        adj = confound_adjusted_models(analytic, acfg)
        adj.to_csv(d_ana / "confound_adjusted_models.csv", index=False)
        batch = batch_confound_report(analytic, acfg)
        cols = primary_columns(analytic)
        X = analytic[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        y = analytic["label"].to_numpy(int)
        res = nested_cv_scores(X, y, acfg)
        met = metrics_from_probs(y, res["prob_mean"], acfg,
                                 oof_per_repeat=res["oof_prob_per_repeat"])  # B10
        covs = [c for c in acfg.covariates if c in analytic.columns]
        base = {"majority_class": {"balanced_accuracy": 0.5, "roc_auc": 0.5}}
        if covs:
            Xd = analytic[covs].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            if np.isfinite(Xd).any():
                rd = nested_cv_scores(Xd, y, acfg)
                base["demographic"] = metrics_from_probs(y, rd["prob_mean"], acfg,
                                                         rd["oof_prob_per_repeat"])
                rc = nested_cv_scores(np.hstack([X, Xd]), y, acfg)
                base["eeg_plus_demographic"] = metrics_from_probs(
                    y, rc["prob_mean"], acfg, rc["oof_prob_per_repeat"])
            else:
                base["demographic"] = {"error": "covariaveis sem valores finitos (ver B9)"}
        # B1: observado e nulo sob o MESMO estimador
        obs_matched = matched_observed_auc(X, y, acfg)
        perm = permutation_auc_test(X, y, obs_matched, acfg)

        # B8: TODAS as janelas com o mesmo numero de repeticoes da primaria
        reps_sens = int(acfg.sensitivity_outer_repeats or acfg.outer_repeats)
        sens = {}
        for k, dfw in feat_by_window.items():
            m = dfw.merge(metadata, on="subject_id", how="inner")
            m = m.loc[m["subject_id"].isin(analytic["subject_id"])]
            if len(m) < policy.min_subjects_for_ml or m["label"].nunique() < 2:
                continue
            ck = confirmatory_group_tests(m, acfg)
            Xk = m[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            yk = m["label"].to_numpy(int)
            rk = nested_cv_scores(Xk, yk, acfg, repeats=reps_sens)
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
            c = sign_agreement_with_ci(gref, pd.Series(v["effects"]),
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
                  indent=2, ensure_ascii=False, default=_json_default)

    manifest = {"version": __version__, "changelog": list(CHANGELOG_V24),
                "environment": environment_manifest(), "config": asdict(cfg),
                "policy": asdict(policy), "selection_criterion": selection_criterion,
                "input_files": [{"name": p.name, "sha256": sha256_file(p)} for p in files]}
    json.dump(summary, open(out / "run_summary.json", "w"), indent=2,
              ensure_ascii=False, default=_json_default)
    json.dump(manifest, open(out / "manifest.json", "w"), indent=2,
              ensure_ascii=False, default=_json_default)
    return {"summary": summary, "manifest": manifest, "features": merged,
            "analytic": analytic, "qc": qc_df, "qc_distribution": qc_dist,
            "fs_verification": fs_df, "settling": settle_df,
            "features_by_window": feat_by_window, "failures": fail_df}

# ==================================================================================
# BLOCO 10/11 - AUTOTESTES
# ==================================================================================

def _synth_recording(fs: float, seconds: float, alpha_hz: float = 10.0,
                     settle_s: float = 120.0, seed: int = 0) -> np.ndarray:
    """Gera sinal sintetico de 3 canais com transitorio inicial de amplitude."""
    rng = np.random.default_rng(seed)
    n = int(fs * seconds)
    t = np.arange(n) / fs
    env = 1.0 + 3.0 * np.exp(-t / max(settle_s / 3.0, 1.0))  # acomodacao decrescente
    base = (np.sin(2 * np.pi * alpha_hz * t) + 0.5 * rng.standard_normal(n))
    x = np.stack([env * (base + 0.1 * rng.standard_normal(n)) for _ in range(3)])
    return np.round(x * 1000.0)


def selftest_v24(verbose: bool = True) -> Dict[str, bool]:
    """Autotestes com respostas conhecidas, incluindo os novos testes B1-B11."""
    import tempfile
    res: Dict[str, bool] = {}
    cfg = RunConfig()
    fs = cfg.acquisition.fs

    # --- basicos herdados ---
    res["nome_com_task"] = parse_modma_filename("02010002_still.txt") == ("02010002", "still")
    res["nome_sem_task"] = parse_modma_filename("02010001.txt") == ("02010001", TASK_UNSPECIFIED)
    a, n = fix_integer_wraparound(np.array([[10.0, 4294967295.0]]), 32)
    res["wraparound"] = bool(n == 1 and a[0, 1] == -1.0)
    res["js_identico_zero"] = abs(jensen_shannon([.5, .5], [.5, .5])) < 1e-12
    res["js_disjunto_um"] = abs(jensen_shannon([1, 1e-12], [1e-12, 1]) - 1.0) < 1e-3
    res["cramers_v_colinear"] = abs(cramers_v(["a", "a", "b", "b"], [1, 1, 0, 0]) - 1.0) < 1e-9
    res["bh_fdr_monotono"] = bool(np.all(np.diff(np.sort(bh_fdr([0.001, 0.02, 0.3, 0.7]))) >= 0))

    # --- B2: verificacao empirica de fs ---
    t = np.arange(int(fs * 60)) / fs
    x50 = np.tile(1000 * (np.sin(2 * np.pi * 10 * t) + 3.0 * np.sin(2 * np.pi * 50 * t)), (3, 1))
    v50 = verify_sampling_rate(x50, cfg.acquisition)
    res["B2_detecta_pico_50hz"] = bool(v50["fs_supported"]
                                       and abs(v50["line_peak_hz_under_assumed_fs"] - 50) <= 0.6)
    x43 = np.tile(1000 * (np.sin(2 * np.pi * 10 * t) + 3.0 * np.sin(2 * np.pi * 43 * t)), (3, 1))
    res["B2_rejeita_pico_inesperado"] = bool(not verify_sampling_rate(x43,
                                                                     cfg.acquisition)["fs_supported"])
    res["B2_fs_declarado_suposicao"] = bool(cfg.acquisition.fs_is_assumption)

    # --- A5/B3: acomodacao e referencia externa ---
    x = _synth_recording(fs, 700.0, settle_s=150.0, seed=1)
    filt, fdiag = preprocess_continuous(x, cfg.acquisition, cfg.preproc)
    prof = block_profile(filt, fs, cfg.temporal.block_seconds)
    ref = cohort_terminal_profile([prof], cfg.temporal)
    res["B3_referencia_coorte_normalizada"] = bool(ref is not None
                                                   and abs(float(np.sum(ref)) - 1.0) < 1e-6)
    st = subject_settling_time(prof, cfg.temporal, cohort_ref_profile=ref)
    res["A5_detecta_transitorio"] = bool(np.isfinite(st["settling_time_s"])
                                         and st["settling_time_s"] > 30.0)
    res["B3_reporta_referencia"] = st["reference_used"] == "cohort_terminal_blocks"
    res["B3_flag_deriva_global"] = "drift_spans_record" in st
    y = _synth_recording(fs, 700.0, settle_s=0.5, seed=2)
    fl2, _ = preprocess_continuous(y, cfg.acquisition, cfg.preproc)
    st2 = subject_settling_time(block_profile(fl2, fs, cfg.temporal.block_seconds),
                                cfg.temporal, cohort_ref_profile=ref)
    res["A5_estavel_acomoda_cedo"] = bool(st2["settling_time_s"] <= st["settling_time_s"])

    # --- B4: cap explicito ---
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

    # --- B9: parser de sexo ---
    s_txt, r_txt = parse_sex_column(pd.Series(["M", "F", "male", "Female"]))
    res["B9_texto"] = list(s_txt) == [1.0, 0.0, 1.0, 0.0]
    s_num, r_num = parse_sex_column(pd.Series([1, 2, 1, 2]))
    res["B9_numerico_1M_2F"] = (list(s_num) == [1.0, 0.0, 1.0, 0.0]
                                and r_num["coding_detected"] == "numeric_1M_2F")
    res["B9_relata_nao_mapeado"] = parse_sex_column(pd.Series(["?", "?"]))[1]["n_mapped"] == 0

    # --- B7: concordancia com IC ---
    g1 = pd.Series({"a": 0.4, "b": -0.3, "c": 0.2, "d": -0.5})
    ag = sign_agreement_with_ci(g1, g1, n_boot=200, seed=1)
    res["B7_concordancia_total"] = abs(ag["sign_agreement"] - 1.0) < 1e-12
    res["B7_tem_ic_e_p"] = ("ci95" in ag and np.isfinite(ag["p_binomial_vs_chance"]))
    res["B7_discordancia"] = abs(sign_agreement_with_ci(g1, -g1, n_boot=200,
                                                        seed=1)["sign_agreement"]) < 1e-12

    # --- B1: estimadores pareados ---
    res["B1_reps_pareados"] = bool(AnalysisConfig().permutation_matched_repeats >= 1)
    res["B8_repeticoes_iguais"] = AnalysisConfig().sensitivity_outer_repeats is None

    # --- B5: distribuicao de QC e trava de congelamento ---
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

    # --- B11: coorte canonica ---
    res["B11_55_ids_unicos"] = len(set(EXPECTED_SUBJECT_IDS_55)) == 55
    aud = audit_expected_cohort(list(EXPECTED_SUBJECT_IDS_55[:54]), CohortConfig())
    res["B11_audita_faltantes"] = (len(aud["missing_ids"]) == 1 and not aud["unexpected_ids"])
    res["B11_threshold_configuravel"] = abs(AnalysisConfig().decision_threshold - 0.5) < 1e-12

    # --- features e barreiras ---
    ep = epoch_signal(extract_window(filt, fs, 300.0, 240.0), fs, 4.0)
    f, psd = compute_psd(ep, fs, cfg.features)
    good, diag = epoch_rejection_mask(ep, f, psd, cfg.qc, cfg.features)
    res["qc_epoca_maioria_boa"] = bool(good.mean() > 0.5)
    arrs = epoch_feature_arrays(ep[good][:6], fs, cfg.features, cfg.acquisition.channel_names)
    row_f = aggregate_features(arrs, cfg.acquisition.channel_names, cfg.features)
    res["nove_features_primarias"] = sum(k.startswith("p__") for k in row_f) == 9
    res["alfa_detectado"] = bool(abs(np.nanmedian(arrs["alpha_peak_hz"]) - 10.0) <= 1.5)
    res["barreira_primaria"] = False
    try:
        assert_primary_only(["p__rel_alpha", "e__faa"])
    except RuntimeError:
        res["barreira_primaria"] = True
    res["barreira_circular"] = False
    try:
        assert_no_circular_features(["p__rel_alpha", "phq9"])
    except RuntimeError:
        res["barreira_circular"] = True

    # --- execucao ponta a ponta em dados sinteticos ---
    with tempfile.TemporaryDirectory() as td:
        root = Path(td); eegd = root / "eeg"; eegd.mkdir()
        for i, sid in enumerate(["02010001", "02010002", "02020004", "02020007"]):
            np.savetxt(eegd / ("%s_still.txt" % sid),
                       _synth_recording(fs, 700.0, settle_s=120.0, seed=i).T, fmt="%d")
        np.savetxt(eegd / "02020018_still.txt", _synth_recording(fs, 72.8, seed=9).T, fmt="%d")
        pd.DataFrame({"subject id": [2010001, 2010002, 2020004, 2020007, 2020018],
                      "type": ["MDD", "MDD", "HC", "HC", "HC"],
                      "age": [30, 40, 35, 45, 50], "gender": [1, 2, 1, 2, 1],
                      "education（years）": [12, 16, 14, 18, 12]}).to_csv(root / "meta.csv",
                                                                        index=False)
        got = run_v24(eegd, root / "meta.csv", root / "out",
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

# ==================================================================================
# BLOCO 11/11 - CLI
# ==================================================================================

def main_v24() -> None:
    """CLI: python depressao_eeg_v24.py --eeg-dir DIR --metadata FILE --output DIR
    --selection-criterion "texto" [--mode research]
    """
    ap = argparse.ArgumentParser(description="MODMA EEG v24: correcoes B1-B11 de auditoria metodologica")
    ap.add_argument("--eeg-dir"); ap.add_argument("--metadata"); ap.add_argument("--output")
    ap.add_argument("--selection-criterion", default="")
    ap.add_argument("--mode", default="integration", choices=["integration", "research"])
    ap.add_argument("--temporal-mode", default="derived", choices=["derived", "fixed"])
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        ok = selftest_v24(True)
        raise SystemExit(0 if all(ok.values()) else 1)
    missing = [k for k in ("eeg_dir", "metadata", "output") if not getattr(args, k)]
    if missing:
        raise SystemExit("Faltam argumentos obrigatorios: %s" % missing)
    cfg = RunConfig(temporal=TemporalProtocolConfig(mode=args.temporal_mode))
    result = run_v24(args.eeg_dir, args.metadata, args.output, args.selection_criterion,
                     EvidencePolicy(mode=args.mode), cfg)
    print(json.dumps(result["summary"], indent=2, ensure_ascii=False, default=_json_default)[:4000])


if __name__ == "__main__":
    main_v24()
