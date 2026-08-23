
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
import re

_NAME_RE = re.compile(r"^(?P<sid>\d{8})(?:_(?P<task>[A-Za-z]+))?\.txt$")

TASK_UNSPECIFIED = "unspecified"
