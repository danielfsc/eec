
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

