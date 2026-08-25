# Guia Didático Completo — depressao_eeg_v24 / main

### Como funciona, do início ao fim, um estudo auditável de EEG frontal de três canais aplicado ao conjunto MODMA de depressão maior

*Manual escrito para leitores sem experiência em programação. Cada termo técnico é explicado quando aparece pela primeira vez e retomado no glossário final. Versão do pipeline: 24.0.0. Estado: 52 autotestes internos executados e aprovados.*

---

## Parte I — Resumo executivo

Existe uma pergunta simples por trás deste projeto: é possível olhar para a atividade elétrica do cérebro, registrada com apenas três eletrodos na testa, e distinguir quem tem depressão maior de quem não tem?

A resposta que este trabalho oferece não é "sim" nem "não". É uma terceira resposta, mais útil e mais honesta: com os dados disponíveis, essa pergunta não pode ser respondida — e nós conseguimos demonstrar exatamente por quê, com números.

O motivo central é o seguinte. Nos 55 registros disponíveis, o número de identificação de cada participante começa com um prefixo que indica em qual lote ele foi gravado. Acontece que todos os 26 participantes do lote que começa com **0201** são pacientes com depressão, e todos os 29 dos lotes **0202** e **0203** são controles saudáveis. A separação é perfeita. Isso significa que qualquer diferença técnica entre os lotes — um equipamento levemente diferente, uma sala com mais interferência elétrica, um operador com técnica distinta, uma época do ano diferente — produz exatamente o mesmo padrão nos dados que uma diferença cerebral real produziria. As duas explicações são matematicamente indistinguíveis. Nenhuma técnica estatística resolve isso, porque não sobra variação nos dados a partir da qual separar as duas causas.

Some-se a isso um segundo problema: com 26 pessoas em um grupo e 28 no outro, o estudo só consegue detectar diferenças muito grandes. Diferenças do tamanho que a literatura de EEG costuma encontrar em depressão simplesmente não apareceriam — e se aparecessem, viriam com tamanho exagerado, por um fenômeno estatístico conhecido.

E um terceiro: os três eletrodos estão tão próximos uns dos outros que gravam quase o mesmo sinal. A correlação mediana entre eles é **0,983**, onde **1,0** significaria sinais idênticos. Na prática, não temos três medidas independentes; temos aproximadamente uma.

Diante disso, o programa foi construído com uma filosofia incomum: em vez de tentar extrair o máximo de desempenho dos dados, ele foi projetado para tornar impossível esconder as limitações. O código recusa-se a produzir qualquer resultado estatístico antes que decisões críticas sejam registradas e datadas. Ele calcula e imprime o grau de confusão entre lote e diagnóstico. Ele deriva parâmetros dos próprios dados em vez de usar convenções arbitrárias. E ele guarda uma assinatura digital de cada arquivo de entrada, para que qualquer pessoa possa verificar depois que os dados não mudaram.

O produto final, portanto, não é um classificador de depressão. É um procedimento auditável e um mapa explícito do que este conjunto de dados permite e do que ele não permite concluir.

---

## Parte II — Visão geral: o caminho de um dado, do arquivo ao resultado

Antes de entrar em detalhes, é útil ter em mente o percurso completo. Imagine que você acompanha um único participante ao longo de todo o programa.

### II.1 O ponto de partida

Tudo começa com um arquivo de texto simples, chamado por exemplo `02010002_still.txt`. Dentro dele há centenas de milhares de linhas, cada uma com três números inteiros separados por espaço. Cada linha é um instante no tempo; cada número é a leitura de um dos três eletrodos naquele instante. "still" indica que o participante estava em repouso.

Esses números não são microvolts. São contagens de um conversor analógico-digital — o componente que transforma a voltagem contínua captada pelo eletrodo em um número inteiro que o computador consegue armazenar. Quantos microvolts vale cada unidade de contagem não está documentado publicamente. Essa é a primeira coisa que o pipeline declara abertamente, e ela tem uma consequência de projeto: todas as medidas realmente importantes do estudo foram escolhidas para não depender desse fator desconhecido.

### II.2 As nove etapas

| Etapa | O que acontece | Por que importa |
| :---: | :--- | :--- |
| **1** | **Descoberta e leitura** | O programa varre a pasta, identifica arquivos com nome válido, lê os números e verifica se o formato bate com o esperado.<br> *Um arquivo corrompido ou com formato diferente é isolado, não contamina o lote.* |
| **2** | **Correção de números negativos** | Alguns valores aparecem como números gigantescos (perto de 4,29 bilhões). São, na verdade, números negativos gravados de forma diferente. O programa os reinterpreta.<br> *Sem essa correção, o sinal teria degraus artificiais enormes — e os resultados sairiam errados sem que nada parecesse errado.* |
| **3** | **Verificação da taxa de amostragem** | O programa procura no sinal a assinatura da rede elétrica (50 Hz na China) para checar se a suposição sobre a velocidade de gravação é plausível.<br> *Se essa suposição estiver errada, todas as frequências e todos os tempos do estudo estarão errados por um fator constante.* |
| **4** | **Filtragem** | Remove componentes muito lentos e muito rápidos, e retira a interferência da rede elétrica quando ela é detectada.<br> *Isola a faixa de frequências onde está a atividade cerebral de interesse.* |
| **5** | **Perfil temporal** | O registro é dividido em blocos e o programa mede, bloco a bloco, quão forte e como está distribuído o sinal ao longo do tempo.<br> *Revela que o sinal muda sistematicamente nos primeiros minutos — o participante e os eletrodos estão se acomodando.* |
| **6** | **Derivação da janela de análise** | A partir desses perfis, o programa calcula quando o sinal estabiliza e escolhe automaticamente qual trecho de cada registro será analisado.<br> *A janela não é escolhida por convenção nem por conveniência: ela emerge dos dados, sem olhar quem é paciente e quem é controle.* |
| **7** | **Controle de qualidade** | Trechos e participantes com sinal ruim são identificados por critérios fixados de antemão, e o programa reporta quantos cada critério exclui.<br> *Impede que a decisão sobre quem fica de fora seja tomada depois de ver os resultados.* |
| **8** | **Extração de características** | Do trecho escolhido, o programa calcula nove medidas numéricas por participante.<br> *Transforma centenas de milhares de números em nove valores comparáveis entre pessoas.* |
| **9** | **Análise** | Só se autorizado: compara os grupos, testa um classificador, mede o confundimento de lote e avalia robustez.<br> *Produz os números do estudo — sempre acompanhados de suas incertezas e limitações.* |

### II.3 As duas travas de segurança

O programa tem dois estados de operação, e a diferença entre eles é o coração do desenho metodológico.

No **modo integração**, que é o padrão, o programa lê os dados, avalia a qualidade, deriva a janela, extrai as características — e para. Ele se recusa a calcular qualquer coisa que se pareça com um resultado: nada de acurácia, nada de p-valor, nada de curva de desempenho. A ideia é que você possa examinar exaustivamente os dados sem que exista a tentação de ajustar decisões olhando para o resultado final.

No **modo pesquisa**, a análise completa é executada. Mas ele não abre automaticamente: é preciso primeiro declarar formalmente, com data e justificativa registradas no código, que os critérios de qualidade foram fixados. Se você tentar entrar em modo pesquisa sem esse registro, o programa interrompe a execução e explica o que falta. Essa é uma trava deliberada contra um problema real e comum na ciência: ajustar critérios de exclusão até que o resultado fique bonito.

### II.4 O que sai no final

O programa cria quatro pastas e dois arquivos de registro. Nada é impresso apenas na tela — tudo fica gravado em disco, em formatos que qualquer pessoa pode abrir em uma planilha.

| Pasta ou arquivo | Conteúdo |
| :--- | :--- |
| `01_integration` | Características extraídas de cada participante, e a lista de tudo que falhou, com o motivo e em qual etapa. |
| `02_quality` | Métricas de qualidade por participante, distribuição empírica de cada critério, verificação da taxa de amostragem, perfis temporais e tempos de estabilização. |
| `03_analysis` | Só existe em modo pesquisa: comparações entre grupos, modelos ajustados, desempenho do classificador, teste de permutação. |
| `04_sensitivity` | As mesmas análises repetidas em outros trechos do registro, para checar se os achados dependem do trecho escolhido. |
| `run_summary.json` | Resumo geral em formato de texto estruturado: contagens, protocolo derivado, limitações declaradas. |
| `manifest.json` | Versões de todos os programas usados e a assinatura digital de cada arquivo de entrada. |

---

## Parte III — O estudo: pergunta, premissas, hipóteses

### III.1 A pergunta de pesquisa

A pergunta principal não é *"o EEG detecta depressão?"*. É:

> **Em que medida este conjunto específico de dados permite distinguir um sinal cerebral associado ao diagnóstico de depressão de uma variação que vem do lote de gravação, da qualidade do registro ou das escolhas feitas por quem analisa?**

Essa formulação é deliberada. Perguntar *"o EEG detecta depressão?"* convida a uma resposta em forma de porcentagem de acerto, que soa objetiva mas pode ser inteiramente produzida por artefatos. Perguntar *"o que estes dados permitem concluir?"* força o estudo a medir primeiro seus próprios limites.

Três perguntas secundárias derivam dela:
1. É possível determinar, a partir dos próprios dados e sem olhar quem é paciente, o momento em que cada registro se estabiliza — substituindo a prática comum de simplesmente descartar os primeiros trinta segundos?
2. As diferenças entre os grupos, se existirem, permanecem quando analisamos um trecho diferente do mesmo registro?
3. Um modelo que usa apenas idade, sexo e escolaridade acerta tanto quanto o modelo que usa o EEG? Se sim, o EEG não está acrescentando informação.

### III.2 As premissas — o que estamos assumindo sem poder provar

Toda análise repousa sobre suposições. A diferença entre uma análise honesta e uma enganosa costuma ser se essas suposições foram declaradas. Aqui estão as nossas.

| Premissa | Situação | O que acontece se estiver errada |
| :--- | :--- | :--- |
| **O sinal foi gravado a 250 amostras por segundo** | **SUPOSIÇÃO**, com verificação empírica pela frequência da rede elétrica. A documentação disponível descreve repouso de cerca de 5 minutos, mas os arquivos têm cerca de 20 minutos — a fonte é contradita pelos próprios dados em um ponto adjacente, e por isso não é tratada como autoridade. | Todas as frequências e todos os tempos ficam multiplicados por um fator constante. Nenhum resultado sobrevive, e o erro não gera mensagem de erro: os números continuam plausíveis. |
| **As três colunas do arquivo são, nesta ordem, Fp1, Fpz e Fp2** | **SUPOSIÇÃO** baseada na documentação do aparelho. | Medidas por canal e qualquer noção de assimetria ficam invertidas ou embaralhadas. |
| **Todos os 55 arquivos vêm do mesmo experimento** | **PREMISSA DECLARADA**, verificada indiretamente pela consistência de formato de cada arquivo. | Misturar protocolos diferentes introduziria variação que seria interpretada como diferença entre pessoas. |
| **Cada arquivo corresponde a exatamente um participante** | **VERIFICADO** pelo programa a cada execução. | Duplicatas inflariam artificialmente o tamanho da amostra e quebrariam a independência entre observações. |
| **Os rótulos de diagnóstico na planilha estão corretos** | **ASSUMIDO**; não temos como verificar. | Erros de rótulo reduzem qualquer diferença real e empurram os resultados para a ausência de efeito. |
| **A mediana da coorte corresponde a um valor típico de EEG frontal em repouso** | **ÂNCORA DECLARADA**, usada apenas para relatar amplitudes de forma descritiva. | Todos os valores em microvolts escalam junto. Como o fator é único para todos, nenhuma comparação entre grupos é afetada. |

### III.3 O que exatamente será testado

A tentação, em estudos assim, é calcular dezenas de medidas e reportar as que deram certo. Isso se chama pesca de resultados, e produz achados que não se replicam. Para evitá-lo, o conjunto de medidas foi congelado antes de qualquer análise, e dividido em duas famílias com estatutos diferentes.

A **família confirmatória** tem exatamente nove medidas. São elas que entram nos testes estatísticos e no classificador. Todas foram escolhidas por serem adimensionais — isto é, proporções ou frequências em hertz, que não dependem do fator de conversão desconhecido.

| Medida | O que significa, em linguagem comum |
| :--- | :--- |
| **Potência relativa em delta (1–4 Hz)** | Que fração da energia do sinal está nas oscilações mais lentas. Aumenta com sonolência e com movimento dos olhos. |
| **Potência relativa em teta (4–8 Hz)** | Fração nas oscilações lentas-intermediárias. Associada a estados de sonolência e, em parte da literatura, a estados afetivos. |
| **Potência relativa em alfa (8–13 Hz)** | Fração no ritmo característico do repouso relaxado, especialmente de olhos fechados. |
| **Potência relativa em beta (13–30 Hz)** | Fração nas oscilações rápidas, associadas a alerta e também à tensão muscular. |
| **Entropia espectral** | Quão espalhada é a energia entre as frequências. Valor alto significa um espectro plano, sem ritmo dominante; valor baixo significa um pico bem definido. |
| **Frequência do pico alfa** | Em que frequência exata está o pico do ritmo alfa, quando ele existe. Varia entre pessoas e muda com idade e estado. |
| **Inclinação do componente aperiódico** | Todo EEG tem um "fundo" que decai com a frequência, sem oscilação. A inclinação desse fundo é considerada um indicador do equilíbrio entre excitação e inibição no córtex. |
| **Mobilidade de Hjorth** | Uma medida de quão rápido o sinal varia, calculada diretamente no tempo, sem usar espectro. |
| **Complexidade de Hjorth** | Quão irregular é essa variação. Valor próximo de 1 indica um sinal parecido com uma onda pura. |

A **família exploratória** inclui entropia de permutação, dimensão fractal de Higuchi, complexidade C0 e assimetria frontal alfa. Elas são calculadas e reportadas, mas em arquivo separado, e existe no código uma verificação que interrompe a execução caso alguma delas tente entrar no modelo confirmatório. São descritivas: geram hipóteses para estudos futuros, não conclusões neste.

*Uma observação sobre a assimetria frontal alfa, porque ela é popular na literatura de depressão:* a versão validada compara os eletrodos F3 e F4, que ficam mais atrás na cabeça. Aqui temos Fp1 e Fp2, na testa. Além disso, a medida depende de qual ponto foi usado como referência elétrica, e essa informação não consta da documentação. Por essas duas razões ela permanece exploratória.

### III.4 As hipóteses

Cinco hipóteses foram registradas antes da execução. Quatro delas são formuladas na direção de "não haverá efeito detectável" — uma escolha deliberada, justificada pela análise de poder que vem a seguir.

| ID | Hipótese | Como será testada |
| :---: | :--- | :--- |
| **H1** | O prefixo do identificador prevê o diagnóstico com força alta. | Tabela cruzada e V de Cramér, uma medida de associação entre duas variáveis categóricas que vai de 0 (nenhuma) a 1 (perfeita). |
| **H2** | Os registros apresentam um transitório inicial, com queda sistemática de amplitude ao longo dos primeiros minutos. | Correlação entre amplitude e tempo, calculada participante a participante. |
| **H3** | Nenhuma das nove medidas confirmatórias difere entre os grupos após correção para múltiplos testes. | Teste de Mann-Whitney, tamanho de efeito de Hedges g e controle de taxa de falsas descobertas. |
| **H4** | O classificador baseado em EEG não supera de forma relevante o classificador baseado apenas em dados demográficos. | Validação cruzada aninhada, comparando os dois modelos com intervalos de confiança. |
| **H5** | O teste de permutação não rejeita a hipótese de ausência de sinal. | Embaralhamento dos rótulos milhares de vezes, com o mesmo procedimento de estimação. |

### III.5 A análise de poder — por que esperamos não encontrar nada

Poder estatístico é a probabilidade de um estudo detectar um efeito que realmente existe. É o oposto do risco de deixar passar. Um estudo com poder de 80% encontra o efeito em quatro de cada cinco tentativas; um com 20% de poder o encontra em uma de cada cinco.

O poder depende de três coisas: o tamanho da amostra, o tamanho do efeito real e o nível de exigência do teste. Com 26 pessoas em um grupo e 28 no outro, os números são estes:

| Situação | Tamanho do efeito (Hedges g) | Poder |
| :--- | :---: | :---: |
| **Menor efeito que conseguimos detectar de forma confiável** | 0,78 | 80% |
| **O mesmo, mas exigindo correção para as nove medidas testadas** | 1,02 | 80% |
| **Efeito médio, comum em pesquisa comportamental** | 0,50 | 44% |
| **Efeito pequeno-médio, faixa realista para EEG de repouso em depressão** | 0,30 | 19% |

Para dar sentido a esses números: um Hedges g de 0,78 significa que a diferença entre as médias dos dois grupos é quase oitenta por cento de um desvio padrão. Isso é uma diferença grande — o tipo de coisa que às vezes se enxerga a olho nu em um gráfico. Efeitos dessa magnitude são raros em EEG de repouso aplicado a transtornos psiquiátricos.

A conclusão a extrair disso, e que foi registrada antes de qualquer resultado, é desconfortável mas importante: **se este estudo encontrar um efeito estatisticamente significativo, é mais provável que ele seja uma flutuação amostral favorável do que um efeito real de magnitude modesta.** Esse fenômeno tem nome — *erro de magnitude*, ou *erro tipo M* — e é uma consequência matemática de estudos com poder baixo. Em amostras pequenas, apenas as estimativas exageradas conseguem cruzar o limiar de significância. Por isso, um achado significativo aqui deve ser tratado com mais suspeita, não com mais entusiasmo.

---

## Parte IV — Como o código funciona, função por função

Esta parte percorre o programa na ordem em que ele foi organizado. Você não precisa saber programar para acompanhar: cada função é apresentada como uma pequena máquina, com uma entrada, uma saída e uma razão de existir.

Antes de começar, três conceitos que aparecerão o tempo todo:
- **Função** é um trecho de programa com nome, que recebe informação, faz uma coisa específica e devolve um resultado. É a unidade básica de organização.
- **Configuração** é um conjunto de parâmetros — números e opções — reunidos em um só lugar. Neste programa todas as configurações são imutáveis: uma vez criadas, não podem ser alteradas durante a execução. Isso garante que ninguém mude um limiar no meio do caminho.
- **Docstring** é o texto explicativo que acompanha cada função dentro do próprio código. Neste projeto, as docstrings também registram as limitações metodológicas, de modo que a justificativa viaje junto com o código.

### IV.1 Constantes e listas de referência

| Nome | Função |
| :--- | :--- |
| `EXPECTED_SUBJECT_IDS_55` | A lista dos 55 identificadores que se espera encontrar. Permite dizer não apenas "faltou um arquivo", mas exatamente qual faltou. |
| `PRIMARY_FEATURES` | Os nomes das nove medidas confirmatórias. Serve de referência para as verificações de segurança. |
| `CLINICAL_SCALE_TOKENS` | Fragmentos de nome de escalas clínicas (phq, gad, psqi e outras). Usado para impedir que uma escala de sintomas seja usada como preditor — o que seria circular, já que o diagnóstico deriva delas. |
| `BLOCKED_ANALYSES` | A lista explícita do que o modo integração se recusa a fazer. |
| `OPEN_LIMITATIONS` | Limitações que o código não resolve, gravadas em todo relatório produzido. |
| `QC_THRESHOLD_MAP` | Mapeia cada métrica de qualidade ao parâmetro que a controla e ao sentido da comparação. Permite gerar automaticamente o relatório de impacto dos critérios. |

### IV.2 As configurações

Um bloco inteiro do programa contém apenas parâmetros. Reuni-los assim tem uma vantagem prática: para saber exatamente como uma execução foi feita, basta ler esse bloco. E tem uma vantagem de auditoria: o programa calcula uma impressão digital criptográfica de toda a configuração, um código curto que muda se qualquer parâmetro mudar. Duas execuções com a mesma impressão digital usaram exatamente os mesmos parâmetros.

| Configuração | O que controla |
| :--- | :--- |
| `AcquisitionConfig` | Taxa de amostragem, nomes e ordem dos canais, número de bits do formato. Registra explicitamente que a taxa é uma suposição. |
| `SchemaConfig` | O contrato que todo arquivo deve cumprir: três colunas, valores finitos, valores inteiros, duração dentro de faixa plausível. |
| `TemporalProtocolConfig` | As tolerâncias e limites da derivação da janela de análise. |
| `PreprocConfig` | Filtros: faixa de frequências, ordem, quanto cortar nas bordas, duração das fatias de análise. |
| `QCConfig` | Todos os critérios de qualidade, por trecho e por participante. |
| `QCFreezeConfig` | O registro formal do congelamento dos critérios: se foram congelados, quando e com base em quê. |
| `ScaleInferenceConfig` | Parâmetros da inferência descritiva de escala e de referência elétrica. |
| `CohortConfig` | A definição da coorte-alvo, incluindo a lista dos identificadores esperados. |
| `FeatureConfig` | Como cada característica é calculada: faixas de frequência, parâmetros do espectro, número mínimo de trechos válidos. |
| `AnalysisConfig` | O plano estatístico: número de repetições, grade de hiperparâmetros, número de permutações, nível de significância. |
| `EvidencePolicy` | O contrato de evidência: modo de operação e requisitos mínimos de tamanho de amostra. |
| `RunConfig` | Reúne todas as anteriores e verifica coerência global — por exemplo, que nenhuma frequência de filtro ultrapasse metade da taxa de amostragem, o que produziria resultados sem sentido. |

### IV.3 Leitura dos arquivos e dos metadados

- `parse_modma_filename` — Recebe o nome de um arquivo e extrai dele o identificador do participante e a tarefa. Aceita tanto nomes com sufixo (`02010002_still.txt`) quanto sem (`02010001.txt`), porque ambos os formatos existem no conjunto. Se um único arquivo tivesse o nome recusado, todo o lote poderia ser silenciosamente reduzido.
- `fix_integer_wraparound` — Esta função corrige um problema sutil e perigoso. Computadores guardam números inteiros em uma quantidade fixa de bits. Quando um número negativo é gravado sem indicar o sinal, ele reaparece na leitura como um número positivo enorme: o valor menos um vira 4.294.967.295. A função detecta esses valores e os traduz de volta para negativos, contando quantas amostras foram corrigidas.<br>*Por que isso é perigoso e não apenas chato:* sem a correção, o sinal teria saltos de mais de quatro bilhões de unidades. Os filtros processariam esses saltos normalmente, e as medidas resultantes continuariam sendo números finitos, de aparência razoável. O erro não se anuncia — ele se disfarça de resultado.
- `load_modma_txt` — Lê o arquivo inteiro, confere o formato contra o contrato definido em `SchemaConfig`, aplica a correção acima e devolve um objeto organizado com o sinal e um relatório de leitura. Se algo não bate, gera um erro descritivo em vez de prosseguir com dados suspeitos.
- `canonical_subject_id` — Normaliza identificadores para oito dígitos com zeros à esquerda. Isso resolve um problema clássico e traiçoeiro: planilhas eletrônicas costumam converter "02010002" no número 2010002, comendo o zero inicial. Se isso acontece, o cruzamento entre EEG e metadados falha em silêncio — não dá erro, apenas devolve uma tabela vazia ou incompleta.
- `discover_eeg_files` — Varre a pasta e separa arquivos com nome válido dos demais. Os inválidos vão para uma lista de quarentena, com o motivo registrado. Nenhum arquivo problemático interrompe o processamento dos outros.
- `load_modma_metadata` — Lê a planilha com diagnóstico, idade, sexo e escolaridade. Identifica as colunas mesmo que os nomes variem, converte "MDD" e "HC" em 1 e 0, e recusa a execução se algum rótulo não puder ser interpretado — em vez de assumir um valor por conta própria.
- `parse_sex_column` — Interpreta a coluna de sexo, que pode vir de várias formas: como texto ("M", "F", "male", "Female") ou como número (1 e 2). A função detecta a codificação, informa qual detectou, conta quantas linhas conseguiu mapear e avisa quando não conseguiu mapear nenhuma.<br>Vale explicar por que uma função inteira para algo tão simples. Se a coluna de sexo for lida incorretamente, ela vira uma coluna de valores ausentes. Os modelos que ajustam para idade, sexo e escolaridade então falham ou perdem a covariável — e isso pode acontecer sem nenhuma mensagem de erro. O relatório explícito torna a falha visível.
- `audit_expected_cohort` — Compara os identificadores efetivamente lidos com a lista esperada e devolve quais faltam e quais apareceram sem estar previstos. É o insumo do fluxograma de participantes.

### IV.4 Análise espectral e verificação da taxa de amostragem

Um parênteses conceitual, porque tudo a seguir depende disso. Qualquer sinal que varia no tempo pode ser descrito de duas maneiras equivalentes: como uma sequência de valores ao longo do tempo, ou como uma soma de oscilações de diferentes frequências. A segunda descrição chama-se **espectro**. É como decompor um acorde musical nas notas que o compõem. O EEG é analisado quase sempre no espectro, porque as faixas de frequência têm significado fisiológico conhecido.

- `psd_continuous` — Calcula a densidade espectral de potência: quanta energia o sinal tem em cada frequência. Usa o método de Welch, que divide o sinal em pedaços, calcula o espectro de cada um e faz a média. Isso produz uma estimativa mais estável, ao custo de um pouco de resolução em frequência.
- `band_share` — Calcula que fração da energia total está dentro de uma faixa de frequência. O resultado é uma proporção entre zero e um, portanto adimensional: não depende do fator de conversão desconhecido. É a operação básica de quase todas as medidas do estudo.
- `line_noise_ratio` — Mede quanta energia está concentrada em torno de 50 Hz. A rede elétrica irradia nessa frequência (na China e na maior parte do mundo; nos Estados Unidos são 60 Hz) e contamina registros de EEG. Esta medida é a evidência que decide se o filtro específico para rede elétrica será aplicado.
- `verify_sampling_rate` — Uma verificação independente e engenhosa. Sob a suposição de que o sinal foi gravado a 250 amostras por segundo, a interferência de rede deve aparecer exatamente em 50 Hz. A função procura o pico dominante entre 35 e 70 Hz e classifica o resultado:

| Resultado | Interpretação | O que fazer |
| :--- | :--- | :--- |
| **Pico na frequência esperada** | Evidência consistente com a suposição. | Não é prova, mas corrobora. Prossiga. |
| **Pico em frequência inesperada** | Alerta sério. A função calcula qual seria a taxa de amostragem se aquele pico fosse mesmo a rede. | Investigar antes de qualquer análise. |
| **Nenhum pico detectado** | Inconclusivo — o aparelho pode ter filtro interno de rede. | Não refuta a suposição, mas a deixa sem corroboração. Declare isso no relatório. |

### IV.5 Filtragem e segmentação

- `preprocess_continuous` — O coração do pré-processamento. Aplica, nesta ordem: o filtro de rede elétrica quando há evidência que o justifique; um filtro passa-banda entre 1 e 40 Hz; o corte de dois segundos em cada extremidade; e a remoção da média.<br>Cada uma dessas escolhas tem uma razão. O filtro é aplicado ao sinal contínuo, antes de fatiá-lo, porque filtrar cada fatia separadamente criaria descontinuidades artificiais nas bordas. O corte das extremidades remove o transitório do próprio filtro — filtros digitais produzem uma perturbação no início e no fim, que contaminaria as primeiras e últimas fatias. E o filtro é aplicado em duas passadas, para frente e para trás, o que elimina distorção de fase: sem esse cuidado, diferentes frequências sairiam deslocadas no tempo umas em relação às outras.
- `extract_window` — Recorta o trecho de análise definido pelo protocolo temporal.
- `epoch_signal` — Divide o trecho em fatias de quatro segundos. As fatias são disjuntas, sem sobreposição. Fatias sobrepostas compartilhariam dados e deixariam de ser observações independentes, o que faria os testes estatísticos parecerem mais precisos do que realmente são.

### IV.6 O protocolo temporal derivado dos dados

Esta é a parte metodologicamente mais original do pipeline, e merece uma explicação com calma.

**O problema.** Um registro de EEG em repouso não é homogêneo. Nos primeiros minutos, o participante ainda está se acomodando, e o gel condutor entre o eletrodo e a pele ainda está estabilizando sua condutividade elétrica. O sinal muda sistematicamente. Nestes 55 registros, a amplitude cai visivelmente: a correlação entre amplitude e tempo é negativa na maioria, com valor mediano em torno de **−0,55** e **31 dos 55** abaixo de **−0,5**.

A prática habitual é descartar os primeiros trinta segundos e analisar um trecho fixo. Mas trinta segundos é um número escolhido por tradição, não medido. Se a acomodação demora dois minutos, esse trecho está justamente dentro do transitório.

**A solução implementada.** Em vez de convencionar, o programa mede.

- `jensen_shannon` — Compara dois espectros e devolve um número entre zero e um indicando quão diferentes eles são. Zero significa idênticos; um significa que não compartilham nada. É a régua usada para dizer se o conteúdo em frequência mudou.
- `block_profile` — Divide o registro em blocos e mede, em cada um, a amplitude e a distribuição de energia entre as faixas de frequência. O resultado é um retrato de como o sinal evolui ao longo dos vinte minutos.
- `cohort_terminal_profile` — Calcula um espectro de referência externo: a mediana, entre todos os participantes, dos blocos finais de cada registro. É calculated uma única vez, para toda a coorte, e sem olhar quem é paciente.
- *Por que a referência precisa ser externa, e não interna:* Para dizer que um registro "estabilizou", é preciso comparar o começo com algum estado de referência. O caminho intuitivo seria usar a segunda metade do próprio registro. Mas quando existe uma deriva que atravessa todo o registro, essa segunda metade também está se movendo. Comparar o início com um alvo em movimento faz o método declarar estabilidade cedo demais — e, por construção, torna impossível detectar deriva global. Usar a coorte como referência resolve isso para o espectro, porque proporções de energia são adimensionais e comparáveis entre pessoas.<br>*Um limite honesto:* a referência de amplitude continua sendo interna a cada participante. Amplitude em contagens depende do quão bem o eletrodo fez contato com a pele daquela pessoa naquele dia, e não é comparável entre indivíduos sem calibração. Consequência: deriva que abrange o registro inteiro não é corrigida — ela é sinalizada, e o tempo de estabilização desses participantes deve ser lido como um valor mínimo, não exato.
- `subject_settling_time` — Percorre os blocos de trás para frente e encontra o primeiro instante a partir do qual todos os blocos seguintes estão simultaneamente dentro de duas tolerâncias: a de amplitude e a de forma espectral. Esse instante é o tempo de estabilização daquele participante. Devolve também qual referência foi usada e se há indício de deriva atravessando todo o registro.
- `derive_temporal_protocol` — Transforma os tempos individuais em uma única janela comum a todos. O início é um quantil dos tempos individuais, arredondado para a grade de blocos. A duração é o maior valor que ainda cabe em todos os registros elegíveis. Participantes com registro curto demais são excluídos, e a lista de excluídos é reportada.<br>Esta função também expõe algo que muitos programas escondem. Existe um limite máximo administrativo para o quanto se pode descartar do início. Se o cálculo produzir um valor maior que esse limite, o programa não corta em silêncio: ele registra que houve truncamento, emite um aviso e marca o protocolo como não válido sob a regra declarada. A razão é direta — se o protocolo foi truncado, a janela pode ter voltado para dentro do transitório, que é exatamente o problema que a derivação existia para evitar.
- `within_window_drift` — Verifica se, mesmo dentro da janela escolhida, ainda há mudança sistemática. Participantes com deriva residual acima do tolerado são excluídos da coorte analítica.

### IV.7 Controle de qualidade

A qualidade é avaliada em três níveis, todos cegos ao diagnóstico.

- `epoch_rejection_mask` — **Nível da fatia.** Marca fatias problemáticas usando critérios relativos ao próprio participante: amplitude muito acima ou muito abaixo da mediana das fatias daquela pessoa, saltos bruscos, canal totalmente parado, ou excesso de energia em frequências muito baixas — a assinatura típica de uma piscada, especialmente forte em eletrodos na testa.<br>Todos os critérios são razões, não valores absolutos. Isso os torna independentes da escala e, portanto, imunes ao fator de conversão desconhecido. Mas traz uma limitação que precisa ser dita: critérios internos ao participante não conseguem identificar um registro que seja ruim do início ao fim. Se tudo está igualmente ruim, nada se destaca como atípico. É para isso que existem os outros dois níveis.
- `nonstationarity_cv` — Mede quanto a potência do sinal oscila entre blocos. Um valor alto indica um registro instável.
- `subject_qc_metrics` — **Nível do participante.** Calcula, para o registro inteiro, um conjunto de indicadores: ruído de rede, índice ocular, razão muscular, instabilidade, fração de amostras congeladas, saturação do conversor e correlação entre canais.
- `cohort_qc_decision` — **Nível da coorte.** Combina os limiares absolutos com uma detecção de valores atípicos baseada em desvio absoluto mediano — uma medida de dispersão que, ao contrário do desvio padrão, não é distorcida pelos próprios valores extremos que se quer detectar. Produz, para cada participante, se ele passa e por quais motivos eventualmente não passa.
- `qc_distribution_report` — Talvez a função mais importante para a integridade do estudo, embora não calcule nada sofisticado. Para cada métrica de qualidade, ela monta uma tabela com os valores realmente observados — mínimo, quartis, mediana, máximo — ao lado do limiar vigente, e informa quantos participantes aquele limiar excluiria.<br>*A razão de existir:* um limiar que nunca foi comparado com a distribuição real é um número no escuro. Ele pode excluir ninguém, e então é decorativo; ou pode excluir quarenta pessoas, e então ele é o estudo inteiro. Descobrir isso depois de ver os resultados abre a porta para ajustá-lo até que o resultado agrade. Esta tabela obriga a decisão a ser tomada antes, olhando apenas a distribuição.
- `assert_qc_thresholds_frozen` — **A trava.** Se o modo pesquisa for solicitado sem que o congelamento tenha sido declarado, o programa interrompe a execução e explica exatamente o que precisa ser feito. E o registro de congelamento recusa ser criado sem uma data: congelamento sem data não é congelamento.

### IV.8 Extração das características

- `compute_psd` — Calcula o espectro de cada fatia válida.
- `relative_band_powers` — Converte o espectro nas quatro potências relativas.
- `spectral_entropy` — Mede o espalhamento do espectro. Trata a distribuição de energia como se fosse uma distribuição de probabilidade e calcula sua entropia, normalizada entre zero e um.
- `alpha_peak` — Encontra a frequência do pico alfa. Faz isso em duas etapas: primeiro remove a tendência de fundo do espectro, deixando apenas os picos aparentes; depois procura o pico mais alto na faixa de 7 a 13 Hz, exigindo altura mínima. Se nenhum pico satisfaz o critério, devolve valor ausente em vez de escolher um número arbitrário.<br>Essa recusa em inventar um valor é uma decisão de projeto importante. Uma versão ingênua simplesmente pegaria a frequência de maior energia na faixa — e devolveria um número mesmo quando não existe pico algum, contaminando a análise com valores sem significado.
- `aperiodic_slope` — Ajusta uma reta ao espectro em escala logarítmica, para estimar a inclinação do fundo aperiódico. O ajuste é robusto: identifica os picos oscilatórios, exclui-os, reajusta, e repete algumas vezes. Sem isso, o pico alfa puxaria a reta e distorceria a estimativa.
- `hjorth_parameters` — Calcula mobilidade e complexidade diretamente na série temporal, usando as variâncias do sinal e de suas derivadas. São medidas rápidas e independentes do espectro.
- `perm_entropy`, `higuchi_fd`, `c0_complexity` — Três medidas de complexidade não linear, todas exploratórias. A primeira olha padrões de ordenação entre valores vizinhos; a segunda estima quanto o traçado do sinal "preenche" o plano; a terceira compara o sinal com uma versão dele mesmo com as componentes regulares removidas.
- `epoch_feature_arrays` — Aplica todas as medidas a todas as fatias válidas, de todos os canais.
- `aggregate_features` — Resume para um valor por participante, usando a mediana. A mediana é preferida à média porque resiste melhor a fatias atípicas que tenham escapado do controle de qualidade.
- `assert_primary_only` e `assert_no_circular_features` — Duas verificações de segurança que interrompem a execução se violadas. A primeira impede que uma medida exploratória entre no modelo confirmatório. A segunda impede que escalas clínicas ou o próprio diagnóstico sejam usados como preditores — o que produziria um resultado espetacular e completamente vazio, já que o diagnóstico deriva dessas escalas.

### IV.9 Estatística

- `hedges_g` — Calcula o tamanho do efeito: quão distantes estão as médias dos dois grupos, expressa em desvios padrão. É corrigido para amostras pequenas. Diferentemente de um p-valor, o tamanho do efeito não depende do tamanho da amostra, e por isso é comparável entre estudos.
- `bootstrap_ci` — Estima intervalos de confiança pelo método de reamostragem. A ideia é simples e poderosa: sorteia-se repetidamente, com reposição, um novo conjunto de participantes a partir do conjunto original, recalcula-se a estatística milhares de vezes, e observa-se o quanto ela varia. A faixa que contém 95% dos valores obtidos é o intervalo de confiança. A reamostragem é feita por participante, porque o participante é a unidade independente do estudo.
- `bh_fdr` — Corrige para múltiplos testes. Quando se testam nove medidas, a chance de pelo menos um resultado parecer significativo por acaso é bem maior do que 5%. O procedimento de Benjamini-Hochberg controla a proporção esperada de falsos positivos entre os resultados declarados significativos.
- `cramers_v` — Mede a associação entre duas variáveis categóricas, de 0 a 1. É a função que quantifica o confundimento entre lote e diagnóstico.
- `confirmatory_group_tests` — Aplica, às nove medidas, o teste de Mann-Whitney, calcula os tamanhos de efeito com intervalos e aplica a correção para múltiplos testes. O teste de Mann-Whitney compara as distribuições sem supor que sejam normais, o que é prudente em amostras pequenas.
- `confound_adjusted_models` — Ajusta modelos incluindo idade, sexo e escolaridade, para verificar se alguma diferença observada persiste depois de levar em conta essas variáveis.
- `sign_agreement_with_ci` — Avalia se os achados se mantêm quando mudamos o trecho analisado. Devolve a proporção de medidas que apontam na mesma direção, o intervalo de confiança dessa proporção, o p-valor de um teste binomial exato contra o acaso, e a correlação entre os tamanhos de efeito.<br>Note que a função não emite nenhum veredito do tipo "os resultados são robustos". Isso é intencional. Com nove medidas, obter sete ou mais concordâncias por puro acaso tem probabilidade de **0,090** — quase um em onze. Um selo de robustez que aparece em nove por cento das execuções sem sinal algum não é informação, é ruído com aparência de conclusão. A função entrega números com incerteza e deixa a interpretação para quem lê.

### IV.10 Aprendizado de máquina

- `_build_pipeline` — Monta a sequência de processamento do classificador: preencher valores ausentes, padronizar as escalas e ajustar uma regressão logística com penalização. O ponto crucial é que essa sequência inteira é tratada como uma unidade e ajustada apenas com os dados de treino.<br>*Por que isso é essencial:* se você padronizasse as variáveis usando todos os dados antes de dividir em treino e teste, informação do teste vazaria para o treino, e o desempenho apareceria melhor do que realmente é. Esse erro chama-se **vazamento de dados** e é uma das causas mais comuns de resultados que não se replicam.
- `nested_cv_scores` — Executa a validação cruzada aninhada. Vale explicar em duas camadas.
  - *Validação cruzada* é o procedimento de dividir os participantes em grupos, treinar o modelo com todos menos um grupo, testar nesse grupo, e repetir alternando qual grupo fica de fora. Assim cada pessoa é testada por um modelo que nunca a viu.
  - *Aninhada* significa que existe um segundo nível de divisão dentro do treino, usado exclusivamente para escolher os ajustes internos do modelo. Sem esse segundo nível, a escolha dos ajustes olharia indiretamente para os dados de teste, e o desempenho seria otimista. Todo o procedimento é repetido várias vezes com divisões diferentes, para reduzir a dependência de um sorteio particular.
- `metrics_from_probs` — Converte as previsões em métricas: área sob a curva ROC, acurácia balanceada, sensibilidade, especificidade e outras. Os intervalos de confiança incorporam duas fontes de incerteza: a de quais participantes estão na amostra e a de como eles foram divididos entre as partições.<br>Sobre a área sob a curva ROC, já que é a métrica mais citada: ela representa a probabilidade de que, tomando ao acaso um paciente e um controle, o modelo atribua maior probabilidade ao paciente. Vale 0,5 quando o modelo não sabe nada, e 1,0 quando separa perfeitamente. Com 26 contra 28 participantes, um valor verdadeiro de 0,70 vem acompanhado de um intervalo de aproximadamente **0,56 a 0,84** — uma faixa que abriga simultaneamente "quase nenhum sinal" e "clinicamente útil". Reportar o valor pontual sem esse intervalo é enganoso.
- `matched_observed_auc` e `permutation_auc_test` — Juntas, executam o teste de permutação. A lógica do teste é perguntar: se eu embaralhar os rótulos, destruindo qualquer relação verdadeira, com que frequência obtenho um desempenho tão bom quanto o observado? Repete-se o embaralhamento centenas de vezes, constrói-se a distribuição dos desempenhos sob ausência de sinal, e verifica-se onde o valor observado cai nessa distribuição.<br>Existe uma exigência técnica que é fácil violar e cara de errar: o desempenho observado e os desempenhos embaralhados precisam ser medidos exatamente do mesmo jeito, com o mesmo número de repetições. Se o observado usa uma estimativa mais suavizada e os embaralhados usam uma mais ruidosa, as duas distribuições não são comparáveis, o p-valor perde interpretação, e o viés é justamente na direção de parecer significativo. Por isso o programa calcula uma versão do observado deliberadamente pareada ao procedimento nulo, e registra no relatório quantas repetições foram usadas.
- `batch_confound_report` — Mede a associação entre o prefixo do identificador e o diagnóstico. Se o prefixo separa perfeitamente os grupos, marca todos os resultados como confundidos e insere um aviso obrigatório em todos os relatórios.

### IV.11 Orquestração e verificação

- `extract_features_for_window` — Executa, para um participante e um trecho, toda a cadeia: recortar, verificar deriva, fatiar, aplicar controle de qualidade, calcular e agregar as características.
- `run_v24` — A função principal, que coordena tudo na ordem correta. Um detalhe de projeto importante: falhas em um arquivo não interrompem o processamento dos demais. Elas são registradas com o motivo e a etapa em que ocorreram, e o lote prossegue. As falhas são contadas separadamente por etapa, para que um participante que falhou em três trechos diferentes não seja contado como três participantes perdidos.
- `environment_manifest` e `sha256_file` — Registram as versões de todos os programas usados e uma assinatura digital única de cada arquivo de entrada. A assinatura é um código longo que muda completamente se um único caractere do arquivo mudar. Permite provar, meses depois, que os dados analisados são exatamente os mesmos.
- `selftest_v24` — Executa **52 verificações automáticas** com respostas conhecidas de antemão. Alimenta o programa com sinais sintéticos construídos para ter propriedades específicas — um ritmo alfa em exatamente 10 Hz, um transitório inicial de duração conhecida, uma interferência de rede em 50 Hz — e confere se o programa encontra o que deveria encontrar. Também testa as travas: verifica se o modo pesquisa realmente se recusa a abrir sem congelamento, e se as barreiras contra medidas exploratórias e circulares realmente disparam.<br>Todos os 52 testes passam na versão atual. Isso significa que o programa faz o que diz fazer em condições controladas. Não significa que os resultados sobre dados reais estejam corretos — nenhum teste de software pode garantir isso.
- `main_v24` — A interface de linha de comando, que permite executar o programa digitando um comando no terminal.

---

## Parte V — Como executar

Você precisa de Python instalado e de cinco bibliotecas. Em um terminal, digite:

```bash
pip install numpy pandas scipy scikit-learn statsmodels
```
Ou pode rodar:
```bash
pip install -r requirements.txt
```
Passo 1 — Sempre comece pelos autotestes. Se algo estiver errado no ambiente, você descobre aqui e não no meio da análise.
```
python main.py --selftest
```
>  resultado esperado: RESULTADO: 52/52

Passo 2 — Execução em modo integração. Obrigatória, e não consome nada: nenhuma inferência é feita.
```bash
python main.py \
  --eeg-dir  ./modma_eeg_3ch \
  --metadata ./meta.xlsx \
  --output   ./saida_integracao \
  --selection-criterion "55 registros de repouso do experimento de 3 eletrodos do MODMA, lote completo, sem seleção adicional." \
  --mode integration
  ````

O critério de seleção é texto livre e obrigatório. Ele descreve como este lote de arquivos foi escolhido. Sem essa informação é impossível avaliar viés de seleção, e o programa se recusa a rodar.

Passo 3 — Inspecionar e congelar. Abra 02_quality/qc_distribution.csv em qualquer programa de planilha. Para cada métrica, olhe a distribuição observada e a coluna que informa quantos participantes o limiar vigente excluiria. Decida os limiares olhando apenas isso — nunca o diagnóstico. Depois registre a decisão no código:
```bash
cfg = RunConfig(
    qc=QCConfig(max_ocular_index=0.80, max_muscle_ratio=0.60),
    qc_freeze=QCFreezeConfig(
        thresholds_frozen=True,
        freeze_date="2026-08-22",
        freeze_source="qc_distribution.csv da execução de integração de 22/08/2026; "
                      "limiares fixados na distribuição empírica, cegos ao diagnóstico."))
```
Passo 4 — Execução em modo pesquisa, que agora está liberado.

```bash
python main.py --eeg-dir ./data --metadata ./meta.xlsx \
  --output ./saida_pesquisa --selection-criterion "..." --mode research
  ````

## Parte VI — Como os resultados devem ser apresentados
### VI.1 A ordem de apresentação
A ordem em que os resultados aparecem comunica uma mensagem, queira você ou não. Se a acurácia vem primeiro, o leitor entende que o estudo é sobre desempenho. Por isso a ordem recomendada aqui é deliberadamente invertida:

    1. Fluxograma de participantes — de 55 arquivos até a coorte final, com cada exclusão explicada.
    2. Achados estruturais — confundimento de lote, redundância entre canais, deriva temporal, análise de poder.
    3. Qualidade e protocolo — distribuição das métricas, limiares congelados, verificação da taxa de amostragem, janela derivada.
    4. Comparações entre grupos — as nove medidas, com tamanhos de efeito e intervalos.
    5. Classificação — desempenho, comparadores, teste de permutação.
    6. Robustez temporal — concordância entre trechos, com incerteza.

Colocar os achados estruturais antes dos resultados de desempenho não é modéstia retórica. É a única ordem em que os números de desempenho podem ser lidos corretamente.

### VI.2 Como apresentar cada tipo de resultado
|Resultado|Forma correta de apresentar|Erro comum a evitar|
|:---:|:---:|:---:|
|Fluxograma|Diagrama com todas as etapas e contagens, mais o teste de exclusão diferencial entre grupos.|Reportar apenas o número final da coorte.|
|Confundimento de lote|Tabela cruzada completa e o valor de V de Cramér, no resumo do artigo.|Mencionar apenas na seção de limitações.|
|Poder estatístico|Tabela com o efeito detectável, conduzida antes dos resultados.|Calcular poder depois, a partir do efeito observado — prática sem validade.|
|Tamanhos de efeito|Sempre com intervalo de confiança, em unidades padronizadas.|Reportar só o p-valor.|
|Área sob a curva ROC|Com o intervalo que inclui a variabilidade da validação cruzada.|Valor pontual isolado, ou intervalo que ignora as partições.|
|Comparadores|Lado a lado: classe majoritária, só demografia, só EEG, EEG mais demografia.|Reportar apenas o modelo de EEG.|
|Teste de permutação|Com o número de repetições explicitado e nota de que não deve ser comparado com a estimativa suavizada.|Misturar as duas estimativas na mesma frase.|
|Robustez|Proporção de concordância com intervalo e p binomial.|Afirmar "resultados robustos" sem quantificar.|

## VI.3 Como discutir

A discussão deve responder três perguntas, nesta ordem.

**Primeira:** o que estes dados podem responder? Eles permitem caracterizar propriedades técnicas do sinal frontal em repouso, quantificar deriva temporal, estimar tempos de estabilização e demonstrar um procedimento de análise auditável. Permitem, sobretudo, mostrar como a estrutura de um conjunto de dados pode inviabilizar a pergunta clínica sem que isso apareça em nenhuma métrica de desempenho.

**Segunda:** o que eles não podem responder? Não permitem afirmar que qualquer diferença encontrada seja de origem cerebral. Com o lote separando perfeitamente os grupos, a hipótese "os grupos foram gravados em condições diferentes" explica os dados exatamente tão bem quanto a hipótese "os grupos têm eletrofisiologia diferente". Nenhum dado deste conjunto pode arbitrar entre as duas.

**Terceira:** o que isso implica para a literatura? Se a colinearidade entre identificador e diagnóstico é propriedade do conjunto e não do nosso recorte, então trabalhos publicados sobre esses dados que não a reportem estão sujeitos à mesma ambiguidade. Uma recomendação concreta: que estudos futuros publiquem a tabela cruzada entre identificador e grupo, e que revisores a solicitem.

## VI.4 Que conclusões podem e não podem ser tiradas

|Se o resultado for...|Conclusão legítima|Conclusão ilegítima|
|:---:|:---:|:---:|
|Nenhuma medida significativa após correção|Não detectamos diferenças, com poder limitado a efeitos grandes. Ausência de evidência não é evidência de ausência.|"O EEG frontal não distingue depressão."|
|Uma ou duas medidas significativas|Um sinal candidato, com magnitude provavelmente inflada e indistinguível de efeito de lote. Requer replicação independente.|"Identificamos um biomarcador."|
|Área sob a curva alta|Um classificador que separa os grupos — que também são lotes distintos. Compatível com um excelente detector de lote.|"O modelo detecta depressão com X% de acurácia."|
|Modelo demográfico igual ou melhor|O EEG não acrescenta informação além da demografia nesta amostra.|Ignorar o comparador e reportar só o EEG.|
|Concordância entre trechos|Estabilidade interna ao mesmo conjunto de participantes.|Tratar como validação externa.|
|Permutação não significativa|Não há evidência de sinal além do acaso, com o poder disponível.|"Provamos que não há sinal."|

## VI.5 Dez frases que não devem ser escritas
    • "O EEG discrimina depressão de controles" — o lote também discrimina, e perfeitamente.

    • Comparar o p-valor da permutação com a estimativa suavizada de desempenho.

    • Citar qualquer frequência ou tempo sem antes checar a verificação da taxa de amostragem.

    • Chamar de "derivada dos dados" uma janela que foi truncada por limite administrativo.

    • Reportar o intervalo da área sob a curva ignorando a variabilidade das partições.

    • Interpretar concordância entre trechos como validação externa.

    • Apresentar os modelos ajustados por idade, sexo e escolaridade como se atenuassem o efeito de lote.
    
    • Tratar exclusões de qualidade como neutras sem reportar o teste de exclusão diferencial.

    • Calcular poder estatístico depois, usando o efeito observado.

    • Descrever a assimetria frontal alfa como achado, quando ela é exploratória e medida em sítio não validado.

## VI.6 Checklist antes de submeter
- [ ] Os 52 autotestes passaram nesta máquina, com as versões registradas no manifesto.
- [ ] A execução em modo integração foi concluída e a distribuição de qualidade foi inspecionada.
- [ ] Os limiares foram congelados, com data e justificativa registradas.
- [ ] A verificação da taxa de amostragem foi revisada e seu estado consta do texto.
- [ ] O protocolo temporal é válido, ou a limitação está declarada explicitamente.
- [ ] O fluxograma está completo, com o teste de exclusão diferencial.
- [ ] A permutação usou estimadores pareados, e o número de repetições consta do relatório.
- [ ] O intervalo da área sob a curva inclui a variabilidade da validação cruzada.
- [ ] Confundimento de lote e poder estatístico aparecem no resumo, não apenas nas limitações.
- [ ] O manifesto foi arquivado com a assinatura digital de todos os arquivos de entrada.

## Parte VII — As três limitações que nenhuma linha de código resolve

Estas limitações são do desenho do estudo, não do programa. Nenhuma melhoria de software as corrige. Elas precisam constar do resumo de qualquer publicação, não apenas do rodapé.

|Limitação|Número verificado|Consequência|
|:---:|:---:|:---:|
|Confundimento de lote| V de Cramér = 1,00. Prefixo 0201 corresponde a 26 de 26 pacientes; prefixos 0202 e 0203 a 29 de 29 controles.| Efeito de lote e efeito clínico são matematicamente inseparáveis. Ajustar por idade, sexo e escolaridade não toca neste problema.|
|Poder insuficiente|Efeito detectável com 80% de poder: g = 0,78 (1,02 com correção). Poder de 44% para g = 0,50 e de 19% para g = 0,30.|Efeitos de magnitude realista são indetectáveis. Achados significativos virão com magnitude inflada.|
|Redundância espacial|Correlação mediana entre canais = 0,983; 38 de 55 participantes acima de 0,95.|Os três eletrodos carregam essencialmente um sinal. Medidas de assimetria não são sustentáveis.|

Há ainda quatro limitações menores, declaradas em todo relatório: a taxa de amostragem carece de corroboração documental plena; o fator de conversão para microvolts é inferido, não calibrado; com três eletrodos frontais não é possível separar artefatos por decomposição em componentes; e nenhum conjunto de validação externa foi executado — o código está preparado para isso, mas preparação não é validação.

## Parte VIII — Glossário

|Termo|Significado neste estudo|
|:---:|:---:|
|Acurácia balanceada|Média entre a taxa de acerto nos pacientes e nos controles. Preferível à acurácia simples quando os grupos têm tamanhos diferentes.|
|Adimensional|Medida que é uma proporção ou razão, sem unidade. Não depende de calibração de amplitude.|
|Área sob a curva ROC|Probabilidade de o modelo dar nota maior a um paciente do que a um controle, tomados ao acaso. 0,5 é o acaso; 1,0 é a separação perfeita.|
|Artefato|Qualquer coisa registrada pelo eletrodo que não venha do cérebro: piscadas, músculo, movimento do cabo, rede elétrica.|
|Bootstrap| Método de estimar incerteza por reamostragem repetida dos próprios dados.
|Característica (feature)|Um número que resume alguma propriedade do sinal e serve de entrada para a análise.|
|Cego ao rótulo|Cálculo feito sem acesso ao diagnóstico. Impede que a escolha de parâmetros favoreça um resultado.|
|Coorte analítica|O grupo final de participantes que efetivamente entra na análise, após todas as exclusões.|
|Confundimento|Situação em que duas causas possíveis variam juntas e não podem ser separadas pelos dados.|
|Congelamento de limiares|Ato datado de fixar os critérios de qualidade após ver sua distribuição e antes de ver os resultados.|
|Deriva|Mudança lenta e sistemática do sinal ao longo do tempo.|
|Desvio absoluto mediano|Medida de dispersão resistente a valores extremos, usada para detectar atípicos.|
|Entropia|Medida de imprevisibilidade. Alta significa muito espalhado; baixa, muito concentrado.|
|Época (ou fatia)|Um segmento curto do registro, aqui de quatro segundos, tratado como unidade de cálculo.|
|Erro tipo M|Superestimação da magnitude de um efeito, consequência de estudos com poder baixo.|
|Espectro|Descrição de um sinal em termos das frequências que o compõem.|
|Filtro passa-banda|Filtro que preserva apenas as frequências dentro de uma faixa.|
|Fluxograma de participantes|Diagrama que mostra quantos participantes entraram, quantos saíram em cada etapa e por quê.|
|Hedges g|Tamanho do efeito padronizado, corrigido para amostras pequenas.|
|Hertz (Hz)|Ciclos por segundo. Um ritmo de 10 Hz oscila dez vezes por segundo.|
|Impressão digital criptográfica|Código curto derivado de um conteúdo, que muda completamente se o conteúdo mudar.|
|Intervalo de confiança|Faixa de valores compatíveis com os dados, dado um nível de confiança.|
|Janela de análise|O trecho do registro efetivamente usado para calcular as características.|
|Medida confirmatória|Uma das nove medidas pré-declaradas que entram nos testes formais.|
|Medida exploratória|Medida descritiva que gera hipóteses, mas nunca entra no teste confirmatório.|
|Notch|Filtro estreito que remove uma frequência específica, aqui a da rede elétrica.|
|P-valor|Probabilidade de observar um resultado tão ou mais extremo, se não houvesse efeito real.|
|Poder estatístico|Probabilidade de detectar um efeito que existe de verdade.|
|Potência relativa|Fração da energia total do sinal que está em uma faixa de frequência.|
|Regressão logística|Modelo estatístico que estima a probabilidade de pertencer a um grupo.|
|Reprodutibilidade|Propriedade de uma análise poder ser repetida por outra pessoa com o mesmo resultado.|
|Semente aleatória|Número fixo que torna sorteios do computador repetíveis.|
|Tamanho do efeito|Magnitude da diferença, expressa em desvios padrão. Não depende do tamanho da amostra.|
|Taxa de amostragem|Quantas medidas por segundo o aparelho registrou.|
|Taxa de falsas descobertas|Proporção esperada de falsos positivos entre os resultados declarados significativos.|
|Tempo de estabilização|Instante a partir do qual o registro deixa de mudar sistematicamente.|
|Validação cruzada|Procedimento de testar o modelo em dados que ele não viu durante o treino.|
|Validação cruzada aninhada|Versão com dois níveis: um para estimar desempenho, outro para escolher os ajustes do modelo.|
|Validação externa|Testar o modelo em um conjunto de dados totalmente independente. Não foi feita aqui.|
|Vazamento de dados|Erro em que informação do teste influencia o treino, inflando o desempenho.|
|V de Cramér|Medida de associação entre duas variáveis categóricas, de 0 a 1.|
