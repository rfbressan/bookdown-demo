# Métodos de Calibração {#calibracao} 



Neste post iremos mostrar as diferenças existentes entre "interpolação", "suavização" e "parametrização" de superfícies de volatilidade.

Como já apresentado em posts anteriores, existem diversas formas de interpolar, extrapolar, parametrizar e calibrar smiles de volatilidade. Exsitem vantagens e desvantagens para cada método. 

## Calibração versus Interpolação

Uma forma simples de gerar um smile de volatilidade a partir de dados observados no mercado é a **interpolação** destes dados. Diversas formas de interpolação existem, sendo talvez a mais conhecida a spline cúbica. Não é a proposta deste artigo detalhar os procedimentos de interpolação, restando saber que em tal procedimento é gerada uma função contínua em partes (piecewise) que **passa** por todos os pontos observados. 

Uma interpolação força a passagem da função interpolada em todos os seus pontos de referência, como se estivesse ligando os pontos em um desenho a mão livre. Portanto, nestes pontos o erro da interpolação é zero por definição, entretanto em pontos intermediários podem surgir erros, inclusive aqueles que possibilitam a existência de arbitragem entre strikes de um mesmo smile[^1].

Em contraposição a métodos de interpolação, podem ser derivados métodos de suavização (smoothing) ou então a parametrização do smile de volatilidade. Seja suavização, ou parametrização, estes métodos não forçam a passagem da função que representa o smile pelos pontos de mercado, mas buscam minimizar alguma função perda com relação aos desvios em relação a estes pontos ao mesmo tempo em que buscam "suavizar" o smile, para que este não apresente variações bruscas entre os strikes ou alterações de convexidade na curva de preços, que não são condizentes com a teoria de precificação de derivativos.

Um método paramétrico, como o SVI, Heston, SABR ou Volatilidade Local, busca ajustar às volatilidades implícitas observadas através dos preços das opções sendo praticados no mercado uma determinada função, que possui parâmetros em sua definição que por sua vez determinam a forma desta função. Ao se ajustar os parâmetros, pode-se adequar a função para ficar "o mais próxima possível" dos dados observados, sem necessariamente, no entanto, passar por todos estes pontos.

A figura abaixo tenta mostrar as diferenças entre uma interpolação spline cúbica, uma suavização e uma parametrização SVI. Enquanto que a interpolação liga todos os pontos marcados, a suavização e a parametrização não necessariamente passam sobre estes pontos mas fornecem uma curva mais "suave", sem trocas de convexidade, o que gera oportunidades de arbitragem e probabilidades negativas de ocorrência de determinados preços para o ativo subjacente, que ferem os princípios de precificação de opções. Os dados utilizados neste e nos próximos artigos sobre superfícies de volatlidade foram obtidos do site [ivolatility.com](http://www.ivolatility.com/doc/usa/IV_Raw_Delta_surface.csv) na forma de amostra gratuita fornecida livremente. O ativo subjacente é o ETF [`IWM`](https://www.ishares.com/us/products/239710/ishares-russell-2000-etf) para a data de 21/09/2017. 

<div class="figure">
<img src="06-metodos-calibracao_files/figure-epub3/diferencas-1.png" alt="Diferentes métodos de ajuste de dados a um smile."  />
<p class="caption">(\#fig:diferencas)Diferentes métodos de ajuste de dados a um smile.</p>
</div>

Pode-se verificar como os métodos SVI e a suavização não passam sobre todos os pontos marcados, com a suavização tendo dificuldade com a curvatura nos valores mais altos de _moneyness_ e a SVI possuindo uma inclinação mais branda na asa esquerda do smile.

## Spline cúbica

Este método é possivelmente um dos mais flexíveis e conhecidos de interpolação de dados univariados existente, embora também exista sua versão bi-dimensional. Uma spline nada mais é que "uma curva definida matematicamente por dois ou mais pontos de controle"[^2]. 

No caso da spline cúbica, esta é uma função polinomial de grau 3 definida em cada subintervalo demarcados pelos pontos de controle, no caso de interpolação são todos nós. Ou seja, considere um segmento entre dois pontos consecutivos $[c, d]\in S$ a spline é uma função cúbica com seus parâmetros calculados pelo algoritmo de ajuste. Para o próximo intervalo de pontos dentro do domínio da função, um novo polinômio de grau 3 é ajustado, sendo que nos pontos de nós uma restrição de igualdade entre as derivadas nos dois segmentos é aplicada para garantir a suavidade da função interpolada como um todo.

Assim, uma spline cúbica é uma função contínua, suave e diferenciável até a segunda ordem. Entretanto, suas derivadas, apesar de contínuas, podem não ser suaves, especialmente aquela de segunda ordem que pode apresentar pontos de "ruptura". Esta característica de uma spline cúbica a torna pouco atrativa para a inferência de distribuições de probabilidade a partir de dados de volatilidade ou mesmo dos preços de opções.

<div class="figure">
<img src="./images/cubic_spline.png" alt="Cada segmento de uma spline cúbica é um polinômio de grau 3 diferente."  />
<p class="caption">(\#fig:cubic-spline)Cada segmento de uma spline cúbica é um polinômio de grau 3 diferente.</p>
</div>


## Suavização

A técnica de suavização é muito semelhante a interpolação, inclusive o método spline também é aplicado, com algumas modificações de forma que nem todos os pontos fornecidos serão nós. 

Na spline de suavização (ou aproximação), os pontos fornecidos são separados entre os nós, onde a função deve passar e pontos de controle, que são utilizados para controlar a curvatura da função nestes pontos.

Estas suavizações são principalmente utilizadas quando se possui muitas observações sujeitas a ruídos, de forma que uma interpolação entre todos os pontos seria tanto impraticável quanto sem sentido. O que se deseja, portanto, é uma função **aproximada** que melhor descreva o processo sob análise.

Um ponto em comum entre estas técnicas é o parâmetro de suavização, ausente, na interpolação, que controla a "suavidade" da função estimada.

<div class="figure">
<img src="06-metodos-calibracao_files/figure-epub3/suavizacao-1.png" alt="Menor parâmetro de suavização gera granularidade na curva."  />
<p class="caption">(\#fig:suavizacao)Menor parâmetro de suavização gera granularidade na curva.</p>
</div>

## Parametrização

E por fim as técnicas de parametrização. Nesta categoria estão diversos conhecidos modelos de superfícies de volatilidade implícita, dentre eles os modelos de @Heston1993, Volatilidade Local de @Dupire1994 e SVI de @Gatheral2004.

Em comum, estes modelos tentam parametrizar a superfície, e por conseguinte o smile de volatilidade, de acordo com alguma função, em geral não-linear, que possui características condizentes com a teoria de precificão de derivativos e também a observação empírica das superfícies.

Por exemplo, a parametrização _raw_ da SVI possui a seguinte forma para a **variância total**[^3] :

$$ w(k) = a + b\left(\rho(k-m)+\sqrt{(k-m)^2 + \sigma^2}\right)$$

que fornece um espaço de cinco parâmetros $\chi_B=\{a, b, \rho, m, \sigma\}$ que definem o smile e devem, portanto, serem calibrados a partir de dados observados no mercado.

O procedimento de calibração consiste em encontrar o conjunto de parâmetros que minimizam uma função perda entre a volatilidade prevista pelo modelo e os dados de mercado, enquanto satisfazem algumas restrições adicionais, como "ausência de arbitragem", suavidade, etc. Trata-se, via de regra, de problemas de otimização não-linear com restrições de inequalidade também não-lineares.

### Função perda

A função perda, ou função de calibração pode ser definida de diversas maneiras, de forma geral, para uma determinada maturidade, ela toma a forma:

$$L=\sum\limits_{i=1}^n\lambda_i||\hat w(k_i)-w_{imp}(k_i)||$$
onde $||\cdot||$ é alguma medida de norma, sendo a mais conhecida o quadrado das diferenças, dando origem a minimização do erro quadrático médio (RMSE). Para este smile sendo calibrado existem $n$ strikes ($k_i$) e suas volatilidades implícitas observadas são $w_{imp}(k_i)$. A resposta do modelo para um determinado strike é $\hat w(k_i)$ e $\lambda_i$ são os pesos dados na função perda para cada um destes strikes.

Os pesos $\lambda_i$ são utilizados para ponderar as observações das volatilidades mais importantes para o cálculo, onde se deseja que a curva ajustada possua  um menor erro. Em geral, estes pesos são calculado como inversamente proporcionais: 

  - ao quadrado dos _spreads bid-ask_, para dar mais importância às opções mais líquidas
  - ao quadrado da grega vega calculada a partir do modelo BSM

### Otimizadores

Os otimizadores são os algoritmos pelos quais o problema de minimização posto é resolvido. Se a função perda é convexa, e ela deve ser construída de forma a ser, mesmo que não estritamente, então ela possui um ou mais pontos de mínimo onde o gradiente desta função é igual a zero. O que os otimizadores fazem é buscar o conjunto de parâmetros que minimizam a função perda e atendem as restrições impostas simultaneamente. Os otimizadores podem ser classificados em dois grandes grupos, globais e locais. 

Algoritmos locais dependem de uma estimativa inicial dos parâmetros para começarem a busca pelo mínimo. Seguindo uma regra utilizando-se o gradiente da função ou alguma heurística, estes otimizadores caminham em direção ao ponto de mínimo mais próximo da estimativa inicial, daí o nome "local". Como desvantagem destes otimizadores é a mais evidente é que se a função perda for altamente não-linear, com diversos pontos de mínimo local, este otimizador pode ficar preso em um destes pontos sem nunca, no entanto, encontrar o mínimo global. Eles são, portanto muito sensíveis à estimativa inicial dos parâmetros.

Por sua vez, otimizadores globais buscam mapear todo o espaço factível para os parâmetros e encontrar o ponto mínimo da função perda dentro deste espaço. Estes algoritmos não dependem de estimativas iniciais, uma vez que tentarão avaliar o espaço completo. São utilizados quando o problema de minimização é não-linear e possui múltiplos pontos de mínimo local. Estes algoritmos usam alguma forma de heurística para encontrar a região onde o mínimo global está localizado, mas são, em geral, ineficientes em apontar rapidamente onde este ponto de mínimo se encontra com precisão. Por esta razão, é frequente a utilização de otimizadores globais com um posterior refinamento de sua solução por algum algoritmo local.

Abaixo apresentamos alguns exemplos mais comuns de otimizadores, tanto locais quanto globais:

- **Gauss-Newton**: Este método é utilizado para encontrar as raízes de alguma função. Para encontrar o ponto de mínimo da função perda, precisa-se encontrar as raízes do gradiente desta função, portanto o método de Newton em otimização faz uso da função gradiente. Este é um método de otimização local. 

- **Levenberg-Marquardt**: Método muito utilizado para problemas não-lineares, ele parte de uma modificação ao método de Gauss-Newton ao introduzir um fator de amortecimento calculado iterativamente.

- **L-BFGS-B**: BFGS é um método conhecido como quasi-Newton, onde não é necessário calcular a Hessiana do problema, ela é aproximada a partir do próprio gradiente. É bastante utilizado para resolver problemas não-lineares e em sua versão L-BFGS-B pode lidar com restrições do tipo _box_, intervalo dos parâmetros é fixo. 

- **Nelder-Mead**: Este é um método livre do uso de gradiente, já que usa uma heurística para construir um simplex e a partir deste "procurar" por um mínimo. Bastante utilizado quando a função objetivo pode não ser diferenciável. Faz uso de um simplex inicial, que pode ser grande o suficiente para encampar o mínimo global, entretanto, não se classifica como um otimizador global.

- **Algoritmo Genético**: Este método utiliza conceitos da seleção natural para gerar os resultados da otimização. É um otimizador global, no sentido que independe de uma estimativa inicial de parâmetros e faz uma busca por todo o espaço factível. Em um algoritmo genético, uma população aleatória inicial de parâmetros é criada e a partir desta, as gerações evoluem conforme mutações e _cross-over_ de características e é avaliado o _fitness_ de cada conjunto de parâmetros até um deles ser considerado adequado.

- **Evolução Diferencial**: É um método de otimização global, assim como o Algoritmo Genético e o Enxame de Partículas. Sua diferença reside no fato de que sua população inicial é constantemente avaliada e deslocada de posição. Se o agente obtiver uma situação melhor (menor valor para a função perda) na nova posição, esta agora faz parte da população. Desta forma os agentes, antes espalhados pelo espaço factível dos parâmetros, tendem a convergir para um ponto com o menor valor da função perda.

- **Enxame de Partículas**: Do inglês, _Particle Swarm Optimization - PSO_ este método é semelhante ao DE _(Differential Evolution)_ porém as partículas (o equivalente dos agentes no DE) matém  informações sobre a posição da melhor partícula até então, de forma a fazer com que as partículas tendam para a melhor solução.

## Conclusão

Dependendo do objetivo da aplicação, superfícies de volatilidade podem ser interpoladas, suavizadas ou parametrizadas. A parametrização tem recebido especial interesse pois pode, ao mesmo tempo que garante uma superfície livre de arbitragem estática se devidamente construída, ajustar-se muito bem aos dados observados e gerar distribuições neutras ao risco implícitas factíveis.

Para gerar uma superfície parametrizada, primeiramente é necessário um modelo teórico com propriedades desejáveis e que se ajuste aos dados de mercado quando calibrado. Escolhido este modelo paramétrico, passa-se a calibração do mesmo onde exsitem diversas opções de escolha entre otimizadores. Ao final do processo teremos um modelo de superfície devidamente parametrizado com valores que melhor se ajustam segundo alguma função perda escolhida.

Com a superfície de volatilidade calibrada, as aplicações possíveis incluem a precificação de derivativos, gerenciamento de risco, simulações de Monte Carlo, análises de stress, entre outras.

## Referências

[^1]: Veja mais detalhes no artigo anterior: [Smile de volatilidade - parte 2](/2019/01/25/smile-de-volatilidade-parte-2)
[^2]: Definição retirada de https://pt.wikipedia.org/wiki/Spline
[^3]: A variância total é definida pelo produto entre a variância implícita e o tempo para expiração, $w=\sigma^2_{imp}\cdot\tau$.
