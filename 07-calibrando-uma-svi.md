# Calibrando uma SVI {#svi}



Neste capítulo iremos mostrar como fazer uma calibração de um smile SVI baseado nos trabalhos de [@Gatheral2004] e [@DeMarco2009]. Escolheremos apenas uma fatia da superfície de volatilidade, fixando o tempo para expiração (maturidade) e coletando as volatilidades implícitas para diversos strikes.

Como já apresentado em posts anteriores, existem diversas formas de interpolar, extrapolar, parametrizar e calibrar smiles de volatilidade. Exsitem vantagens e desvantagens para cada método. Neste post iremos fixar nossa atenção no modelo paramétrico de smile chamado SVI - Stochastic Volatility Inspired - uma forma que une a "simplicidade" de um modelo paramétrico com o poder de adesão aos dados de mercado dos modelos de volatilidade estocástica (i.e. Heston, SABR e afins).

Chamar o modelo SVI de simples é puro eufemismo, ele é um modelo poderoso, com fundamento teórico avançado e diversos detalhes para sua calibração.

## Modelo SVI

Este modelo foi apresentado por [Jim Gatheral](https://mfe.baruch.cuny.edu/jgatheral/) na conferência _Global Derivatives & Risk Management 2004_ e foi bem recebido pelos profissionais de mercado interessados em superfícies de volatilidade para _equities_, principalmente.

O modelo possui duas propriedades que são as razões para sua popularidade. Ele satisfaz a fórmula do momento de @Lee2004, que é um resultado independente de modelo que especifica os limites assintóticos para um _smile_ de volatilidade implícita. Portanto, o modelo SVI é válido para extrapolação além da região central dos dados disponíveis. Além disso, afirma-se que o modelo SVI é relativamente fácil de calibrar para dados de mercado, de modo que a superfície de volatilidade implícita correspondente é livre de arbitragem de calendário. As condições que garantem a ausência de arbitragem de borboleta forma resolvidas em um segundo artigo por @Gatheral2014.

No SVI é possível se estabelecer condições explícitas em seus parâmetros, de modo que o modelo não gere preços onde oportunidades de [arbitragem estáticas](#arbestatica) possam ocorrer. A calibração para dados reais de mercado requer algoritmos de otimização não-linear e pode ser bastante demorada. Mais recentemente, um método para calibração que usa a estrutura inerente do modelo para reduzir as dimensões do problema de otimização foi desenvolvido em @DeMarco2009.

A parametrização conhecida como _RAW_ do SVI é apresentada na equação \@ref(eq:rawsvi), seguindo a notação já introduzida [anteriormente](#smile), portanto, estamos modelando a **variância total** implícita para um determinado prazo. Para diferentes maturidades, teremos diferentes conjuntos de parâmetros.

\begin{equation}
w(k) = a + b\left(\rho(k-m)+\sqrt{(k-m)^2 + \sigma^2}\right)
(\#eq:rawsvi)
\end{equation}

onde: $a \in \mathbb R$, $b \geq 0$, $|\rho| < 1$, $m \in \mathbb R$, $\sigma > 0$, e $a+b \sigma\sqrt{1 − \rho^2} \geq 0$ para garantir que $\min w(k)>0, \, \forall k \in \mathbb R$.

## Restrições de não-arbitragem

Antes de demonstrar a restrição imposta aos parâmetros $b$ e $\rho$ em função dos limites de inclinação das asas do _smile_, vamos derivar as expressões para $w\prime(k)$ e $w\prime\prime(k)$ que nos serão úteis na demonstração.

A expressão para $w\prime(k)$ é bastante simples:

\begin{equation}
w\prime(k) = b \left[\rho + \frac{(k-m)}{\sqrt{(k-m)^2+\sigma^2}}\right]
(\#eq:wk)
\end{equation}

Derivando novamente a equação \@ref(eq:wk) em relação a $k$ teremos uma expressão ainda mais simples, mesmo que após alguma manipulação algébrica um tanto tediosa[^71], e resulta em:

\begin{equation}
w\prime\prime(k)=\frac{b\sigma^2}{[(k-m)^2+\sigma^2]^{3/2}}
(\#eq:wkk)
\end{equation}

onde, se considerarmos $b>0$ temos que $w\prime\prime(k)>0, \,\forall k\in \mathbb R$, ou seja, o _smile_ de volatilidade definido pela equação \@ref(eq:rawsvi) é **estritamente** convexo.

@Rogers2010 definiram os limites possíveis para a inclinação das asas em função do tempo para expiração, provando que o _smile_ tende a ficar mais horizontal a medida que o prazo aumenta. Este limite pode ser escrito da seguinte forma e é uma **condição necessária** para a ausência de arbitragem:

\begin{equation}
|w\prime(k)|\leq \frac{4}{\tau} \qquad \forall k \in \mathbb R, \quad \forall \tau \in (0, \infty)
(\#eq:rogers)
\end{equation}

Sendo o _smile_ convexo, suas máximas inclinações ocorrem quando $k\rightarrow \pm \infty$. Portanto, deve-se avaliar a restrição dada pela equação \@ref(eq:rogers) nestes limites da seguinte maneira:

\begin{align}
\lim\limits_{k\rightarrow\infty}w\prime(k)&=b(1+\rho)\geq 0\\
\lim\limits_{k\rightarrow-\infty}w\prime(k)&=-b(1-\rho)\leq 0
\end{align}

que satisfazendo estas duas relações ao mesmo tempo em que se restringe os parâmetros $b$ e $\rho$ através da inequalidade de Rogers e Tehranchi nos garante o seguinte resultado para um SVI **livre de arbitragem de travas**.

\begin{equation}
b(1+|\rho|)\leq\frac{4}{\tau}
(\#eq:trava)
\end{equation}

Para garantir que a superfície gerada está livre de arbitragem do tipo borboleta deve-se primeiramente definir uma função[^72] $g: \mathbb R\rightarrow \mathbb R$, tal que:

\begin{equation}
g(k)=\left(1-\frac{kw\prime(k)}{2w(k)}\right)^2-\frac{w\prime(k)^2}{4}\left(\frac{1}{w(k)}+\frac{1}{4}\right)+\frac{w\prime\prime(k)}{2}
(\#eq:g)
\end{equation}

e seguir o lema :

**Lema 1** Uma fatia da superfície de volatilidade está livre de arbitragem do tipo borboleta se, e somente se, $g(k) \geq 0$ para todo $k \in \mathbb R$ e $\lim\limits_{k\rightarrow+\infty}d_1(k)=-\infty$.

Infelizmente, a natureza altamente não-linear da função $g(k)$ impossibilita a derivação de restrições gerais aos parâmetros do SVI. A forma mais simples de eliminar arbitragem do tipo borbobleta é incluir a restrição $g(k) \geq 0$ na função perda e proceder com a calibração dos parâmetros.

## Reparametrização Quasi-explicit

Um dos problemas mais marcantes com a calibração do SVI dado pela equação \@ref(eq:rawsvi) é sua natureza altamente não-linear que gera inúmeros pontos de mínimo locais. Mesmo em um ambiente simulado, um típico otimizador de [mínimos quadrados](#calibracao) como Levenberg-Marquardt não consegue chegar ao mínimo global, onde a função perda é igual a zero. A solução encontrada é dependente dos valores iniciais inputados ao otimizador e a robustez do conjunto de parâmetros encontrados não é garantida.

Uma forma de contornar este problema pode ser a utilização de otimizadores globais, como algortimos genéticos, em um primeiro estágio e então o refinamento desta solução através de um otimizador local (LM, por exemplo).

Outra forma, adotada em @DeMarco2009 é a reparametrização da equação \@ref(eq:rawsvi) de forma que esta possa ser tratada como um prolema linear. Para tanto, considere a seguinte troca de variáveis:

\begin{equation}
y = \frac{k-m}{\sigma}
\end{equation}

então a parametrização _RAW_ do SVI se torna.

\begin{equation}
w(y) = a + b\sigma\left(\rho y + \sqrt{y^2 + 1}\right)
\end{equation}

Definindo agora as seguintes variáveis reparametrizadas é possível reduzir a dimensão de parâmetros de um SVI de 5 para apenas 3:

\begin{align}
c = &b\sigma\\
d = &\rho b \sigma
\end{align}

\begin{equation}
w(y)=a+dy+c\sqrt{y^2+1}
(\#eq:quasiexplicit)
\end{equation}

Portanto, para um par fixo de $(m, \sigma)$ nosso problema reduzido é:

\begin{equation}
P_{m, \sigma}:=\min\limits_{a, c, d \in D}f_y(a, c, d)
(\#eq:reduzido)
\end{equation}

onde $f_y(\cdot)$ é a função objetivo da reparametrização, e é dada pela seguinte equação:

\begin{equation}
f_y(a, c, d)=\sum_{i=1}^{n}\left[w(y_i)-\tilde w_i\right]^2
(\#eq:fy)
\end{equation}

onde $\tilde w_i$ é a variância total observada correspondente ao _moneyness_ $k_i$.

O domínio $D$ dos parâmetros $\{a, c, d\}$ é encontrado a partir do limite imposto por \@ref(eq:trava).

\begin{equation}
	D =  
	\begin{cases}
  	0 \leq c \leq 4\sigma\\
  	|d| \leq c \quad \text{e}\quad |d| \leq 4\sigma - c\\
  	0 \leq a \leq \max\{\tilde w_i\}\\
	\end{cases}
	(\#eq:D)
\end{equation}

O problema reduzido dado pela equação \@ref(eq:reduzido), é um típico problema de **mínimos quadrados** com restrições lineares. Este problema, por ser convexo, admite uma única solução interior (se existente) que será o mínimo global para este problema e é encontrada através do gradiente da função objetivo igualando-o a zero, $\nabla f_y = 0$. Esta equação gera um sistema linear nos parâmetros $a, c, d$ que pode ser explicitamente resolvido. Caso a solução encontrada para este problema esteja contida no domínio $D$, esta solução é interior e é o mínimo desejado, caso contrário, deve-se percorrer o perímetro do domínio e encontrar o menor valor da função objetivo que será uma solução de canto.

Seja $(a^*, c^*, d^*)$ a solução de \@ref(eq:reduzido) e $(a^*, b^*, \rho^*)$ os correspondentes parâmetros originais recuperados, então o problema completo de calibração é:

\begin{equation}
P:=\min\limits_{m, \sigma}\sum_{i=1}^n (w_*(k_i)-\tilde w_i)^2
(\#eq:completo)
\end{equation}

onde $w_*(k)=a^*+b^*\left(\rho^*(k-m)+\sqrt{(k-m)^2 + \sigma^2}\right)$.

O problema completo, \@ref(eq:completo) é um problema em apenas duas dimensões, $(m, \sigma)$ e não-linear, que deve ser abordado através de algum tipo de otimizador global.

### Solução explícita do problema reduzido {#pq}

Nesta seção apresentaremos a solução, sem considerar as restrições impostas em \@ref(eq:D) para o problema reduzido em \@ref(eq:reduzido), algo omitido em @DeMarco2009. Esta seção é opcional para o leitor atento que já percebeu a semelhança entre o problema reduzido e um típico problema de regressão linear múltipla.

Para encontrar o conjunto de parâmetros $(a^*, c^*, d^*)$ que representam os valores ótimos na equação \@ref(eq:reduzido), devemos resolver o seguinte sistema de equações:

\begin{equation}
\nabla f_y = \left[
	\begin{array}{c}
		\partial f_y / \partial a\\ 
		\partial f_y / \partial d\\
		\partial f_y / \partial c
	\end{array} 
\right] = \boldsymbol{0}
(\#eq:gradiente)
\end{equation}

Cada uma das derivadas parciais da equação acima quando igualadas a zero dão origem ao sistema linear apresentado abaixo:


\begin{equation}
\scriptsize
\begin{bmatrix}
&n &\sum y_i &\sum\sqrt{y_i^2+1}\\
&\sum y_i &\sum y_i^2 &\sum(y_i\sqrt{y_i^2+1})\\
&\sum\sqrt{y_i^2+1} &\sum(y_i\sqrt{y_i^2+1}) &\sum(y_i^2+1)\\
\end{bmatrix}
\cdot
\begin{bmatrix}
a \\
d \\
c
\end{bmatrix}
=
\begin{bmatrix}
\sum\tilde w_i \\
\sum \tilde w_i y_i \\
\sum(\tilde w_i\sqrt{y_i^2+1})
\end{bmatrix}
(\#eq:linear)
\end{equation}

Portanto, o problema reduzido pode ser resolvido através de um sistema linear de ordem 3 sob restrições também lineares.

## Algoritmo

A otimização para encontrar os parâmetros ótimos de um SVI dadas observações de mercado e a técnica de calibração **Quasi-explicit** de @DeMarco2009 pode ser resumida nos seguintes passos:

1. Definir valores iniciais para os parâmetros $(m, \sigma)$,

2. Iniciar algum otimizador global com estes parâmetros e resolver o problema completo \@ref(eq:completo)

    2.1 Dentro da otimização global, resolver o problema reduzido \@ref(eq:reduzido) para os parâmetros $(m, \sigma)$ dados,
  
3. Na convergência do problema completo do passo 2, otimizar uma última vez o problema reduzido,

4. Recuperar os parâmetros $(a, b, \rho, m, \sigma)$

A escolha dos otimizadores fica a cargo pessoal, sendo sugerido testar vários para o mesmo problema. Eventualmente, para um determinado _smile_ um otimizador pode se mostrar melhor que outro que vinha sendo utilizado em outras ocasiões.

Nos testes realizados pelo [Clube de Finanças](http://clubedefinancas.com.br), entre os otimizadores globais para o problema completo utilizamos [**Algoritmos Genéticos**](https://cran.r-project.org/package=GA), [**Nelder-Mead restrito**](https://www.rdocumentation.org/packages/stats/versions/3.5.2/topics/constrOptim) e um [**Método não-linear generalizado**](https://cran.r-project.org/package=Rsolnp). Para o problema reduzido, apesar de ser linear, o método de Nelder-Mead restrito se mostrou tão eficiente quanto e de mais fácil implementação. Se o objetivo for fazer uma calibração direta, dos cinco parâmetros ao mesmo tempo, uma combinação de otimizador global em primeiro estágio e o método de Levenberg-Marquardt restrito para refinamento da solução é o ideal. 

## Resultados

A seguir apresentamos um _smile_ de referência para a calibração, obtido de [ivolatility.com](http://www.ivolatility.com/doc/usa/IV_Raw_Delta_surface.csv) e então partimos para diferentes técnicas de calibração de um SVI. Os códigos em [R](https://cran.r-project.org/) também estão apresentados ao longo do texto para melhor compreensão e estudo do leitor.

Os dados utilizados neste exemplo estão apresentados na tabela \@ref(tab:dados) abaixo. Esta é uma típica apresentação de um _slice_ de superfície, ou seja, dados para um _smile_ apenas. As principais variáveis são: a data em que os dados foram coletados (date), o preço de fechamento do ativo (stock_price), o prazo para expiração em dias (period), e medidas de moneyness como delta, _strike_ e o próprio _moneyness_, além é claro da volatilidade implícita (iv) retirada do mercado.

<table class="table table-striped" style="font-size: 14px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:dados)Dados reais para exemplo de calibração de uma SVI.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> date </th>
   <th style="text-align:left;"> symbol </th>
   <th style="text-align:left;"> exchange </th>
   <th style="text-align:right;"> stock_price_for_iv </th>
   <th style="text-align:right;"> period </th>
   <th style="text-align:right;"> delta </th>
   <th style="text-align:right;"> moneyness </th>
   <th style="text-align:right;"> strike </th>
   <th style="text-align:right;"> iv </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 148.41 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 147.49 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 146.80 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 146.21 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 145.69 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 145.19 </td>
   <td style="text-align:right;"> 0.10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 144.69 </td>
   <td style="text-align:right;"> 0.10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 144.18 </td>
   <td style="text-align:right;"> 0.10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 143.66 </td>
   <td style="text-align:right;"> 0.10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 143.12 </td>
   <td style="text-align:right;"> 0.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> 142.53 </td>
   <td style="text-align:right;"> 0.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> 141.88 </td>
   <td style="text-align:right;"> 0.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> 141.13 </td>
   <td style="text-align:right;"> 0.12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> -0.02 </td>
   <td style="text-align:right;"> 140.26 </td>
   <td style="text-align:right;"> 0.13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> -0.03 </td>
   <td style="text-align:right;"> 139.16 </td>
   <td style="text-align:right;"> 0.13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> -0.04 </td>
   <td style="text-align:right;"> 137.66 </td>
   <td style="text-align:right;"> 0.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2017-09-21 </td>
   <td style="text-align:left;"> IWM </td>
   <td style="text-align:left;"> NYSEArca </td>
   <td style="text-align:right;"> 143.73 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> -0.06 </td>
   <td style="text-align:right;"> 135.32 </td>
   <td style="text-align:right;"> 0.16 </td>
  </tr>
</tbody>
</table>

Esta tabela poderia conter (de fato contém no arquivo original) outros períodos de expiração, e neste caso uma das colunas de _moneyness_ começa a se repetir, no caso seria o delta pois baixamos uma tabela de volatilidades implícitas por delta. Assim, em uma tabela simples em formato [`tidy`](http://vita.had.co.nz/papers/tidy-data.html) é possível armazenar informações de uma superfície inteira, a qual de outra forma necessitaria de um arranjo em 3 dimensões.

Ressaltamos aqui que a unidade de volatilidade implícita está em percentuais **ao ano**, equanto que nosso período é de dias corridos. É necessário harmonizar estas medidas de forma que, para volatiliades dadas em percentual ao ano, o período também seja dado em **anos**. Logo nosso $\tau = 30/365$, ou seja, 0.08219. 

Demonstraremos aqui os resultados para a calibração de uma _RAW SVI_ pelos métodos "Direto", "GA", "Quasi-NM" e "Quasi-PQ", abaixo explicados.

O método "Direto" é uma calibração direta através de um algoritmo de **Levenberg-Marquardt** da equação \@ref(eq:rawsvi), ou seja, não existe reparametrização _Quasi-explicit_ e o problema resolvido é não-linear em 5 dimensões. São realizadas 10 calibrações com estimativas iniciais dos parâmetros aleatórias, mas dentro de seus respectivos domínios. A melhor solução, aquela com o menor valor para a função objetivo, é selecionada. Todos os outros métodos utilizam a reparametrização, ocorrendo variações apenas nos algoritmos de otimização utilizados nos problemas reduzido e completo.

A calibração "GA" faz uso do otimizador global de **algoritmos genéticos** para o problema completo, ou seja, para estimar o par $(m, \sigma)$ que corresponde ao mínimo global. Após, o problema reduzido é resolvido através do algoritmo de **Nelder-Mead**. Este método é robusto, pois o algoritmo genético tem grande probabilidade de encontrar a região onde se encontra o mínimo global e não ficar preso localmente. Entretanto a robustez ocorre as expensas do tempo de computação.

Os métodos ditos "Quasi" diferem entre si na resolução do problema reduzido. Enquanto "PQ" remete a **programação quadrática** e faz uso da resolução do sistema linear apresentado na equação \@ref(eq:linear) com as restrições impostas por \@ref(eq:D), o método "Quasi-NM" utiliza o método de **Nelder-Mead** com restrições para a resolução deste mesmo problema reduzido. Em ambos os métodos, o problema completo é resolvido com um algoritmo de Nelder-Mead com 50 reinicializações das estimativas iniciais dos parâmetros $(m, \sigma)$, o que causa algum impacto no tempo de computação destes métodos.


```r
smile <- dados %>% 
  mutate(tau = period / 365) %>% 
  select(moneyness, iv, tau)

par_names <- factor(c("a", "b", "$\\rho$", "m", "$\\sigma$"),
                    levels = c("a", "b", "$\\rho$", "m", "$\\sigma$"))
k <- smile$moneyness
w <- smile$iv^2 * smile$tau

init_direct <- proc.time()
par_direct <- svi_fit_direct(k, w)
end_direct <- proc.time()
time_direct <- end_direct - init_direct 

init_ga <- proc.time()
par_ga <- svi_fit_ga(k, w)
end_ga <- proc.time()
time_ga <- end_ga - init_ga 

init_quasipq <- proc.time()
par_quasipq <- svi_fit_quasi(k, w, inner = "quadprog")
end_quasipq <- proc.time()
time_quasipq <- end_quasipq - init_quasipq 

init_quasinm <- proc.time()
par_quasinm <- svi_fit_quasi(k, w)
end_quasinm <- proc.time()
time_quasinm <- end_quasinm - init_quasinm 

iv_direct <- sqrt(svi_fun(par_direct$par[[1]], k) / smile$tau)
iv_ga <- sqrt(svi_fun(par_ga$par[[1]], k) / smile$tau)
iv_quasipq <- sqrt(svi_fun(par_quasipq$par[[1]], k) / smile$tau)
iv_quasinm <- sqrt(svi_fun(par_quasinm$par[[1]], k) / smile$tau)

plot_tbl <- tibble(k = k,
              Direct = iv_direct,
              GA = iv_ga,
              QuasiPQ = iv_quasipq,
              QuasiNM = iv_quasinm,
              observed = smile$iv) %>% 
  gather(key = method, value = iv, -c(k, observed))

par_tbl <- bind_rows(par_direct, par_ga, par_quasipq, par_quasinm) %>% 
  select(method, par) %>% 
  mutate(method = c("Direct", "GA", "QuasiPQ", "QuasiNM")) %>% 
  unnest() %>% 
  mutate(names = rep(par_names, 4)) %>% 
  spread(method, par) %>% 
  select(names, Direct, GA, QuasiPQ, QuasiNM) %>% 
  mutate(names = as.character(names)) %>% 
  mutate_at(vars(Direct:QuasiNM), arred)

rmse_tbl <- bind_rows(par_direct, par_ga, par_quasipq, par_quasinm) %>% 
  select(method, par) %>% 
  mutate(method = c("Direct", "GA", "QuasiPQ", "QuasiNM")) %>% 
  unnest() %>% 
  group_by(method) %>% 
  summarise(RMSE = rmse(par, k, w)) %>% 
  spread(method, RMSE) %>% 
  mutate(names = "RMSE") %>% 
  select(names, Direct, GA, QuasiPQ, QuasiNM) %>% 
  mutate_at(vars(Direct:QuasiNM), format, digits = 3, scientific = TRUE)

time_tbl <- tibble(method = c("Direct", "GA", "QuasiPQ", "QuasiNM"),
                   time = rbind(time_direct, time_ga, 
                                time_quasipq, time_quasinm)[, 3]) %>% 
  spread(method, time) %>% 
  mutate(names = "Tempo") %>% 
  select(names, Direct, GA, QuasiPQ, QuasiNM) %>% 
  mutate_at(vars(Direct:QuasiNM), arred)

frame_tbl <- bind_rows(par_tbl, rmse_tbl, time_tbl)
```

Abaixo é apresetanda uma tabela com os valores estimados para os parâmetros da SVI, o RMSE (root mean square error) e o tempo total em segundos para a calibração. Aqui o RMSE é definido como $\sqrt{1/n\sum(w(k_i)-\tilde w_i)^2}$ e nos fornece um valor típico de erro **na variância**.
 

```r
kable(frame_tbl,
      col.names = c("Estimativa", "Direto", "GA", 
                    "QuasiPQ", "QuasiNM"),
      caption = "Parâmetros estimados da calibração, RMSE e tempo de computação em segundos.",
      booktabs = TRUE) %>% 
  kable_styling(bootstrap_options = "striped",
                font_size = 18,
                full_width = FALSE)
```

<table class="table table-striped" style="font-size: 18px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:tabela)Parâmetros estimados da calibração, RMSE e tempo de computação em segundos.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Estimativa </th>
   <th style="text-align:left;"> Direto </th>
   <th style="text-align:left;"> GA </th>
   <th style="text-align:left;"> QuasiPQ </th>
   <th style="text-align:left;"> QuasiNM </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> a </td>
   <td style="text-align:left;"> 0.00000 </td>
   <td style="text-align:left;"> 0.00000 </td>
   <td style="text-align:left;"> 0.00000 </td>
   <td style="text-align:left;"> 0.00001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> b </td>
   <td style="text-align:left;"> 0.01940 </td>
   <td style="text-align:left;"> 0.01959 </td>
   <td style="text-align:left;"> 0.01696 </td>
   <td style="text-align:left;"> 0.01946 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\rho$ </td>
   <td style="text-align:left;"> -0.77790 </td>
   <td style="text-align:left;"> -0.79978 </td>
   <td style="text-align:left;"> -1.00000 </td>
   <td style="text-align:left;"> -0.79674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> m </td>
   <td style="text-align:left;"> -0.00595 </td>
   <td style="text-align:left;"> -0.00777 </td>
   <td style="text-align:left;"> -0.01118 </td>
   <td style="text-align:left;"> -0.00723 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> $\sigma$ </td>
   <td style="text-align:left;"> 0.04947 </td>
   <td style="text-align:left;"> 0.05034 </td>
   <td style="text-align:left;"> 0.06904 </td>
   <td style="text-align:left;"> 0.04998 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RMSE </td>
   <td style="text-align:left;"> 8.73e-06 </td>
   <td style="text-align:left;"> 8.69e-06 </td>
   <td style="text-align:left;"> 9.55e-05 </td>
   <td style="text-align:left;"> 8.69e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tempo </td>
   <td style="text-align:left;"> 0.17000 </td>
   <td style="text-align:left;"> 30.16900 </td>
   <td style="text-align:left;"> 0.18400 </td>
   <td style="text-align:left;"> 13.15000 </td>
  </tr>
</tbody>
</table>

O método Direto, com algoritmo de Levenberg-Marquardt se mostrou muito mais rápido que os demais, principalmente com relação ao algoritmo genético, e com um bom ajuste dado o baixo valor de RMSE. O algoritmo genético é consideravelmente mais lento, entretanto durante as várias calibrações realizadas em testes (e que não estão apresentadas na tabela \@ref(tab:tabela)), este algoritmo sempre se mostrou robusto, com baixo RMSE, diferentemente dos outros métodos que por vezes, denpendendo das estimativas iniciais, podem convergir para um mínimo local.

O gráfico com os ajustes realizados pode ser observado abaixo.

(ref:grafico) Comparação entre diferentes métodos de calibração de uma SVI.


```r
ggplot(plot_tbl, aes(x = k)) + 
  geom_point(aes(y = observed)) +
  geom_line(aes(y = iv, color = method)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "",
       x = "Forward log-moneyness",
       y = "Volatility",
       caption = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_viridis_d() +
  theme_economist_white()
```

<div class="figure">
<img src="07-calibrando-uma-svi_files/figure-epub3/plot-1.png" alt="(ref:grafico)"  />
<p class="caption">(\#fig:plot)(ref:grafico)</p>
</div>

De fato o método direto e o método _Quasi-explicit_ com otimizador global do tipo algoritmo genético se mostram mais adequados para a calibração de um SVI. Enquanto o método direto é muito mais eficiente em termos computacionais, o método _Quasi-explicit_ com GA é mais robusto. Desta forma, deve-se salientar que é necessário que o usuário, ao fazer uma calibração de _smile_ de volatilidade, deve dispor de diferentes métodos de fazê-lo, e a inspeção visual do resultado é **obrigatória** para determinar qual método foi mais eficiente em ajustar a curva aos dados.

## Conclusão

Apesar de neste exemplo ter se mostrado um método efetivo, com bom ajuste e baixo tempo de calibração, o método direto é altamente dependente dos valores iniciais dos parâmetros. Para tornar este método mais robusto, um número maior de reinicializações deve ser feita o que penaliza o tempo de calibração. O método _Quasi-explicit_ com algoritmo genético para encontrar a região de $(m, \sigma)$ onde se encontra o mínimo global se mostrou bastante robusta, entretanto, de convergência lenta. Para ajustar apenas um _smile_ alguns segundos a mais não representam problema. Porém, se imaginarmos que em uma grande instituição financeira são necessárias calibrações de, talvez, milhares de _smiles_ representando inúmeras superfícies de diversos instrumentos, este método pode se mostrar computacionalmente caro.

Já os métodos _Quasi-explicit_ que utilizam um algoritmo de Nelder-Mead para a resolução do problema completo se mostraram muito sensíveis às estimativas iniciais dos parâmetros. Mesmo utilizando 50 reinicialzações do método, diversas vezes o ajuste realizado foi insatisfatório. A resolução através de programação quadrática é rápida, se comparada com NM, entretanto, quando as restrições impostas pela equação \@ref(eq:D) se tornam ativas, este método parece sofrer com algum viés em sua solução.


[^71]: A resolução desta derivada é uma simples regra da divisão, entretanto a simplificação do resultado pede alguma manipulação algébrica. É possível utilizar sistemas de computação simbólica, o qual recomendamos o [SymPy](https://live.sympy.org/)
[^72]: Condições para ausência de arbitragem do tipo borboleta em um SVI estão detalhadas na seção 2.2 do artigo de @Gatheral2014.
