# Introdução as opções {#opcoes}



Neste Capítulo[^11] iremos apresentar os instrumentos financeiros conhecidos como opções. Existem dois tipos de opção. Uma opção de compra dá ao detentor o direito de comprar o ativo subjacente até uma determinada data por um determinado preço. Uma opção de venda dá ao titular o direito de vender o ativo subjacente até uma determinada data por um determinado preço.

Alguns termos importantes para o melhor entendimento deste Capítulo são: “Ativo subjacente”, que é o ativo negociado no contrato, “data de vencimento”, no caso do modelo americano é até a data limite para exercer a opção de compra e no modelo europeu é na data final em que a opção de compra pode ou não ser exercida (serão discutidos mais detalhes sobre estas modalidades posteriormente) e “preço de exercício” (strike), é o valor a ser pago pelo ativo de acordo com o contrato. As opções americanas podem ser exercidas a qualquer momento até a data de vencimento. Por sua vez, as opções europeias podem ser exercidas somente na própria data de vencimento.

As opções europeias são geralmente mais fáceis de analisar do que as opções americanas, e algumas das propriedades de uma opção americana são frequentemente deduzidas daquelas de sua contraparte europeia.

A opção de compra (call) dá ao comprador da opção o **direito de comprar** o ativo subjacente pelo preço de exercício. A opção de venda (put) dá ao comprador da opção o **direito de vender** o ativo subjacente pelo preço combinado no contrato na data futura.

O preço de uma opção de compra diminui à medida que o preço de exercício
aumenta, enquanto o preço de uma opção de venda aumenta à medida que o preço de
exercício aumenta. Ambos os tipos de opção tendem a se tornar mais valiosas à medida que seu tempo até o vencimento aumenta. Na verdade existem seis fatores que afetam o preço de uma opção de ação, por exemplo:

1. O preço atual da ação, $S_0$
2. O preço de exercício, $K$
3. O tempo para expiração, $\tau$
4. A volatilidade do preço das ações, $\sigma$
5. A taxa de juros livre de risco, $r$
6. Os dividendos que se espera sejam pagos, $q$

Observamos que existem quatro tipos de participantes nos mercados de opções:

- Compradores de calls
- Vendedores de calls
- Compradores de puts
- Vendedores de puts

Os compradores são referidos como tendo posições _long_, vendedores são referidos
como tendo posições _short_. Vender uma opção também é conhecido como lançar a opção ou subscrever.

Exemplificando uma operação de compra de call, caso o preço do ativo tenha subido acima do preço de strike o comprador pode exercer sua opção de compra e ele lucrará a partir do momento em que o valor da ação for maior que o strike mais o valor pago pela opção (chamado de prêmio). Caso o preço do ativo tenha caído abaixo do strike, o comprador simplesmente não exerce sua opção, limitando sua perda nessa operação ao prêmio pago. 

\begin{table}[t]

\caption{(\#tab:opcao)Exemplo: Comprando call de X com o preço de exercício (strike) de 10.000}
\centering
\begin{tabular}{rrrrl}
\toprule
Ação & Prêmio Pago & Opção no Exercício & Lucro & Moneyness\\
\midrule
9000 & 200 & 0 & -200 & OTM\\
9500 & 200 & 0 & -200 & OTM\\
10000 & 200 & 0 & -200 & ATM\\
10200 & 200 & 200 & 0 & ITM\\
10500 & 200 & 500 & 300 & ITM\\
\addlinespace
11000 & 200 & 1000 & 800 & ITM\\
\bottomrule
\end{tabular}
\end{table}

Usando a tabela \@ref(tab:opcao) como exemplo é possível ver que o resultado final será maior que R\$ 0 quando o valor do ativo subjacente é maior que R\$10.200 (R\$10.000 de strike + R\$200 de prêmio), e o prejuízo final está limitado ao valor de R\$200. A figura \@ref(fig:call) abaixo mostra o perfil de lucro da operação exemplificada, típico de uma compra de call.

![(\#fig:call)Perfil de lucro típico de uma compra de call.](01-introducao-as-opcoes_files/figure-latex/call-1.pdf) 

No caso de uma operação de uma compra de put, caso o preço do ativo tenha subido acima do strike, não faz sentindo o detentor da opção exercer seu direito, assim sua perda será apenas o valor pago pela opção. Caso o preço do ativo tenha descido abaixo do strike, o comprador da put pode realizar a venda e começará a lucrar a partir do momento em que o strike fique acima do valor do ativo somado ao valor do prêmio pago pela opção. Usando a tabela 1.2 como exemplo é possível ver que o resultado final será de no mínimo -R\$600 caso o valor do ativo subjacente seja igual ou maior que R\$15.000, e o resultado final aumenta conforme o ativo perde o valor, sendo positivo a partir de quando seu valor é de R\$14.400.

## Conceitos in the money, at the money e out the money

Estes termos são usados para se referir a opções que estão com o preço de exercício (strike) do ativo abaixo, acima ou igual ao valor atual do ativo.

- Out the money: Strike do ativo subjacente está abaixo do valor de mercado no caso de calls ou quando o strike do ativo está acima do valor de mercado no caso de puts.   
- At the money: Strike do ativo subjacente é o igual ao valor de mercado.
- In the money: Strike do ativo subjacente está acima do valor de mercado no caso de calls ou quando o strike do ativo está abaixo do valor de mercado no caso de puts.

## Modelos americano e europeu de opções

No modelo americano de opções o comprador pode exercer seu direito de compra ou venda do ativo subjacente a qualquer momento entre o início do contrato e o vencimento dele, enquanto isso no modelo europeu a transação só pode ser realizada na data de vencimento.

## Hedge 

O mercado de opções pode ser usado tanto para hedge (proteção) quanto para especulação. O hedge é feito para limitar as possíveis perdas que um investidor pode ter ao estar com seu patrimônio atrelado a determinado ativo, por exemplo, para um acionista que possui ações de determinada empresa se proteger contra uma possível queda no valor de suas ações, ele pode comprar opções de venda at the money de suas ações para que seu prejuízo máximo seja o prêmio. 

## Travas 

Devido ao mercado de opções nos oferecer diversas possibilidades entre call e put onde você pode estar comprado e/ou vendido irá surgir várias posições a serem assumidas, para nos adequarmos ao quanto estamos dispostos a encarar o risco parar atingirmos o retorno desejado. Essas posições são conhecidas como “Travas”.

Entendendo as travas, existem diversas estratégias, como por exemplo:  Trava de alta, trava de baixa, Long Straddle, Short Straddle, entre outras. Mas afinal qual é o funcionamento delas? Supondo que o leitor espere uma alta do mercado, no entanto acredite que não irá superar determinado ponto ele poderá realizar uma Trava de alta. Onde comprará uma opção de Call a um preço de strike X e vender outra Call com o preço de strike Y, onde obrigatoriamente Y>X. Nesta operação limitaremos o nosso ganho caso o mercado supere as nossas expectativas, no entanto diminuiremos o custo da operação, o custo será o prêmio pago pela opção X menos o prêmio recebido pela opção Y, para facilitar a compreensão observemos o gráfico a seguir:
 
Nesse caso o valor do prêmio da compra foi de R\$30,00 enquanto a da venda foi R\$10,00, assim limitamos nossa perda em R\$20,00, enquanto os preços de strike da compra e da venda da call foram respectivamente R\$250,00 e R\$300,00, fazendo o retorno máximo ser R\$30,00 que é a diferença entre os valores de strike a serem realizados e descontado o valor pago pelo prêmio.

Agora que o leitor já entendeu melhor o conceito da trava, vamos explorar uma mais complexa a Long Butterfly. Aqui é realizado a compra de uma put e call com preços de strike iguais, vendesse uma put com preço inferior e vende alguma call com preço superior as iniciais. Observe que pelo fato de contar com a venda de duas opções nessa estratégia tem um custo de operação reduzido, no entanto o ideal é utilizar em um mercado de pouca volatilidade, dado que se a volatilidade ser alta perdesse a possibilidade de ter um ganho maior, nesse caso recomenda o uso por exemplo de uma Long Straddle. Enfim vamos ao gráfico para facilitar a compreensão da estratégia: 

Teremos então a compra de uma put e call de strike iguais de R\$150,00 a venda de uma put com strike inferior de R\$80, e a venda de uma call com preço superior R\$220. Os valores exatos dos prêmios não nos interessam no momento, porém é importante entender que teremos dois com saldos positivos referente a nossa venda e dois negativos que advém das compras, o resultado será nosso prejuízo máximo, olhando o gráfico nesse caso é de R\$20,00. Pela área de retorno do gráfico podemos ver que nosso risco está reduzido. Onde o pior cenário possível se encontra em o preço do ativo-objeto se aproximar do valor do R\$150,00, que é onde os contratos adquiridos não serão vantajosos em nenhuma ponta. No entanto se o preço se aproximar de R\$220 poderemos exercer nosso direito da compra da call inicial por R\$150,00(você terá o direito de comprar a um preço inferior), o mesmo será valido caso haja uma queda do preço se aproximado do valor de R\$80,00 onde a logica será a mesma só que aqui será o usado o direito da compra da put por R\$150,00(Você terá o direito de vender a um preço superior). Observe que apesar do nosso risco ser reduzido, limita os nossos ganhos, com a venda da call e da put com preços superior e inferior respectivamente.

Espero que o leitor tenha despertado interesse no assunto, com esse conteúdo dominado já saberá o básico sobre opções, fique atento a novas postagens em breve iremos mais a afundo explicando por exemplo o modelo Black Scholes, como as opções são precificadas entre outros materiais.

[^11]: Artigo originalmente escrito por Glauber Naue e Erik Kawano, adaptado por Rafael Bressan para este livro.
