#include <stdio.h>				//	Arquivos de cabeçalho padrão C
#include <stdlib.h>

#include <math.h>

#include "../ParametrosSU3.h"	//	Parametros da simulação
#include "GlobalSU3.h"			//	Definição de variáveis globais

#include "FuncoesMatematicas.h"	//	Funções de aritmética complexa e algebra linear
#include "SU3.h"				//	Funções de cálculo de SU3
#include "RedeSU3.h"			//	Funções de inicialização, cálculo de posições na rede e vizinhos


void CalcularASU3(int posicao[d],int mu, double A[3][3][2]){

//	Função que calcula o campo de gauge a partir dos elos em Uaux, que é a matriz transformada de gauge do programa de fixação.
	double Udag[3][3][2];
	double UmenosUdag[3][3][2];
	double escalar1[2]={0.0,-0.5};
	double escalar2[2]={-(double)1.0/3.0,0.0};

	double traco[2];
	double fator[2];
	double partetraco[3][3][2];

	ConjHermSU3(Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],Udag);
	DiferencaSU3(Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],Udag,UmenosUdag);
	MultiplicacaoEscalarSU3(UmenosUdag,escalar1,A);

	//Tirar parte do traço

	TrSU3(A,traco);

	ProdutoComplexo(traco,escalar2,fator);
	MultiplicacaoEscalarSU3(identidadeSU3,fator,partetraco);// isso pode ser melhorado, toda vez calculando identidade vezes 1/Nc
	SomaSU3(partetraco,A);


}

void DivergenciaASU3(int posicao[d], double divA[3][3][2]){	//	Calcula a divergência na rede do campo de gauge.
	int versormenosmu[d]={0,0,0,0};
	int posicaomenosmu[d];

	double A1[3][3][2];	//	Calcula através da diferença entre As de elos separados por um ponto. 
	double A2[3][3][2];

	double termodivA[3][3][2];

	CopiarSU3(nuloSU3,divA);


	for(int mu=0;mu<d;mu++){

		versormenosmu[mu]=-1;
		SomaVetoresPosicao(posicao,versormenosmu,posicaomenosmu);
		
		CalcularASU3(posicao,mu,A1);		//melhorar esse cálculo... estou calculando A muitas vezes, melhor fazer varredura calculando A da rede inteira e depois só ler
		CalcularASU3(posicaomenosmu,mu,A2);

		DiferencaSU3(A1,A2,termodivA);

		SomaSU3(termodivA,divA);	//	Soma os termos correspondendo as diferentes direções mu.

		versormenosmu[mu]=0;
	}

}