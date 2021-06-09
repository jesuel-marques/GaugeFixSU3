#include<math.h>					//	Arquivo de cabeçalho padrão C

#include "../ParametrosSU3.h"		//	Parâmetros da simulação

#include "GlobalSU3.h"				//	Definição de variáveis globais
#include "FuncoesMatematicas.h"		//	Funções de aritmética complexa e algebra linear

double identidadeSU2[4]={1.0,0,0,0}; //	Identidade no espaço de matrizes SU(2).
double nuloSU2[4]={0.0,0.0,0.0,0.0}; // Elemento nulo no espaço das matrizes SU(2).

void CopiarSU2(double u[4], double v[4]){

//	Copiar a matriz u na matriz v.

	for(int a=0;a<=3;a++)
		v[a]=u[a];
}

void ProdutoSU2(double u[4], double v[4], double uv[4]){ 

//	Função que retorna o produto de duas matrizes SU(2) na representação de Cayley-Klein.

	double prodext[4]={0,0,0,0};

	uv[0]=ProdutoInterno(u,v);	//	A componente 0 é simplesmente dada pelo produto interno.

	ProdutoExterno(u,v,prodext);	//	As componentes 1, 2 e 3 são dadas pelo produto vetorial
									// 	somadas com um termo que as mistura com a componente 0. 
	for(int a=1;a<=3;a++)

		uv[a]=u[a]*v[0]+u[0]*v[a]-prodext[a];

}

void ProdutoSU2Tres(double a[4], double b[4], double c[4], double abc[4]){

//	Função que retorna o produto de três matrizes SU(2), utilizando recursivamente o produto de duas.

	double ab[4];
	ProdutoSU2(a,b,ab);
	ProdutoSU2(ab,c,abc);

}

void ProdutoSU2Quatro(double a[4], double b[4], double c[4], double e[4], double abce[4]){

//	Função que retorna o produto de quatro matrizes SU(2), utilizando recursivamente o produto de duas.

	double abc[4];
	ProdutoSU2Tres(a,b,c,abc);
	ProdutoSU2(abc,e,abce);

}


void PotenciaSU2(double u[4],double unapotencias[10][4]){

//	Calcula as potências de u e as coloca em diferentes linhas de unapotencias.

	for(int a=0;a<=3;a++){
		unapotencias[0][a]=identidadeSU2[a];	// u elevado a 0 é a identidade
		unapotencias[1][a]=u[a];	// u elevadoa 1 é si mesmo
	}

	for(int i=2;i<=9;i++)
		ProdutoSU2(unapotencias[i-1],u,unapotencias[i]);
	
}


void SomaSU2(double a[4], double resultado[4]){

//	Função que retorna a soma de duas matrizes SU(2) na representação de Cayley-Klein e coloca o  
//	resultado na segunda matriz.

	for(int i=0;i<=3;i++)

		resultado[i]+=a[i]; 

}

void DiferencaSU2(double a[4], double b[4], double amenosb[4]){

//	Função que retorna a subtração de duas matrizes SU(2) na representação de Cayley-Klein e coloca o 
//	resultado na segunda matriz.

	for(int i=0;i<=3;i++)

		amenosb[i]=a[i]-b[i];

}

double TrSU2(double u[4]){

//	Função que retorna o traço de uma matriz SU(2) na representação de Cayley-Klein.

	return 2.0*u[0];

}

void ConjHermSU2(double u[4],double udag[4]){

//	Função que retorna a matriz conjugada hermitiana de uma matriz SU(2) na representação de Cayley-Klein.

	udag[0]=u[0]; //	Componente 0 é a mesma.

	for(int i=1;i<=3;i++)

		udag[i]=-u[i];	//	Componentes 1, 2 e 3 mudam de sinal.

}

void MultiplicacaoEscalarSU2(double u[4],double num, double numvezesu[4]){

//	Função que retorna a matriz SU(2) na representação de Cayley-Klein multiplicada por um escalar.

	for(int a=0;a<=3;a++)
		numvezesu[a]=num*u[a];	//	Basta multiplicar cada componente.

}

double DeterminanteSU2(double u[4]){

//	Função que retorna o determinante de uma matriz SU(2) na representação de Cayley-Klein.

	double det=0.0;

	for(int i=0;i<=3;i++)

		det+=pow2(u[i]);	//	Determinante é dado pela soma dos quadrados das componentes.

	return det;

}

void ProjecaoSU2(double a[4],double aSU2[4]){

//	Projeta a matrix a no grupo SU(2) de maneira que tenha determinante 1.

	MultiplicacaoEscalarSU2(a,1.0/sqrt(DeterminanteSU2(a)),aSU2);

}