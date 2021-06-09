#include <stdio.h>							//	Arquivos de cabeçalho padrão C
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../ParametrosSU3.h"				//	Parametros da simulação
#include "../ParametrosFixacaoGaugeSU3.h"	//	Parâmetros dos algoritmos
#include "GlobalSU3.h"						//	Definição de variáveis globais
#include "GlobalFixacaoGaugeSU3.h"			//	Definição de variáveis globais específicas da fixação

#include "FuncoesMatematicas.h"				//	Funções de cálculo SU(3) e gerais de rede
#include "SU3.h"							//	Funções de cálculo de SU3
#include "RedeSU3.h"						//	Funções de inicialização, cálculo de posições na rede e vizinhos
#include "MedicoesSU3.h"					//	Funções de medição de quantidades como A_mu(n) e a div A_mu(n)

#include "MedicoesFixacaoGaugeSU3.h"		//	Funções de medição do processo de fixação
#include "AlgoritmosFixacaoGaugeSU3.h"		//	Algoritmos de fixação de calibre de Landau


void CalcularwSU3(int posicao[d],double w[3][3][2]){	

	//	Calcula w(n)=g(n)*h(n), seguindo notação de hep-lat/0301019v2
	
	//	h(n)
	double h[3][3][2];

	//	Para cada mu haverá dois produtor de matrizes a serem feitos

	double prodfrente[3][3][2];	
	double prodtras[3][3][2];

	//	Matrizes necessárias para o cálculo de h(n)

	double gdaggerfrente[3][3][2];
	double gdaggertras[3][3][2];

	double udaggertras[3][3][2];

	//	Posições relevantes para o cálculo de h(n)

	int versormu[d]={0,0,0,0};
	int versormenosmu[d]={0,0,0,0};
	int posicaomaismu[d];
	int posicaomenosmu[d];

	//	Inicializando h(n)=0

	CopiarSU3(nuloSU3,h);

	//	Cálculo de h(n)

	for(int mu=0;mu<d;mu++){

		versormu[mu]=1;
		versormenosmu[mu]=-1;

		SomaVetoresPosicao(posicao,versormu,posicaomaismu);
		SomaVetoresPosicao(posicao,versormenosmu,posicaomenosmu);
		
		ConjHermSU3(G[posicaomaismu[0]][posicaomaismu[1]][posicaomaismu[2]][posicaomaismu[3]],gdaggerfrente);
		ConjHermSU3(G[posicaomenosmu[0]][posicaomenosmu[1]][posicaomenosmu[2]][posicaomenosmu[3]],gdaggertras);

		ConjHermSU3(U[posicaomenosmu[0]][posicaomenosmu[1]][posicaomenosmu[2]][posicaomenosmu[3]][mu],udaggertras);

		ProdutoSU3(U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],gdaggerfrente,prodfrente);
		SomaSU3(prodfrente,h);

		ProdutoSU3(udaggertras,gdaggertras,prodtras);
		SomaSU3(prodtras,h);

		versormu[mu]=0;
		versormenosmu[mu]=0;

	}


	//	w(x)=g(x)*h(x)

	ProdutoSU3(G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],h,w);

}

		
void AtualizarUauxSU3(){
	
	//	Função que faz o cálculo para atualização do U, que deve ser feito antes de calcular A

	//	Matrizes necessárias para o cálculo de Uaux

	double uatualizado[3][3][2];

	double gdaggerposicaomaismu[3][3][2];

	//	Posições relevantes para o cálculo de Uaux

	int versormu[d]={0,0,0,0};
	int posicao[d];
	int posicaomaismu[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
		for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
			for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
				for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
					for(int mu=0;mu<d;mu++){

						//	U'_mu(x)=g(x)U_mu(x)gdagger(x+mu)

						versormu[mu]=1;

						SomaVetoresPosicao(posicao,versormu,posicaomaismu);

						ConjHermSU3(G[posicaomaismu[0]][posicaomaismu[1]][posicaomaismu[2]][posicaomaismu[3]],gdaggerposicaomaismu);
						
						ProdutoSU3Tres(G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],gdaggerposicaomaismu,uatualizado);

						//	Projetar no grupo SU(3) por segurança
						ProjecaoSU3(uatualizado,Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu]);

						versormu[mu]=0;
					}
}

void AtualizarUauxLocalSU3(int posicao[d]){

	//	Matrizes necessárias para o cálculo de Uaux

	double uatualizado[3][3][2];
	
	int versormu[d]={0,0,0,0};
	int posicaomaismu[d];
	double gdaggerposicaomaismu[3][3][2];

	int versormenosmu[d]={0,0,0,0};
	int posicaomenosmu[d];
	double gdagger[3][3][2];

	for(int mu=0;mu<d;mu++){

		//	U'_mu(x)=g(x)U_mu(x)gdagger(x+mu)

		versormu[mu]=1;

		SomaVetoresPosicao(posicao,versormu,posicaomaismu);

		ConjHermSU3(G[posicaomaismu[0]][posicaomaismu[1]][posicaomaismu[2]][posicaomaismu[3]],gdaggerposicaomaismu);
		
		ProdutoSU3Tres(G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],gdaggerposicaomaismu,uatualizado);

		//	Projetar no grupo SU(3) por segurança
		ProjecaoSU3(uatualizado,Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu]);

		versormu[mu]=0;

		//	U'_mu(x-mu)=g(x-mu)U_mu(x-mu)gdagger(x)

		versormenosmu[mu]=-1;

		SomaVetoresPosicao(posicao,versormenosmu,posicaomenosmu);

		ConjHermSU3(G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],gdagger);

		ProdutoSU3Tres(G[posicaomenosmu[0]][posicaomenosmu[1]][posicaomenosmu[2]][posicaomenosmu[3]],U[posicaomenosmu[0]][posicaomenosmu[1]][posicaomenosmu[2]][posicaomenosmu[3]][mu],gdagger,uatualizado);
		
		//	Projetar no grupo SU(3) por segurança
		ProjecaoSU3(uatualizado,Uaux[posicaomenosmu[0]][posicaomenosmu[1]][posicaomenosmu[2]][posicaomenosmu[3]][mu]);
		
		versormenosmu[mu]=0;
	}	

}

int FixarGaugeSU3(int config,int Algoritmo, double tolerancia){	//	Fixa o calibre de Landau
	//	FILE *arquivoe2e6;	//	Arquivo que guarda medição do parâmetro de fixação, que deverá ser menor que a tolerância

	char AlgoritmoString[20];	//	String com o algoritmo a ser utilizado para a fixação
	
	//char nomearquivoe2e6[ComprimentoMaxArquivo]; //	Nome do arquivo de medida de e2 e e6 (definido em hep-lat/0301019v2)

	//char printe2[100];
	//char printe6[100];

	int posicao[d];	//	Vetor de posição, argumento de g(x)

	int	cont=0;	//	Contador de varreduras para fixar o calibre de Landau
	double e2, e6;	//	Parâmetro de fixação, que deverá ser menor que a tolerância, seguindo notação de hep-lat/0301019v2

	double w[3][3][2];	//	w(n)=g(n)*h(n), seguindo notação de hep-lat/0301019v2

	int impressoe2=NAO,impressoe6=NAO;	//	flags de impressão de e2 e e6, ativados após essas quantidades atingirem tolerancia



	//	Modificador do nome do arquivo de medição e2

	if(Algoritmo==LA)
		strcpy(AlgoritmoString,"LA");
	else if(Algoritmo==SR){
		strcpy(AlgoritmoString,"SR");
		//sprintf(nomearquivoe2e6,"e2e6_%s_omega_%.2lf_hits_%d_Config_%d_beta_%.3lf_Nxyz_%d_Nt_%d.txt",AlgoritmoString,omegaSR,maxhits,config,beta,Nxyz,Nt);
	}
	else if(Algoritmo==SRE){
		strcpy(AlgoritmoString,"SRE");
		//sprintf(nomearquivoe2e6,"e2e6_%s_p_%.2lf_hits_%d_Config_%d_beta_%.3lf_Nxyz_%d_Nt_%d.txt",AlgoritmoString,pSRE,maxhits,config,beta,Nxyz,Nt);
	}
	else if(Algoritmo==MD)
		strcpy(AlgoritmoString,"MD");
	

	// printf("%s\n",nomearquivoe2e6);
	// arquivoe2e6=fopen(nomearquivoe2e6,"w+");


	//	Medição da distância inicial para condição de fixação de calibre

	e2=Calculare2SU3();
	e6=Calculare6SU3();


	do{	
		for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
			for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
				for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
					for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++){	

						
							if(Algoritmo!=MD)
								CalcularwSU3(posicao,w);	//	Calculo de w(n)=g(n)*h(n)

																		
							
							switch(Algoritmo){	//	Algoritmos de fixação de gauge de Landau, recebem w(n) como entrada
								
								case MD:
									MaximoDecliveSU3(posicao);
									break;

								case LA:
									LosAlamosSU3(posicao,w);	
									break;

								case SR:
									SobreRelaxacaoSU3(posicao,w);
									break;

								case SRE:
									SobreRelaxacaoEstocasticaSU3(posicao,w);
									break;

								default:
									printf("Definir Algoritmo\n");
									break;

									
							}

						
						
					}

		cont++;
		


		

		if(cont%varredurasparamedicaoe2==0){	//	Pode fazer cálculo de e2 a cada 10 varreduras, por exemplo

			if(Algoritmo!=MD)	//	Atualização não é necessária para MD porque estou fazendo atualização local em cada passo
				AtualizarUauxSU3();	//	Após atualização de g(x), Uaux deve ser atualizado para cálculo de A

			e2=Calculare2SU3();	//	Cálculo de quão longe está da condição de gauge de Landau, e2=0
			e6=Calculare6SU3();

			printf("\n%d e2: %3.2e\t e6: %3.2e\n",cont,e2,e6);

		}
		

		// if(e2<=tolerancia&&impressoe2==NAO){ 
		// 	sprintf(printe2,"%d\t%lf\t",cont,log10(e6));
		// 	impressoe2=SIM;
		// }

		// if(e6<=tolerancia&&impressoe6==NAO){
		// 	sprintf(printe6,"%d\t%lf\t",cont,log10(e2));
		// 	impressoe6=SIM;
		// }

		
	}while(/*e6>tolerancia||*/e2>tolerancia);	//	Enquanto e2 for maior do que a tolerância, repete o processo iterativamente
	
							
	printf("Varreduras necessárias para fixação: %d \n",cont);

	// fprintf(arquivoe2e6,"%d\t",config);
	// fprintf(arquivoe2e6,"%s\t%s\n",printe2,printe6);
	

	// fclose(arquivoe2e6);


	return cont;

}

void ImprimirUGaugeFixadoSU3(char NomeArquivo[ComprimentoMaxArquivo]){	

//	Imprime em um arquivo a configuração atual de elos, após o processo de fixação de calibre de Landau.

	FILE *Arquivo;


	AtualizarUauxSU3(); 

	Arquivo=fopen(NomeArquivo,"w+");

	fwrite(Uaux, sizeof(Uaux),1, Arquivo);

	fclose(Arquivo);
		
}