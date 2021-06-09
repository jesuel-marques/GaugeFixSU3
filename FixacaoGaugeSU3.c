///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Para compilar, usar no terminal	a linha abaixo
// gcc -o FixacaoGaugeSU3 FixacaoGaugeSU3.c Fonte/FuncoesFixacaoGaugeSU3.c Fonte/MedicoesFixacaoGaugeSU3.c Fonte/AlgoritmosFixacaoGaugeSU3.c Fonte/MedicoesSU3.c Fonte/SU2.c Fonte/SU3.c Fonte/RedeSU3.c Fonte/FuncoesMatematicas.c Fonte/ranlxd.c Fonte/ranlux_common.c -lm -O3 -DSSE2 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	beta (acoplamento inverso), o algoritmo de fixação e a tolerância para e2 estão definidos nesse arquivo.	///
//																												///
//	Algoritmos de fixação:																						///
//	SR é sobrerelaxação, SRE é sobrerelaxação estocástica e MD é máximo declive(ou Cornell)						///
//																												///
//	Os parâmetros de tamanho de rede e quantidade de configurações estão em ParametrosSU3.h						///
//	Parâmetros de fixação de calibre estão em ParametrosFixacaoGaugeSU3.h										///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>								//	Arquivos de cabeçalho padrão C
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "ParametrosSU3.h"						//	Parâmetros da simulação
#include "ParametrosFixacaoGaugeSU3.h"	 		//	Parâmetros dos algoritmos
#include "Fonte/GlobalSU3.h"					//	Definição de variáveis globais
#include "Fonte/GlobalFixacaoGaugeSU3.h"		//	Definição de variáveis globais específicas da fixação

#include "Fonte/SU3.h"							//	Funções de cálculo de SU3
#include "Fonte/RedeSU3.h"						//	Funções de inicialização, cálculo de posições na rede e vizinhos

#include "Fonte/FuncoesFixacaoGaugeSU3.h"		//	Funções específicas de fixação de calibre
#include "Fonte/MedicoesFixacaoGaugeSU3.h"		//	Funções de medição específicas de fixação de calibre
#include "Fonte/AlgoritmosFixacaoGaugeSU3.h"	//	Algoritmos de fixação de calibre

#include "Fonte/ranlux.h"						//	Gerador aleatório http://luscher.web.cern.ch/luscher/ranlux/


double beta=6.0;	//	Inverso do acoplamento

double U[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];		//	Matrix que conterá todos os links da rede. 
double Uaux[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];	//	Matrix com o links atualizados pela transformação de calibre. 

double G[Nt][Nxyz][Nxyz][Nxyz][3][3][2];		//	Transformação de calibre sobre U(x) para levar ao calibre de Landau


int main(){

	char nomearquivoConfig[ComprimentoMaxArquivo-200];	//	Nome do arquivo de configuração a ser fixada
	char nomearquivoConfigCompleto[ComprimentoMaxArquivo];
	char nomearquivoConfigGaugeFixo[ComprimentoMaxArquivo];	//	Nome do arquivo de configuração fixada 

	// char nomearquivovarreduras[ComprimentoMaxArquivo];
	// FILE *arquivovarreduras; 

	double tolerancia=pow(10.0,-12);	//	Tolerância da distância para o calibre de Landau

	//	Alteração no nome para registrar o algoritmo que fez a fixação de calibre

	int Algoritmo=SR;
	char AlgoritmoString[20];

	if(Algoritmo==LA)
		strcpy(AlgoritmoString,"LA");

	else if(Algoritmo==SR){
		strcpy(AlgoritmoString,"SR");
		// sprintf(nomearquivovarreduras,"Varreduras_%s_beta_%lf_Nxyz_%d_Nt_%d_passos_%d_wSR_%f.txt",AlgoritmoString,beta,Nxyz,Nt,maxhits,omegaSR);
	}
	else if(Algoritmo==SRE){
		strcpy(AlgoritmoString,"SRE");
		// sprintf(nomearquivovarreduras,"Varreduras_%s_beta_%lf_Nxyz_%d_Nt_%d_passos_%d_pSRE_%f.txt",AlgoritmoString,beta,Nxyz,Nt,maxhits,pSRE);
		rlxd_init(1,983604080);
	}
	else if(Algoritmo==MD){
		strcpy(AlgoritmoString,"MD");
		// sprintf(nomearquivovarreduras,"Varreduras_%s_beta_%lf_Nxyz_%d_Nt_%d_alfaMD_%lf.txt",AlgoritmoString,beta,Nxyz,Nt,alfaMD);

	}

	InicializarUSU3(FRIA); //	Inicializa matrizes identidade, R, S e T, necessárias para cálculos

	//	Carregando Arquivo de Configurações


	for(int config=1;config<=MAXconfigs;config++){

		//	Loop de configurações a serem fixadas

		printf("Config %d \n",config);

		//	Carregando nome do arquivo que contem configuração a ser fixada

		sprintf(nomearquivoConfig,"Config_%d_beta_%.3lf_Nxyz_%d_Nt_%d.txt",config,beta,Nxyz,Nt);
		sprintf(nomearquivoConfigCompleto,"ConfigsDescorrelacionadas/%s",nomearquivoConfig);
		CarregarUSU3(nomearquivoConfigCompleto);

		//	Função que fixa o calibre da configuração

		FixarGaugeSU3(config,Algoritmo,tolerancia);

		// arquivovarreduras=fopen(nomearquivovarreduras,"a+");

		// fprintf(arquivovarreduras,"%d\t%d\n",config,FixarGaugeSU3(config,Algoritmo,tolerancia));

		// fclose(arquivovarreduras);

		//	Gerando nome do arquivo com a configuração já fixada e eventual impressão da configuração

			sprintf(nomearquivoConfigGaugeFixo,"./ConfigsGaugeFixadas/%s_GaugeFixo_%s",AlgoritmoString,nomearquivoConfig);

	 	ImprimirUGaugeFixadoSU3(nomearquivoConfigGaugeFixo);
	
	}

	return 0;
}