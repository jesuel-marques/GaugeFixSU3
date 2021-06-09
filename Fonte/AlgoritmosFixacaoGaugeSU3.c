#include <stdio.h>								//	Arquivo de cabeçalho padrão C

#include "../ParametrosSU3.h"					//	Parametros da simulação 
#include "../ParametrosFixacaoGaugeSU3.h"	 	//	Parâmetros fixação de calibre
#include "GlobalSU3.h"							//	Definição de variáveis globais
#include "GlobalFixacaoGaugeSU3.h"				//	Definição de variáveis globais específicas para fixação de calibre


#include "SU2.h"								//	Funções para cálculo SU2
#include "SU3.h"								//	Funções para cálculos SU3

#include "FuncoesMatematicas.h"					//	Funções de aritmética complexa e algebra linear
#include "FuncoesFixacaoGaugeSU3.h"				//	Funções específicas de fixação de calibre

#include "MedicoesSU3.h"						//	Funções de medição de quantidades como A_mu(n) e a div A_mu(n)
#include "MedicoesFixacaoGaugeSU3.h"			//	Medições de quantidades da fixação de gauge como e2, e6, etc.

#include "ranlux.h"								//	Gerador de números aleatórios  http://luscher.web.cern.ch/luscher/ranlux/, usado no SRE

void MaximoDecliveSU3(int posicao[d]){

	//	Algoritmo de Cornell descrito em hep-lat/0301019v2

	double ialfa[2]={0.0,alfaMD};
	double quaddivA[3][3][2];
	double ialfaquaddivA[3][3][2];
	double Raprojetar[3][3][2];
	double R[3][3][2];
	double g[3][3][2];

	DivergenciaASU3(posicao,quaddivA);

	MultiplicacaoEscalarSU3(quaddivA,ialfa,ialfaquaddivA);
	DiferencaSU3(identidadeSU3,ialfaquaddivA,Raprojetar);
	ProjecaoSU3(Raprojetar,R);
	ProdutoSU3(R,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],g);
	CopiarSU3(g,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]]);

	AtualizarUauxLocalSU3(posicao);
	
}

void BlocoComumLASU3(double w[3][3][2],double A[3][3][2]){

	//	Calcula a matriz de atualização de calibre A a partir de w(n)=g(n).h(n) como no algoritmo de Los Alamos
	//	para SU(3), com divisão da matriz em submatrizes seguindo Cabbibo-Marinari. 

	double Omega[3][3][2];	//	Inicialmente é o complexo conjugado de w, depois complexo conjugado de Rw, depois de SRw, ...
	double Omegadagger[3][3][2];	//	Omega transposto conjugado
	
	double winv[3][3][2];	//	Inversa de w(x)
	double Rw[3][3][2];		//	Primeira atualização para w(x)
	double SRw[3][3][2];	//	Segunda atualização para w(x)	
	double TSRw[3][3][2];	//	Terceira atualização para w(x)


	//	Submatrizes de Omega e projeções para SU(2);

	double wR[4], wRSU2[4];	
	double wS[4], wSSU2[4];	
	double wT[4], wTSU2[4];


	ConjHermSU3(w,Omega);


	for(int hits=1;hits<=maxhits;hits++){	//	Cada hit contém a subdivisão de Cabbibo-Marinari

		ConjHermSU3(Omega,Omegadagger); 

		//	Primeira Submatriz

		wR[0]=(Omega[0][0][RE]+Omega[1][1][RE]);
		wR[1]=(Omega[0][1][IM]+Omega[1][0][IM]);
		wR[2]=(Omega[0][1][RE]-Omega[1][0][RE]);
		wR[3]=(Omega[0][0][IM]-Omega[1][1][IM]);


		ProjecaoSU2(wR,wRSU2);
		
		R[0][0][RE]=wRSU2[0];
		R[0][0][IM]=wRSU2[3];

		R[0][1][RE]=wRSU2[2];
		R[0][1][IM]=wRSU2[1];

		R[1][0][RE]=-wRSU2[2];
		R[1][0][IM]=wRSU2[1];

		R[1][1][RE]=wRSU2[0];
		R[1][1][IM]=-wRSU2[3];


		ProdutoSU3(R,Omegadagger,Rw);
		ConjHermSU3(Rw,Omega);


		//	Segunda Submatriz

		wS[0]=(Omega[0][0][RE]+Omega[2][2][RE]);
		wS[1]=(Omega[0][2][IM]+Omega[2][0][IM]);
		wS[2]=(Omega[0][2][RE]-Omega[2][0][RE]);
		wS[3]=(Omega[0][0][IM]-Omega[2][2][IM]);


		ProjecaoSU2(wS,wSSU2);

		S[0][0][RE]=wSSU2[0];
		S[0][0][IM]=wSSU2[3];

		S[0][2][RE]=wSSU2[2];
		S[0][2][IM]=wSSU2[1];

		S[2][0][RE]=-wSSU2[2];
		S[2][0][IM]=wSSU2[1];

		S[2][2][RE]=wSSU2[0];
		S[2][2][IM]=-wSSU2[3];


	    ProdutoSU3(S,Rw,SRw);
	    ConjHermSU3(SRw,Omega);


	    //	Terceira Submatriz

		wT[0]=(Omega[1][1][RE]+Omega[2][2][RE]);
		wT[1]=(Omega[1][2][IM]+Omega[2][1][IM]);
		wT[2]=(Omega[1][2][RE]-Omega[2][1][RE]);
		wT[3]=(Omega[1][1][IM]-Omega[2][2][IM]);


		ProjecaoSU2(wT,wTSU2);

		T[1][1][RE]=wTSU2[0];
		T[1][1][IM]=wTSU2[3];

		T[1][2][RE]=wTSU2[2];
		T[1][2][IM]=wTSU2[1];

		T[2][1][RE]=-wTSU2[2];
		T[2][1][IM]=wTSU2[1];

		T[2][2][RE]=wTSU2[0];
		T[2][2][IM]=-wTSU2[3];

		ProdutoSU3(T,SRw,TSRw);
		ConjHermSU3(TSRw,Omega);


	}

	Inversa3por3(w,winv);
	ProdutoSU3(TSRw,winv,A);	//	Matriz de atualização é a matriz acumulada no produto de TSR nos hits
								//	multiplica por winv para eliminar w do produto

}

void LosAlamosSU3(int posicao[d],double w[3][3][2]){	

//	Generalização do algoritmo descrito em hep-lat/0301019v2, através do uso de submatrizes pelo truque de Cabbibo-Marinari
//	Recebe como entrada a posição de g(x), que será alterado para que w(x) seja minimizado

	double A[3][3][2];	//	Atualização final, após todos os hits
	double gatualizado[3][3][2];	//	Transformação de calibre anterior

	BlocoComumLASU3(w,A);

	ProdutoSU3(A,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],gatualizado);
	ProjecaoSU3(gatualizado,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]]);

}

void SobreRelaxacaoSU3(int posicao[d],double w[3][3][2]){	

//	Generalização do algoritmo descrito em hep-lat/0301019v2, através do uso de submatrizes pelo truque de Cabbibo-Marinari
//	Recebe como entrada a posição de g(x), que será alterado para que w(x) seja minimizado

	double A[3][3][2];	//	Atualização final, após todos os hits

	double fator1[2]={1.0-omegaSR,0.0};
	double fator2[2]={omegaSR,0.0};

	double omegaSRA[3][3][2];
	double ASR[3][3][2];
	double gatualizado[3][3][2];

	BlocoComumLASU3(w,A);

	//	Após determinação de T, S e R, atualização de g(x)

	MultiplicacaoEscalarSU3(identidadeSU3,fator1,ASR);
	MultiplicacaoEscalarSU3(A,fator2,omegaSRA);

	SomaSU3(omegaSRA,ASR);
	ProdutoSU3(ASR,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],gatualizado);
	ProjecaoSU3(gatualizado,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]]);
	
}


void SobreRelaxacaoEstocasticaSU3(int posicao[d],double w[3][3][2]){	

//	Generalização do algoritmo descrito em hep-lat/0301019v2, através do uso de submatrizes pelo truque de Cabbibo-Marinari
//	Recebe como entrada a posição de g(x), que será alterado para que g(x).w(x) seja minimizado

	double aleatorio[1];	//	Variável para sorteio aleatório

	double A[3][3][2];	//	Atualização final, após todos os hits
	double gatualizado[3][3][2];	//	g atualizado a projetar

	BlocoComumLASU3(w,A);

	ranlxd(aleatorio,1);
	
	if(aleatorio[0]>pSRE){
		ProdutoSU3(A,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],gatualizado);
		ProjecaoSU3(gatualizado,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]]);
	}	
	else{
		ProdutoSU3Tres(A,A,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],gatualizado);
		ProjecaoSU3(gatualizado,G[posicao[0]][posicao[1]][posicao[2]][posicao[3]]);
	}

}