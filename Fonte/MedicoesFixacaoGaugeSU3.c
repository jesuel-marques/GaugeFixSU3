#include <stdio.h>								//	Arquivos de cabeçalho padrão C

#include "../ParametrosSU3.h"					//	Parâmetros da simulação
#include "../ParametrosFixacaoGaugeSU3.h"		//	Parâmetros dos algoritmos
#include "GlobalSU3.h"							//	Definição de variáveis globais
#include "GlobalFixacaoGaugeSU3.h"				//	Definição de variáveis globais específicas da fixação

#include "FuncoesMatematicas.h"					//	Funções de cálculo SU(3) e gerais de rede
#include "SU3.h"								//	Funções de cálculo de SU3
#include "RedeSU3.h"							//	Funções de inicialização, cálculo de posições na rede e vizinhos
#include "MedicoesSU3.h"						//	Funções de medição de quantidades como A_mu(n) e a div A_mu(n)

double Calculare2SU3(){

//	Função calcula e2 (definido na Dissertação Nelson, eq. 4.53, e hep-lat/0301019v2), 
//	utilizado para definir distancia até situação de gauge fixo.
	
	int posicao[d];

	double divA[3][3][2];	//	Quadridivergência de A
	double e2=0.0;	

	double componentesdivA[9];	//	Componentes da quadridivergência projetada nas matrizes de Gell-Mann

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
		for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
			for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
				for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++){

					//	Calcula quadridivergencia e decompõe na álgebra de SU(3)

					DivergenciaASU3(posicao,divA);
					DecomporAlgebraSU3(divA,componentesdivA);
					for(int a=1;a<=8;a++)	//	Soma normalizada dos quadrados das componentes de cor da divergência de A
						e2+=(double)pow2(componentesdivA[a])/(Volume);

				}

	return e2;					
}

double Calculare6SU3(){	

//	Função calcula e6 (definido em hep-lat/0301019v2), 
//	utilizado para definir distancia até situação de gauge fixo.
	
	int posicao[d];
	double A[3][3][2];
	double componentesA[9];
	double Qnu[9][Nxyz];
	double Qnumedia[9];

	double parcial;
	double e6=0.0;

	for(int nu=0;nu<d;nu++){
		

		for(int j=0;j<=8;j++){	//	Inicializando Qnu e Qnumedia
			Qnumedia[j]=0.0;
			for(posicao[nu]=0;posicao[nu]<Nxyz;posicao[nu]++)
				Qnu[j][posicao[nu]]=0.0;	
		}
	
		for(posicao[nu]=0;posicao[nu]<Nxyz;posicao[nu]++){ //Calculo de Qnu[j][x[nu]]

			for(posicao[(nu+1)%d]=0;posicao[(nu+1)%d]<Nxyz;posicao[(nu+1)%d]++)
				for(posicao[(nu+2)%d]=0;posicao[(nu+2)%d]<Nxyz;posicao[(nu+2)%d]++)
					for(posicao[(nu+3)%d]=0;posicao[(nu+3)%d]<Nxyz;posicao[(nu+3)%d]++){
				
						CalcularASU3(posicao,nu,A);
						DecomporAlgebraSU3(A,componentesA);

						for(int j=1;j<=8;j++)
							Qnu[j][posicao[nu]]+=componentesA[j];

					}

					
			for(int j=1;j<=8;j++){	//Calculo de Qnumedia[j]
				if(nu!=0)
					Qnumedia[j]+=(double)Qnu[j][posicao[nu]]/Nxyz;
				else
					Qnumedia[j]+=(double)Qnu[j][posicao[nu]]/Nt;
			}
			
		}

		parcial=0.0;

		for(int j=1;j<=8;j++)
			if(nu!=0){
				for(posicao[nu]=0;posicao[nu]<Nxyz;posicao[nu]++)
					parcial+=pow2(Qnu[j][posicao[nu]]-Qnumedia[j])/pow2(Qnumedia[j]);
			}
			else{
				for(posicao[nu]=0;posicao[nu]<Nt;posicao[nu]++)
					parcial+=pow2(Qnu[j][posicao[nu]]-Qnumedia[j])/pow2(Qnumedia[j]);
			}

		if(nu!=0){
			e6+=(double)parcial/Nxyz;
		}
		else{
			e6+=(double)parcial/Nt;
		}

	
	}

	return (double)e6/(8.0*d);


}

double CalcularESU3(){ 

	//	Funcional a ser minimizado na condição de calibre de Landau (eqs. 4.39 e 4.49 da Dissertação Nelson).

	double E=0.0;
	double traco[2];

	double termo[3][3][2];


	int posicao[d];	//	Posições relevantes para o cálculo
	int versormaismu[d]={0,0,0,0};
	int posicaomaismu[d];

	double gfrentedagger[3][3][2];

	double prod[3][3][2];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
		for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
			for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
				for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
					for(int mu=0;mu<d;mu++){

						versormaismu[mu]=1;

						SomaVetoresPosicao(posicao,versormaismu,posicaomaismu);
						
						ConjHermSU3(G[posicaomaismu[0]][posicaomaismu[1]][posicaomaismu[2]][posicaomaismu[3]],gfrentedagger);

						ProdutoSU3Tres(G[posicao[0]][posicao[1]][posicao[2]][posicao[3]],U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],gfrentedagger,prod);
						DiferencaSU3(identidadeSU3,prod,termo);
						TrSU3(termo,traco);
						E+=traco[RE]/(4.0*3.0*Volume);

						versormaismu[mu]=0;

					}

	return E;
}