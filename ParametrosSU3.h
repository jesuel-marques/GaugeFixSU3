//Parâmetros da simulação

#define Nxyz 8 //Lado da rede na direção espacial
#define Nt 8 //Lado da rede na direção temporal
#define d 4	// Dimensão espaço temporal da rede

#define Volume Nxyz*Nxyz*Nxyz*Nt //	Número de pontos da rede

#define Nc 3.0	//	Número de cores da teoria de gauge

// Paramêtros gerais da Simulação

#define MAXconfigs 4 //	Número máximo de configurações descorrelacionadas a serem geradas

#define ComprimentoMaxArquivo 2000

//	Códigos

#define FRIA 0	//	Possíveis Partidas para Inicialização da Rede
#define QUENTE 1

#define TRAS -1	//	Usados para se referir a posições vizinhas na rede 
#define FRENTE 1

#define RE 0	// 	Usado para se referir a parte real e imaginária
#define IM 1

#define NAO 0	//	Usado como flags para testes
#define SIM 1
