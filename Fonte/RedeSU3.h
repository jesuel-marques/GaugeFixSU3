#ifndef REDESU3_H
#define REDESU3_H

void InicializarUSU3(int Partida);

void CarregarUSU3(char NomeArquivo[ComprimentoMaxArquivo]);

void ImprimirUauxSU3(char NomeArquivo[ComprimentoMaxArquivo]);

void CopiarPosicao(int v[d],int u[d]);

void SomaVetoresPosicao(int v[d], int u[d], int vmaisu[d]);

void PosicoesVizinhas(int posicao[d], int mu, int nu, int posicaomaismu[d], int posicaomaismumaisnu[d],int posicaomaisnu[d], int posicaomaismumenosnu[d], int posicaomenosnu[d]);

void EloVizinhoSU3(int posicao[d],int direcao,int sentido,double u[3][3][2]);

#endif