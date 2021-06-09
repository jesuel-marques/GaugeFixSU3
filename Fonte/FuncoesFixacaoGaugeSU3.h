#ifndef FUNCOESFIXAGAUGESU3_H
#define FUNCOESFIXAGAUGESU3_H


void CalcularwSU3(int posicao[d],double w[3][3][2]);

void AtualizarUauxSU3();
void AtualizarUauxLocalSU3(int posicao[d]);

int FixarGaugeSU3(int config,int Algoritmo, double tolerancia);

void ImprimirUGaugeFixadoSU3(char NomeArquivo[ComprimentoMaxArquivo]);

#endif