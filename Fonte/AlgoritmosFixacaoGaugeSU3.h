#ifndef ALGORITMOSFIXACAOGAUGESU3_H
#define ALGORITMOSFIXACAOGAUGESU3_H

void MaximoDecliveSU3(int posicao[d]);

void BlocoComumLASU3(double w[3][3][2],double A[3][3][2]);

void LosAlamosSU3(int posicao[d],double h[3][3][2]);

void SobreRelaxacaoSU3(int posicao[d],double w[3][3][2]);

void SobreRelaxacaoEstocasticaSU3(int posicao[d],double w[3][3][2]);

#endif