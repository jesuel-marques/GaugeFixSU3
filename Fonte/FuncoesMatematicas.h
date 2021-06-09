#ifndef FUNCOESMATEMATICAS_H
#define FUNCOESMATEMATICAS_H

double Sign(double  x);

double pow2(double x);

double  ProdutoInterno(double  x[4],double  u[4]);

void ProdutoExterno(double  x[4], double  u[4], double  prodext[4]);

void CopiarComplexo(double x[2],double xcopiado[2]);

void ConjugadoComplexo(double x[2],double xestrela[2]);

void SomaComplexo(double x[2],double soma[2]);

void DiferencaComplexo(double x[2],double y[2],double diferenca[2]);

void ProdutoComplexo(double x[2], double y[2],double xy[2]);

void ProdutoComplexoTres(double x[2],double y[2], double w[2],double xyw[2]);

void ProdutoComplexoQuatro(double x[2],double y[2], double w[2], double z[2], double xywz[2]);

double ModuloQuadComplexo(double x[2]);

void DivisaoComplexo(double x[2], double y[2], double xpory[2]);

void Inversa3por3(double a[3][3][2],double ainv[3][3][2]);
#endif