#ifndef SU2_H
#define SU2_H

void CopiarSU2(double u[4], double v[4]);

void MultiplicacaoEscalarSU2(double u[4],double num, double numvezesu[4]);

void SomaSU2(double a[4], double resultado[4]);
void DiferencaSU2(double a[4], double b[4], double amenosb[4]);

void ProdutoSU2(double u[4], double v[4], double uv[4]);
void ProdutoSU2Tres(double a[4], double b[4], double c[4], double abc[4]);
void ProdutoSU2Quatro(double a[4], double b[4], double c[4], double e[4], double abce[4]);

void PotenciaSU2(double u[4],double unapotencias[10][4]);

double TrSU2(double u[4]);
double DeterminanteSU2(double u[4]);

void ConjHermSU2(double u[4],double udag[4]);

void ProjecaoSU2(double a[4],double aSU2[4]);

#endif