#ifndef CONSTRUTIVAS_H_INCLUDED
#define CONSTRUTIVAS_H_INCLUDED
#include "estruturas.h"

const double fatorAmortizacao = 0.25;
const double janelaDeAceitacao = 0.288;

typedef struct tripla{
    int cliente;
    int facilidade;
    int carga;
} Tripla; // <i, j, q> onde i é o cliente, j a facilidade e q a quantidade de carga.


void construtivaAleatoriaGulosa(Solucao&);

void calcularCusto(int& cliente, int& facilidade, int& quantidade, int& custo);
void criarVetorDeIndex(int vetor[],int& tamanho);
void eCompativel(Solucao& s, int& cliente, int& facilidade, bool& compativel);
void montarMatrizSolucao(Solucao& s);
void inserirAbastecimento(Solucao& s, int& cliente, int& facilidade, int& q);

#endif // CONSTRUTIVAS_H_INCLUDED