#ifndef ESTRUTURAS_H_INCLUDED
#define ESTRUTURAS_H_INCLUDED

#define MAX_FAC 852
#define MAX_CLI 2102

typedef struct solucao{
    int matriz[MAX_CLI][4]; //Matriz de quádruplas <f1,f2,q1,q2> de cada cliente
    int fo;
} Solucao;

int num_fac;
int num_cli;
int num_incompat;

int vet_capacidade[MAX_FAC]; //Capacidade máxima de uma facilidade
int vet_custo_fixo[MAX_FAC]; //Custo de abertura de Facilidade
int vet_demanda[MAX_CLI]; //Demanda total de cada cliente
int vet_aux_incompat[MAX_CLI]; //Vetor auxiliar da mat_incompat
int vet_aux_fac_pref[MAX_CLI]; //Vetor auxiliar da matriz fac_pref
int vet_aux_cli_pref[MAX_FAC]; //Vetor auxiliar da matriz cli_pref
int vet_aux_fac_viz[MAX_FAC] = {0}; //Vetor auxiliar da matriz de facilidades vizinhas

int mat_custo[MAX_CLI][MAX_FAC]; //Matriz de custo unitário [Cliente x Facilidade]
int mat_incompat[MAX_CLI][MAX_CLI]; //Matriz incompatibilidade [Cliente x Cliente]
int mat_bin_incompat[MAX_CLI][MAX_CLI]; //Matriz Binária de incompatibilidade [Cliente x Cliente]
int mat_fac_pref[MAX_CLI][MAX_FAC]; //Matriz de preferência de facilidade [Cliente x Facilidade]
int mat_cli_pref[MAX_FAC][MAX_CLI]; //Matriz de preferência de cliente [Facilidade x Cliente]
int mat_fac_viz[MAX_FAC][MAX_FAC]; //Matriz de facilidades vizinhas [Facilidade x Facilidade]



// ----- FUNÇÕES DE UTILIDADE -----
void testar_dados(char* arq);
void ler_dados(char* arq, int& num_pref);
void criar_mat_pref(int& num_pref);
void criar_mat_viz();
void quicksortCentro(int* vetor,int* vetor2,int inicio, int fim);
//----------------------------------

// ----- FUNÇÕES BÁSICAS -----
void calcular_fo(Solucao &s);

#endif // ESTRUTURAS_H_INCLUDED