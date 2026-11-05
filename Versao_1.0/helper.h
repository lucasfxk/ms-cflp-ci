#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#define MAX_FAC 852
#define MAX_CLI 2102

typedef struct solucao{
    int fo; //Função Objetivo
    int mat_cli[MAX_CLI][4]; //Matriz de quádruplas <f1,f2,q1,q2> de cada cliente
    int mat_num_incompat[MAX_CLI][MAX_FAC]; //Matriz quantidade de incompatibilidade de cada cliente para cada facilidade
    int vet_fac_qnt[MAX_FAC]; //Vetor quantidade de demanda atendida por cada facilidade
    int mat_fac_cli[MAX_FAC][MAX_CLI]; // Matriz da lista clientes de cada facilidade
    int vet_aux_fac_cli[MAX_FAC];
} Solucao;

int num_fac;
int num_cli;
int num_incompat;
int vet_capacidade[MAX_FAC]; //Capacidade máxima de uma facilidade
int vet_custo_fixo[MAX_FAC];
int vet_demanda[MAX_CLI]; //Demanda total de cada cliente
int vet_aux_pref[MAX_CLI]; //Vetor auxiliar da matriz fac_pref
int vet_aux_cli_pref[MAX_FAC]; //Vetor auxiliar da matriz cli_pref
int vet_aux_incompat[MAX_CLI];
int mat_incompat[MAX_CLI][MAX_CLI]; //Matriz incompatibilidade [Cliente x Cliente]
int mat_bin_incompat[MAX_CLI][MAX_CLI]; //Matriz Binária de incompatibilidade [Cliente x Cliente]
int mat_custo[MAX_CLI][MAX_FAC]; //Matriz de custo unitário [Cliente x Facilidade]
int mat_fac_pref[MAX_CLI][MAX_FAC]; //Matriz de preferência de facilidade [Cliente x Facilidade]
int mat_cli_pref[MAX_FAC][MAX_CLI]; //Matriz de preferência de cliente [Facilidade x Cliente]

//Pré-Processamento
void quicksortCentro(int* vetor,int* vetor2,int inicio, int fim);
void ler_dados();
void testar_dados(char* arq);
void escrever_solucao(Solucao& s);
void criar_mat_pref();
void calcular_fo(Solucao& s);

//Construtivas
void heu_const_ale(Solucao& s);
void multi_start_greedy_const(Solucao& s);

//Funções de Gerar Vizinho
void changeSupplier(Solucao& s);
void swapSuppliers(Solucao& s);
void clopenFacilities_autoral_v1(Solucao& s);
void gerarVizinho(Solucao& s);
void gerarVizinho_v2(Solucao &s,const double& TempAtual,const int& melhorFo,const double& p);

//Meta-Heurística
void simulated_annealing(Solucao &s);
void cs_gsa(Solucao &s);

//Buscas Locas
void bl_pm_swap(Solucao& s);
void bl_mm_swap(Solucao& s);
void bl_pm_change(Solucao& s);

//CPLEX
void gerarArquivoCplex();
void otimizarCplex();

#endif // HELPER_H_INCLUDED
