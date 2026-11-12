#include <iostream>
#include <stdio.h>
//#include <stdlib.h>
#include <climits>
#include <memory.h>
#include <math.h>
#include <time.h>
#include "helper.h"

#define MAX(X,Y) ((X > Y) ? X : Y)
#define MIN(X,Y) ((X < Y) ? X : Y)
#define BETA 1000
#define ALPHA 10000

//#define DBG_SWAPSUPPLIERS
#define DBG_SA
#define DBG_BL_PM_SWAP
#define DBG_BL_MM_SWAP
#define DBG_BL_PM_CHANGE
#define DBG_GSA

#define CONST_ALE_MODO_1

using namespace std;

//const char arq[] = "cflp-ci_toy.txt";
char arq[100] /*= "wlp01.txt"*/;
//const char arq[] = "wlp01.txt";

const int beta = 1; //Tamanho da lista de preferência
const double SS_bias = 0.28; //No artigo está 0.44 porém considerando a implementação fica 0.28

//Parâmetros Construtiva Gulosa
const double equal_tolerance = 0.288; 
const double greedy_amortize_factor = 0.25; //Parâmetro α (fator de amortização do custo fixo)

/*const double prob_open = 0.16;
const double prob_close = 0.02;
const double prob_open_close = 0.82;*/
//Probabilidade ClopenFacilities Ajustado:
const double prob_open = 0.98;
const double prob_close = 0.84;
const int lambda = /*8*/100 ;

//Melhor FO: 30614 (30828)
//Simulated_annealing constantes
const float a = 0.994; //0.994
const float Tzero = 16.42;
const float Tc = 0.183;
int SAmax;
const int r_max = 10; //Número máximo de reaquecimentos
const int r_i = 15000; //Número de iterações sem melhoria
const double TempChange = 0.500; //Temperatura do Reanealing
const double r_f = Tzero; //Fator de reaquecimento
const double q_a = 0.13; //1- SA normal 2- Fast SA

//Probabilidade de escolha gerarVizinho
const double p_SS = 0.19; //Probabilidade SwapSuppliers
const double p_CF = 0.01; //Probabilidade ClopenFacilities
const double p_CS = 1 - (p_SS + p_CF); //Proababilidade de ChangeSuppliers

//GerarVizinho v2
const double lambda_i = 5;
const double lambda_f = 1;
const int lim_super = 5; //Limite de passos

//GerarVizinho GSA
double qV = 5.0;
double qA = 1.25;
double qT = 2.0;
double T0 = 10.0;
/*
qA: Controla a aceitação:
	maior → aceita mais pioras;
	1 → SA clássico.
    2 → FSA.

qT: Controla o resfriamento:
	maior → temperatura cai lentamente;
	menor → converge mais rápido.

qV: Controla quantos bits são perturbados:
	próximo de 1 → poucos flips;
	maior → movimentos longos.
*/

int num_pref;

int main(int argc, char* argv[])
{
    for(int x = 0; x < argc; x++){
        if(argv[x] == "--qV"){
            qV = atof(argv[x]);
        }else if(argv[x] == "--qA"){
            qA = atof(argv[x]);
        }else if(argv[x] == "--qT"){
            qT = atof(argv[x]);
        }else{
            strcpy(arq,argv[x]);
        }
    }
    //srand(time(NULL));
    ler_dados();
    criar_mat_pref();
    criar_mat_viz();
    SAmax = (num_cli * num_fac)/180;

    /*for(int x=0; x<num_fac; x++){
        printf("\nFac %i: ", x);
        for(int y=0; y<vet_aux_fac_viz[x]; y++){
            printf("%i ",mat_fac_viz[x][y]);
        }
    }*/
    //testar_dados("");

    Solucao s;
    multi_start_greedy_const(s);
    //heu_const_ale(s);
    printf("FO: %i\n",s.fo);
    //escrever_solucao(s);

    //printf("\n\n");
    GSA(s, 10 * sqrt(num_fac));
    //simulated_annealing(s, 10 * sqrt(num_fac));
    //bl_pm_change(s);
    //bl_pm_swap(s);

    printf("FO: %i",s.fo);
    //escrever_solucao(s);

    return 0;
}

void bl_pm_change(Solucao& s){
    int cli, index_cli, index_fac;
    int antiga_fac,nova_fac; //Facilidade antiga e nova
    int vet_idx_cli[num_cli]; //Vetor de index de clientes (aleatoriedade sem repetição)
    for(int x=0; x<num_cli; x++) vet_idx_cli[x] = x;
    int vet_idx_fac[num_fac]; //Vetor de index de facilidades (aleatoriedade sem repetição)
    int cap_resid_a, cap_resid_b, n_carga;
    int parar = 0;
    Solucao temp;

    INICIO_BL_PM_CHANGE : ;
    memcpy(&temp, &s, sizeof(Solucao));

    for(int x=0; x<num_cli; x++){
        index_cli =  (rand() % (num_cli - x)) + x;
        cli = vet_idx_cli[index_cli];
        vet_idx_cli[index_cli] = vet_idx_cli[x];
        vet_idx_cli[x] = cli;

        //Copia as facilidades preferenciais para o vetor auxiliar de index
        for(int y=0; y<vet_aux_pref[cli]; y++) vet_idx_fac[y] = mat_fac_pref[cli][y];
        
        if(temp.mat_cli[cli][1] != -1){ //Possui 2 facilidades
            antiga_fac = temp.mat_cli[cli][1]; /*---LEMBRAR DE TROCAR TODOS---*/
            int base = 0;

            for(int base=0; base<vet_aux_pref[cli]; base++){
                index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                nova_fac = vet_idx_fac[index_fac];
                vet_idx_fac[index_fac] = vet_idx_fac[base];
                vet_idx_fac[base] = nova_fac;
                cap_resid_a = MAX(0,vet_capacidade[temp.mat_cli[cli][0]] - temp.vet_fac_qnt[temp.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                cap_resid_b = MAX(0,vet_capacidade[nova_fac] - temp.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade

                if(nova_fac == temp.mat_cli[cli][0] 
                || nova_fac == temp.mat_cli[cli][1]
                //|| mat_custo[cli][nova_fac] >= mat_custo[cli][s.mat_cli[cli][1]] //Talvez tenha que tirar*
                || temp.mat_num_incompat[cli][nova_fac] != 0
                || cap_resid_b <= 0 
                || (cap_resid_a + cap_resid_b) < temp.mat_cli[cli][3]  //Verifica se a carga atual do 2º fornecedor
                ) continue;                                         //cabe na capacidade residual de f_1 e nova_fac

                if(mat_custo[cli][temp.mat_cli[cli][0]] <= mat_custo[cli][nova_fac]){ //Caso f_1 seja mais barato que nova_fac
                    if(temp.mat_cli[cli][3] <= cap_resid_a){ //Caso caiba tudo em f_1
                        temp.mat_cli[cli][2] += temp.mat_cli[cli][3]; //Passa tudo para f_1
                        temp.vet_fac_qnt[temp.mat_cli[cli][0]] += temp.mat_cli[cli][3];
                        temp.vet_fac_qnt[temp.mat_cli[cli][1]] -= temp.mat_cli[cli][3];

                        //Remove cliente da matriz de lista de clientes por facilidade

                        for(int y=0; y<temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]; y++){
                            if(temp.mat_fac_cli[temp.mat_cli[cli][1]][y] == cli){
                                //Troca o elemento a ser removido com o último e decrementa o contadora
                                temp.mat_fac_cli[temp.mat_cli[cli][1]][y] = temp.mat_fac_cli[temp.mat_cli[cli][1]][temp.vet_aux_fac_cli[temp.mat_cli[cli][1]] - 1];
                                temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]--;
                                break;
                            }
                        }

                        //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                        for(int y=0; y<vet_aux_incompat[cli]; y++) temp.mat_num_incompat[mat_incompat[cli][y]][temp.mat_cli[cli][1]] -= 1;
                        
                        temp.mat_cli[cli][1] = -1; //remove f_2
                        temp.mat_cli[cli][3] = 0;

                    }else{//Caso NÃO caiba tudo em f_1

                        //Quantidade de carga para nova_fac
                        n_carga = temp.mat_cli[cli][3] - cap_resid_a;
                        temp.mat_cli[cli][2] += cap_resid_a;
                        temp.vet_fac_qnt[temp.mat_cli[cli][0]] += cap_resid_a;
                        temp.vet_fac_qnt[temp.mat_cli[cli][1]] -= temp.mat_cli[cli][3];

                        //Remove cliente da matriz de lista de clientes por facilidade
                        for(int y=0; y<temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]; y++){
                            if(temp.mat_fac_cli[temp.mat_cli[cli][1]][y] == cli){
                                //Troca o elemento a ser removido com o último e decrementa o contadora
                                temp.mat_fac_cli[temp.mat_cli[cli][1]][y] = temp.mat_fac_cli[temp.mat_cli[cli][1]][temp.vet_aux_fac_cli[temp.mat_cli[cli][1]] - 1];
                                temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]--;
                                break;
                            }
                        }

                        //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                        for(int y=0; y<vet_aux_incompat[cli]; y++) temp.mat_num_incompat[mat_incompat[cli][y]][temp.mat_cli[cli][1]] -= 1;
                        
                        temp.mat_cli[cli][1] = nova_fac;
                        temp.mat_cli[cli][3] = n_carga;
                        temp.vet_fac_qnt[nova_fac] +=  n_carga;
                        temp.mat_fac_cli[nova_fac][temp.vet_aux_fac_cli[nova_fac]] = cli;
                        temp.vet_aux_fac_cli[nova_fac]++;

                        for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                            temp.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                        }
                    }

                }else{ //Caso nova_fac seja mais barato 

                    //Define a carga de nova_fac
                    n_carga = MIN(temp.mat_cli[cli][3],cap_resid_b);
                    temp.mat_cli[cli][2] += temp.mat_cli[cli][3] - n_carga; //Atribui o restante à f_1
                    temp.vet_fac_qnt[temp.mat_cli[cli][0]] += temp.mat_cli[cli][3] - n_carga;
                    temp.vet_fac_qnt[temp.mat_cli[cli][1]] -= temp.mat_cli[cli][3];

                    //Remove cliente da matriz de lista de clientes por facilidade
                    for(int y=0; y<temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]; y++){
                        if(temp.mat_fac_cli[temp.mat_cli[cli][1]][y] == cli){
                            //Troca o elemento a ser removido com o último e decrementa o contadora
                            temp.mat_fac_cli[temp.mat_cli[cli][1]][y] = temp.mat_fac_cli[temp.mat_cli[cli][1]][temp.vet_aux_fac_cli[temp.mat_cli[cli][1]] - 1];
                            temp.vet_aux_fac_cli[temp.mat_cli[cli][1]]--;
                            break;
                        }
                    }

                    //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                    for(int y=0; y<vet_aux_incompat[cli]; y++) temp.mat_num_incompat[mat_incompat[cli][y]][temp.mat_cli[cli][1]] -= 1;
                    
                    temp.mat_cli[cli][1] = nova_fac;
                    temp.mat_cli[cli][3] = n_carga;
                    temp.vet_fac_qnt[nova_fac] +=  n_carga;
                    temp.mat_fac_cli[nova_fac][temp.vet_aux_fac_cli[nova_fac]] = cli;
                    temp.vet_aux_fac_cli[nova_fac]++;

                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        temp.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }

                    //Transforma nova_fac em f_1
                    int aux;
                    aux = temp.mat_cli[cli][0];
                    temp.mat_cli[cli][0] = temp.mat_cli[cli][1];
                    temp.mat_cli[cli][1] = aux;

                    aux = temp.mat_cli[cli][2];
                    temp.mat_cli[cli][2] = temp.mat_cli[cli][3];
                    temp.mat_cli[cli][3] = aux;
                }

                break;
            }

        }else{ //Possui 1 facilidade
            antiga_fac = rand() % 2;

            if(antiga_fac == 0){ //Caso seja f_a
                
                for(int base=0; base < vet_aux_pref[cli]; base++){
                    index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                    nova_fac = vet_idx_fac[index_fac];
                    vet_idx_fac[index_fac] = vet_idx_fac[base];
                    vet_idx_fac[base] = nova_fac;

                    cap_resid_a = MAX(0,vet_capacidade[temp.mat_cli[cli][0]] - temp.vet_fac_qnt[temp.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                    cap_resid_b = MAX(0,vet_capacidade[nova_fac] - temp.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade


                    if(nova_fac == temp.mat_cli[cli][0] //Verifica qualquer tipo de incompatibilidade
                    || temp.mat_cli[cli][2] > cap_resid_b
                    || temp.mat_num_incompat[cli][nova_fac] != 0
                    ) continue;

                    temp.vet_fac_qnt[temp.mat_cli[cli][0]] -= temp.mat_cli[cli][2];

                    //Remove cliente da matriz de lista de clientes por facilidade
                    for(int y=0; y<temp.vet_aux_fac_cli[temp.mat_cli[cli][0]]; y++){
                        if(temp.mat_fac_cli[temp.mat_cli[cli][0]][y] == cli){
                            //Troca o elemento a ser removido com o último e decrementa o contadora
                            temp.mat_fac_cli[temp.mat_cli[cli][0]][y] = temp.mat_fac_cli[temp.mat_cli[cli][0]][temp.vet_aux_fac_cli[temp.mat_cli[cli][0]] - 1];
                            temp.vet_aux_fac_cli[temp.mat_cli[cli][0]]--;
                            break;
                        }
                    }

                    //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                    for(int y=0; y<vet_aux_incompat[cli]; y++) temp.mat_num_incompat[mat_incompat[cli][y]][temp.mat_cli[cli][0]] -= 1;
                    
                    
                    temp.mat_cli[cli][0] = nova_fac;
                    temp.vet_fac_qnt[nova_fac] += temp.mat_cli[cli][2];
                    temp.mat_fac_cli[nova_fac][temp.vet_aux_fac_cli[nova_fac]] = cli;
                    temp.vet_aux_fac_cli[nova_fac]++;
                    
                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        temp.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }
                    
                    break;
                }
            }else{ //Caso seja f_b (f_b == -1)

                for(int base=0; base < vet_aux_pref[cli]; base++){
                    index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                    nova_fac = vet_idx_fac[index_fac];
                    vet_idx_fac[index_fac] = vet_idx_fac[base];
                    vet_idx_fac[base] = nova_fac;

                    cap_resid_a = MAX(0,vet_capacidade[temp.mat_cli[cli][0]] - temp.vet_fac_qnt[temp.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                    cap_resid_b = MAX(0,vet_capacidade[nova_fac] - temp.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade
                    
                    if(nova_fac == temp.mat_cli[cli][0]
                    || mat_custo[cli][nova_fac] > mat_custo[cli][temp.mat_cli[cli][0]]
                    || temp.mat_num_incompat[cli][nova_fac] != 0
                    || cap_resid_b == 0
                    ) continue;
                    
                    n_carga = MIN(cap_resid_b, temp.mat_cli[cli][2] - 1); //Sempre deixa 1 carga restante pra f_a

                    if(n_carga == 0) continue; //Evita que entre um cliente que fique ocioso (carga 0)

                    //Retira carga de f_a
                    temp.mat_cli[cli][2] -= n_carga;
                    temp.vet_fac_qnt[temp.mat_cli[cli][0]] -= n_carga;

                    //Insere nova_fac
                    temp.mat_cli[cli][1] = nova_fac;
                    temp.mat_cli[cli][3] = n_carga;
                    temp.vet_fac_qnt[nova_fac] += n_carga;
                    temp.mat_fac_cli[nova_fac][temp.vet_aux_fac_cli[nova_fac]] = cli;
                    temp.vet_aux_fac_cli[nova_fac]++;
                    
                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        temp.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }
                    
                    //Transforma nova_fac em f_1
                    int aux;
                    aux = temp.mat_cli[cli][0];
                    temp.mat_cli[cli][0] = temp.mat_cli[cli][1];
                    temp.mat_cli[cli][1] = aux;

                    aux = temp.mat_cli[cli][2];
                    temp.mat_cli[cli][2] = temp.mat_cli[cli][3];
                    temp.mat_cli[cli][3] = aux;

                    break;
                }
            }
        }

        calcular_fo(temp);

        if(temp.fo < s.fo){
            #ifdef DBG_BL_PM_CHANGE
            printf("FO: %i\n",temp.fo);
            #endif
            memcpy(&s, &temp, sizeof(Solucao));
            goto INICIO_BL_PM_CHANGE;
        }else{
            memcpy(&temp, &s, sizeof(Solucao));
        }
    }
}

void bl_pm_swap(Solucao& s){
    int idx_cli_1, idx_cli_2;
    int cli_1, cli_2;
    int vet_aux_1[num_cli], vet_aux_2[num_cli]; //Vetores de índices (para aleatoriedade)
    int pos_1, pos_2; //Facilidade escolhida (0 ou 1)
    int fac_1, fac_2; //índice da Facilidade
    int outra_fac;
    int aux; //Variável auxiliar mesmo...
    int cap_resid_1, cap_resid_2; //Capaciade residual
    int delta, carga_1, carga_2;
    Solucao temp;

    for(int x=0; x<num_cli; x++){ //Preenche vetores auxiliares
        vet_aux_1[x] = x;
    }

    INICIO_BL_PM_SWAP : ;
    memcpy(&temp, &s, sizeof(Solucao));

    for(int index_1=0; index_1 < num_cli; index_1++){ //Escolhe cli_1 (de forma aleatória)
        idx_cli_1 = (rand() % (num_cli - index_1)) + index_1;
        cli_1 = vet_aux_1[idx_cli_1];
        vet_aux_1[idx_cli_1] = vet_aux_1[index_1];
        vet_aux_1[index_1] = cli_1;

        //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
        pos_1 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, temp.mat_cli[cli_1][1] + 1); //Confia que funciona ;D
        fac_1 = temp.mat_cli[cli_1][pos_1];

        //Calcula capacidade residual da facilidade escolhida
        cap_resid_1 = MAX(0,vet_capacidade[fac_1] - temp.vet_fac_qnt[fac_1]);
        
        //Copia os clientes preferenciais
        for(int x=0; x<vet_aux_cli_pref[fac_1]; x++){
            vet_aux_2[x] = mat_cli_pref[fac_1][x];
        }
        for(int index_2=0; index_2 < vet_aux_cli_pref[fac_1]; index_2++){ //Escolhe cli_2 (de forma aleatória)
            idx_cli_2 = (rand() % (vet_aux_cli_pref[fac_1] - index_2)) + index_2;
            cli_2 = vet_aux_2[idx_cli_2];
            vet_aux_2[idx_cli_2] = vet_aux_2[index_2];
            vet_aux_2[index_2] = cli_2;

            //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
            pos_2 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, temp.mat_cli[cli_2][1] + 1);
            fac_2 = temp.mat_cli[cli_2][pos_2];

            //Calcula capacidade residual da facilidade escolhida
            cap_resid_2 = MAX(0,vet_capacidade[fac_2] - temp.vet_fac_qnt[fac_2]);

            if(cli_2 == cli_1
            || fac_1 == fac_2
            || cap_resid_2 < temp.mat_cli[cli_1][pos_1 + 2] - temp.mat_cli[cli_2][pos_2 + 2]
            || cap_resid_1 < temp.mat_cli[cli_2][pos_2 + 2] - temp.mat_cli[cli_1][pos_1 + 2]
            ) continue;

            if(mat_bin_incompat[cli_1][cli_2]){ 
                if(temp.mat_num_incompat[cli_1][fac_2] > 1
                || temp.mat_num_incompat[cli_2][fac_1] > 1
            ) continue;
            }else{
                if(temp.mat_num_incompat[cli_1][fac_2] > 0
                || temp.mat_num_incompat[cli_2][fac_1] > 0
            ) continue;
            }

            //Calcula delta
            carga_1 = temp.mat_cli[cli_1][pos_1 + 2];
            carga_2 = temp.mat_cli[cli_2][pos_2 + 2];
            delta = carga_1 * (mat_custo[cli_1][fac_2] - mat_custo[cli_1][fac_1])
            + carga_2 * (mat_custo[cli_2][fac_1] - mat_custo[cli_2][fac_2]);

            //Realiza a troca
            //No Cliente 1:
            temp.vet_fac_qnt[fac_1] -= temp.mat_cli[cli_1][pos_1 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<temp.vet_aux_fac_cli[fac_1]; y++){
                if(temp.mat_fac_cli[fac_1][y] == cli_1){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    temp.mat_fac_cli[fac_1][y] = temp.mat_fac_cli[fac_1][temp.vet_aux_fac_cli[fac_1] - 1];
                    temp.vet_aux_fac_cli[fac_1]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_1]; y++) temp.mat_num_incompat[mat_incompat[cli_1][y]][fac_1] -= 1;
            
            //Adiciona a fac_2
            if(temp.mat_cli[cli_1][(pos_1 + 1) % 2] != fac_2){ //Verifica se cli_1 já possui fac_2
                temp.mat_cli[cli_1][pos_1] = fac_2;
                temp.vet_fac_qnt[fac_2] += temp.mat_cli[cli_1][pos_1 + 2];
                temp.mat_fac_cli[fac_2][temp.vet_aux_fac_cli[fac_2]] = cli_1;
                temp.vet_aux_fac_cli[fac_2]++;
    
                for(int y=0; y<vet_aux_incompat[cli_1]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    temp.mat_num_incompat[mat_incompat[cli_1][y]][fac_2]++;
                }
            }else{
                temp.mat_cli[cli_1][0] = fac_2;
                temp.vet_fac_qnt[fac_2] += temp.mat_cli[cli_1][pos_1 + 2];
                temp.mat_cli[cli_1][2] += temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = 0;
                temp.mat_cli[cli_1][1] = -1;
            }


            //No Cliente 2:
            temp.vet_fac_qnt[fac_2] -= temp.mat_cli[cli_2][pos_2 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<temp.vet_aux_fac_cli[fac_2]; y++){
                if(temp.mat_fac_cli[fac_2][y] == cli_2){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    temp.mat_fac_cli[fac_2][y] = temp.mat_fac_cli[fac_2][temp.vet_aux_fac_cli[fac_2] - 1];
                    temp.vet_aux_fac_cli[fac_2]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_2]; y++) temp.mat_num_incompat[mat_incompat[cli_2][y]][fac_2] -= 1;
            
            //Adiciona a fac_1
            if(temp.mat_cli[cli_2][(pos_2 + 1) % 2] != fac_1){ //Verifica se cli_2 já possui fac_1
                temp.mat_cli[cli_2][pos_2] = fac_1;
                temp.vet_fac_qnt[fac_1] += temp.mat_cli[cli_2][pos_2 + 2];
                temp.mat_fac_cli[fac_1][temp.vet_aux_fac_cli[fac_1]] = cli_2;
                temp.vet_aux_fac_cli[fac_1]++;
    
                for(int y=0; y<vet_aux_incompat[cli_2]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    temp.mat_num_incompat[mat_incompat[cli_2][y]][fac_1]++;
                }
            }else{
                temp.mat_cli[cli_2][0] = fac_1;
                temp.vet_fac_qnt[fac_1] += temp.mat_cli[cli_2][pos_2 + 2];
                temp.mat_cli[cli_2][2] += temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = 0;
                temp.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_1 e nova f_2 sejam iguais (junta ambas) 
            if(temp.mat_cli[cli_1][0] == temp.mat_cli[cli_1][1]){
                temp.mat_cli[cli_1][2] += temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = 0;
                temp.mat_cli[cli_1][1] = -1;
            }
            if(temp.mat_cli[cli_2][0] == temp.mat_cli[cli_2][1]){
                temp.mat_cli[cli_2][2] += temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = 0;
                temp.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_2 seja mais barata que f_1
            if(temp.mat_cli[cli_1][1] != -1
            && mat_custo[cli_1][temp.mat_cli[cli_1][1]] < mat_custo[cli_1][temp.mat_cli[cli_1][0]]){
                aux = temp.mat_cli[cli_1][0];
                temp.mat_cli[cli_1][0] = temp.mat_cli[cli_1][1];
                temp.mat_cli[cli_1][1] = aux;
                aux = temp.mat_cli[cli_1][2];
                temp.mat_cli[cli_1][2] = temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = aux;
            }
            if(temp.mat_cli[cli_2][1] != -1
            && mat_custo[cli_2][temp.mat_cli[cli_2][1]] < mat_custo[cli_2][temp.mat_cli[cli_2][0]]){
                aux = temp.mat_cli[cli_2][0];
                temp.mat_cli[cli_2][0] = temp.mat_cli[cli_2][1];
                temp.mat_cli[cli_2][1] = aux;
                aux = temp.mat_cli[cli_2][2];
                temp.mat_cli[cli_2][2] = temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = aux;
            }

            /*calcular_fo(temp);
            if(temp.fo < s.fo){
                #ifdef DBG_BL_PM_SWAP
                printf("FO: %i\n",temp.fo);
                #endif
                memcpy(&s, &temp, sizeof(Solucao));
                goto INICIO_BL_PM_SWAP;
            }else{
                memcpy(&temp, &s, sizeof(Solucao));
            }*/
           if(delta < 0){
                temp.fo += delta;
                #ifdef DBG_BL_PM_SWAP
                printf("FO: %i\n",temp.fo);
                #endif
                memcpy(&s, &temp, sizeof(Solucao));
                goto INICIO_BL_PM_SWAP;
            }else{
                memcpy(&temp, &s, sizeof(Solucao));
            }
        }
    }
}

void bl_mm_swap(Solucao& s){
    int idx_cli_1, idx_cli_2;
    int cli_1, cli_2;
    int vet_aux_1[num_cli], vet_aux_2[num_cli]; //Vetores de índices (para aleatoriedade)
    int pos_1, pos_2; //Facilidade escolhida (1 ou 2)
    int fac_1, fac_2; //índice da Facilidade
    int outra_fac;
    int aux; //Variável auxiliar mesmo...
    int cap_resid_1, cap_resid_2; //Capaciade residual
    Solucao temp,melhor;
    memcpy(&melhor, &s, sizeof(Solucao));
    int melhorou;

    for(int x=0; x<num_cli; x++){ //Preenche vetores auxiliares
        vet_aux_1[x] = x;
    }

    INICIO_BL_MM_SWAP : ;
    memcpy(&temp, &s, sizeof(Solucao));
    melhorou = 0;

    for(int index_1=0; index_1 < num_cli; index_1++){ //Escolhe cli_1 (de forma aleatória)
        idx_cli_1 = (rand() % (num_cli - index_1)) + index_1;
        cli_1 = vet_aux_1[idx_cli_1];
        vet_aux_1[idx_cli_1] = vet_aux_1[index_1];
        vet_aux_1[index_1] = cli_1;

        //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
        pos_1 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, temp.mat_cli[cli_1][1] + 1); //Confia que funciona ;D
        fac_1 = temp.mat_cli[cli_1][pos_1];

        //Calcula capacidade residual da facilidade escolhida
        cap_resid_1 = MAX(0,vet_capacidade[fac_1] - temp.vet_fac_qnt[fac_1]);
        
        //Copia os clientes preferenciais
        for(int x=0; x<vet_aux_cli_pref[fac_1]; x++){
            vet_aux_2[x] = mat_cli_pref[fac_1][x];
        }
        for(int index_2=0; index_2 < vet_aux_cli_pref[fac_1]; index_2++){ //Escolhe cli_2 (de forma aleatória)
            idx_cli_2 = (rand() % (vet_aux_cli_pref[fac_1] - index_2)) + index_2;
            cli_2 = vet_aux_2[idx_cli_2];
            vet_aux_2[idx_cli_2] = vet_aux_2[index_2];
            vet_aux_2[index_2] = cli_2;

            //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
            pos_2 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, temp.mat_cli[cli_2][1] + 1);
            fac_2 = temp.mat_cli[cli_2][pos_2];

            //Calcula capacidade residual da facilidade escolhida
            cap_resid_2 = MAX(0,vet_capacidade[fac_2] - temp.vet_fac_qnt[fac_2]);

            if(cli_2 == cli_1
            || fac_1 == fac_2
            || cap_resid_2 < temp.mat_cli[cli_1][pos_1 + 2] - temp.mat_cli[cli_2][pos_2 + 2]
            || cap_resid_1 < temp.mat_cli[cli_2][pos_2 + 2] - temp.mat_cli[cli_1][pos_1 + 2]
            ) continue;

            if(mat_bin_incompat[cli_1][cli_2]){ 
                if(temp.mat_num_incompat[cli_1][fac_2] > 1
                || temp.mat_num_incompat[cli_2][fac_1] > 1
            ) continue;
            }else{
                if(temp.mat_num_incompat[cli_1][fac_2] > 0
                || temp.mat_num_incompat[cli_2][fac_1] > 0
            ) continue;
            }

            //Realiza a troca
            //No Cliente 1:
            temp.vet_fac_qnt[fac_1] -= temp.mat_cli[cli_1][pos_1 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<temp.vet_aux_fac_cli[fac_1]; y++){
                if(temp.mat_fac_cli[fac_1][y] == cli_1){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    temp.mat_fac_cli[fac_1][y] = temp.mat_fac_cli[fac_1][temp.vet_aux_fac_cli[fac_1] - 1];
                    temp.vet_aux_fac_cli[fac_1]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_1]; y++) temp.mat_num_incompat[mat_incompat[cli_1][y]][fac_1] -= 1;
            
            //Adiciona a fac_2
            if(temp.mat_cli[cli_1][(pos_1 + 1) % 2] != fac_2){ //Verifica se cli_1 já possui fac_2
                temp.mat_cli[cli_1][pos_1] = fac_2;
                temp.vet_fac_qnt[fac_2] += temp.mat_cli[cli_1][pos_1 + 2];
                temp.mat_fac_cli[fac_2][temp.vet_aux_fac_cli[fac_2]] = cli_1;
                temp.vet_aux_fac_cli[fac_2]++;
    
                for(int y=0; y<vet_aux_incompat[cli_1]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    temp.mat_num_incompat[mat_incompat[cli_1][y]][fac_2]++;
                }
            }else{
                temp.mat_cli[cli_1][0] = fac_2;
                temp.vet_fac_qnt[fac_2] += temp.mat_cli[cli_1][pos_1 + 2];
                temp.mat_cli[cli_1][2] += temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = 0;
                temp.mat_cli[cli_1][1] = -1;
            }


            //No Cliente 2:
            temp.vet_fac_qnt[fac_2] -= temp.mat_cli[cli_2][pos_2 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<temp.vet_aux_fac_cli[fac_2]; y++){
                if(temp.mat_fac_cli[fac_2][y] == cli_2){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    temp.mat_fac_cli[fac_2][y] = temp.mat_fac_cli[fac_2][temp.vet_aux_fac_cli[fac_2] - 1];
                    temp.vet_aux_fac_cli[fac_2]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_2]; y++) temp.mat_num_incompat[mat_incompat[cli_2][y]][fac_2] -= 1;
            
            //Adiciona a fac_1
            if(temp.mat_cli[cli_2][(pos_2 + 1) % 2] != fac_1){ //Verifica se cli_2 já possui fac_1
                temp.mat_cli[cli_2][pos_2] = fac_1;
                temp.vet_fac_qnt[fac_1] += temp.mat_cli[cli_2][pos_2 + 2];
                temp.mat_fac_cli[fac_1][temp.vet_aux_fac_cli[fac_1]] = cli_2;
                temp.vet_aux_fac_cli[fac_1]++;
    
                for(int y=0; y<vet_aux_incompat[cli_2]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    temp.mat_num_incompat[mat_incompat[cli_2][y]][fac_1]++;
                }
            }else{
                temp.mat_cli[cli_2][0] = fac_1;
                temp.vet_fac_qnt[fac_1] += temp.mat_cli[cli_2][pos_2 + 2];
                temp.mat_cli[cli_2][2] += temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = 0;
                temp.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_1 e nova f_2 sejam iguais (junta ambas) 
            if(temp.mat_cli[cli_1][0] == temp.mat_cli[cli_1][1]){
                temp.mat_cli[cli_1][2] += temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = 0;
                temp.mat_cli[cli_1][1] = -1;
            }
            if(temp.mat_cli[cli_2][0] == temp.mat_cli[cli_2][1]){
                temp.mat_cli[cli_2][2] += temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = 0;
                temp.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_2 seja mais barata que f_1
            if(temp.mat_cli[cli_1][1] != -1
            && mat_custo[cli_1][temp.mat_cli[cli_1][1]] < mat_custo[cli_1][temp.mat_cli[cli_1][0]]){
                aux = temp.mat_cli[cli_1][0];
                temp.mat_cli[cli_1][0] = temp.mat_cli[cli_1][1];
                temp.mat_cli[cli_1][1] = aux;
                aux = temp.mat_cli[cli_1][2];
                temp.mat_cli[cli_1][2] = temp.mat_cli[cli_1][3];
                temp.mat_cli[cli_1][3] = aux;
            }
            if(temp.mat_cli[cli_2][1] != -1
            && mat_custo[cli_2][temp.mat_cli[cli_2][1]] < mat_custo[cli_2][temp.mat_cli[cli_2][0]]){
                aux = temp.mat_cli[cli_2][0];
                temp.mat_cli[cli_2][0] = temp.mat_cli[cli_2][1];
                temp.mat_cli[cli_2][1] = aux;
                aux = temp.mat_cli[cli_2][2];
                temp.mat_cli[cli_2][2] = temp.mat_cli[cli_2][3];
                temp.mat_cli[cli_2][3] = aux;
            }

            calcular_fo(temp);
            if(temp.fo < melhor.fo){
                printf("temp: %i\n",temp.fo);
                memcpy(&melhor, &temp, sizeof(Solucao));
                melhorou = 1;
            }
            memcpy(&temp, &s, sizeof(Solucao));
        }
    }

    if(melhorou){
        #ifdef DBG_BL_MM_SWAP
        printf("FO: %i\n",melhor.fo);
        #endif
        memcpy(&s, &melhor, sizeof(Solucao));
        goto INICIO_BL_MM_SWAP;
    }

    memcpy(&s, &melhor, sizeof(Solucao));
}

void clopenFacilities_autoral_v1(Solucao& s){
    int var_aleatoria = rand() % 100;
    int index_fac, fac_close = -1, fac_open, index_nova_fac, nova_fac;
    int index_cli, cli, pos, aux;
    int p_abrir, p_fechar; //Probabilidade de abrir ou fechar
    int abrir = 0, fechar = 0; // 1- ativo  0- desativado
    int vet_idx[num_fac]; //Vetor auxiliar dos índices das facilidades 
    int vet_idx_nova[num_fac]; //Vetor auxiliar dos índices das novas facilidades 
    int vet_idx_cli[num_cli];
    int escolheu = 0;

    for(int x=0; x<num_fac; x++){ //Inicializa vetores index
        vet_idx[x] = x;
    }
    
    while(abrir == 0 && fechar == 0){
        p_abrir = rand() % 100;
        p_fechar = rand() % 100;

        if(p_abrir <= prob_open * 100) abrir = 1;
        if(p_fechar <= prob_close * 100) fechar = 1;
    }

    if(fechar){ //Fechar facilidade
        for(int index=0; index < num_fac; index++){ //Seleciona uma facilidade aleatória para fechar
            index_fac = (rand() % (num_fac - index)) + index;
            fac_close = vet_idx[index_fac];
            vet_idx[index_fac] = vet_idx[index];
            vet_idx[index] = fac_close;

            if(s.vet_fac_qnt[fac_close] == 0) continue; //Caso já esteja fechada

            for(int idx_cli = s.vet_aux_fac_cli[fac_close] - 1; idx_cli >= 0; idx_cli--){ //Remove cada cliente que é abastecido por fac
                cli = s.mat_fac_cli[fac_close][idx_cli];

                if(s.mat_cli[cli][1] == fac_close){ //Caso fac esteja em f_2
                    pos = 1;
                }else{ // Caso fac esteja em f_1
                    pos = 0;
                }

                //Remove cliente da matriz de lista de clientes por facilidade
                s.vet_aux_fac_cli[fac_close]--;

                //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][fac_close] -= 1;

                s.vet_fac_qnt[fac_close] -= s.mat_cli[cli][pos + 2];

                for(int x=0; x<vet_aux_pref[cli]; x++) vet_idx_nova[x] = mat_fac_pref[cli][x]; //Copia facilidades preferenciais

                //Escolhe nova_fac
                escolheu = 0;
                for(int idx=0; idx < vet_aux_pref[cli]; idx++){
                    index_nova_fac = (rand() % (vet_aux_pref[cli] - idx)) + idx;
                    nova_fac = vet_idx_nova[index_nova_fac];
                    vet_idx_nova[index_nova_fac] = vet_idx_nova[idx];
                    vet_idx_nova[idx] = nova_fac;

                    if(vet_capacidade[nova_fac] - s.vet_fac_qnt[nova_fac] < s.mat_cli[cli][pos + 2]
                    || s.mat_num_incompat[cli][nova_fac] > 0
                    || nova_fac == fac_close
                    ) continue;

                    //Adiciona a nova_fac
                    if(s.mat_cli[cli][(pos + 1) % 2] != nova_fac){ //Verifica se cli já possui nova_fac
                        s.mat_cli[cli][pos] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                        s.vet_aux_fac_cli[nova_fac]++;
            
                        for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                            s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                        }
                    }else{
                        s.mat_cli[cli][0] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_cli[cli][2] += s.mat_cli[cli][3];
                        s.mat_cli[cli][3] = 0;
                        s.mat_cli[cli][1] = -1;
                    }

                    escolheu = 1;
                    break;
                }

                if(!escolheu){ //Caso não tenha encontrado uma facilidade válida, insere aleatoriamente
                    do{
                        index_nova_fac = rand() % vet_aux_pref[cli];
                        nova_fac = vet_idx_nova[index_nova_fac];
                    }while(nova_fac == fac_close);

                    //Adiciona a nova_fac
                    if(s.mat_cli[cli][(pos + 1) % 2] != nova_fac){ //Verifica se cli já possui nova_fac
                        s.mat_cli[cli][pos] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                        s.vet_aux_fac_cli[nova_fac]++;
            
                        for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                            s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                        }
                    }else{
                        s.mat_cli[cli][0] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_cli[cli][2] += s.mat_cli[cli][3];
                        s.mat_cli[cli][3] = 0;
                        s.mat_cli[cli][1] = -1;
                    }
                }

                //Verificação caso f_2 seja mais barata que f_1
                if(s.mat_cli[cli][1] != -1
                && mat_custo[cli][s.mat_cli[cli][1]] < mat_custo[cli][s.mat_cli[cli][0]]){
                    aux = s.mat_cli[cli][0];
                    s.mat_cli[cli][0] = s.mat_cli[cli][1];
                    s.mat_cli[cli][1] = aux;
                    aux = s.mat_cli[cli][2];
                    s.mat_cli[cli][2] = s.mat_cli[cli][3];
                    s.mat_cli[cli][3] = aux;
                }
            }

            break;
        }
    }
    
    if(abrir){ //Abrir facilidade
        for(int index=0; index < num_fac; index++){
            index_fac = (rand() % (num_fac - index)) + index;
            fac_open = vet_idx[index_fac];
            vet_idx[index_fac] = vet_idx[index];
            vet_idx[index] = fac_open;
        
            if(s.vet_fac_qnt[fac_open] > 0
            || fac_open == fac_close
            ) continue; //Caso já esteja aberta

            for(int x=0; x<vet_aux_cli_pref[fac_open]; x++) vet_idx_cli[x] = mat_cli_pref[fac_open][x];

            for(int index=0; index<vet_aux_cli_pref[fac_open] && (vet_capacidade[fac_open] - s.vet_fac_qnt[fac_open]) > 0; index++){ //Verifica se os clientes preferenciais melhoram com a fac_open
                index_cli = (rand() % (vet_aux_cli_pref[fac_open] - index)) + index;
                cli = vet_idx_cli[index_cli];
                vet_idx_cli[index_cli] = vet_idx_cli[index];
                vet_idx_cli[index] = cli;

                pos = MIN(1,s.mat_cli[cli][1] + 1);

                if(mat_custo[cli][s.mat_cli[cli][pos]] < mat_custo[cli][fac_open]
                || (vet_capacidade[fac_open] - s.vet_fac_qnt[fac_open]) < s.mat_cli[cli][pos + 2]
                || s.mat_num_incompat[cli][fac_open] > 0
                ) continue;

                //Remove cliente da matriz de lista de clientes por facilidade
                for(int y=0; y<s.vet_aux_fac_cli[s.mat_cli[cli][pos]]; y++){
                    if(s.mat_fac_cli[s.mat_cli[cli][pos]][y] == cli){
                        //Troca o elemento a ser removido com o último e decrementa o contadora
                        s.mat_fac_cli[s.mat_cli[cli][pos]][y] = s.mat_fac_cli[s.mat_cli[cli][pos]][s.vet_aux_fac_cli[s.mat_cli[cli][pos]] - 1];
                        s.vet_aux_fac_cli[s.mat_cli[cli][pos]]--;
                        break;
                    }
                }

                //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][s.mat_cli[cli][pos]] -= 1;

                s.vet_fac_qnt[s.mat_cli[cli][pos]] -= s.mat_cli[cli][pos + 2];
                
                //Adiciona fac_open
                s.mat_cli[cli][pos] = fac_open;
                s.vet_fac_qnt[fac_open] += s.mat_cli[cli][pos + 2];
                s.mat_fac_cli[fac_open][s.vet_aux_fac_cli[fac_open]] = cli;
                s.vet_aux_fac_cli[fac_open]++;
    
                for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    s.mat_num_incompat[mat_incompat[cli][y]][fac_open]++;
                }

                //Verificação caso f_2 seja mais barata que f_1
                if(s.mat_cli[cli][1] != -1
                && mat_custo[cli][s.mat_cli[cli][1]] < mat_custo[cli][s.mat_cli[cli][0]]){
                    aux = s.mat_cli[cli][0];
                    s.mat_cli[cli][0] = s.mat_cli[cli][1];
                    s.mat_cli[cli][1] = aux;
                    aux = s.mat_cli[cli][2];
                    s.mat_cli[cli][2] = s.mat_cli[cli][3];
                    s.mat_cli[cli][3] = aux;
                }
            }

            break;
        }
    }

    calcular_fo(s);
}

void clopenFacilities(Solucao& s){
    Solucao temp;
    int index_close_fac,index_open_fac,index_nova_fac;
    int close_fac,open_fac,nova_fac;
    int melhorCustoParcial, melhorCustoTotal;
    int melhorFacCustoParcial = -1, melhorFacCustoTotal = -1;
    int cli;
    int pos;

    int vet_idx_fac[num_fac];
    int vet_idx_nova_fac[num_fac];
    int vet_idx_fac_viz[num_fac];

    //Criação do vetor auxiliar de index
    for(int x=0; x<num_fac; x++) vet_idx_fac[x] = x;

    //Escolhe uma facilidade para ser fechada
    for(int idx_close_fac=0; idx_close_fac<num_fac; idx_close_fac++){
        index_close_fac = (rand() % (num_fac - idx_close_fac)) + idx_close_fac;
        close_fac = vet_idx_fac[index_close_fac];
        vet_idx_fac[index_close_fac] = vet_idx_fac[idx_close_fac];
        vet_idx_fac[idx_close_fac] = close_fac;

        //Copia o vetor de facilidades vizinhas
        for(int x=0; x<vet_idx_fac_viz[close_fac]; x++) vet_idx_fac_viz[x] = vet_aux_fac_viz[x];

        //Escolhe uma facilidade vizinha para ser aberta
        for(int idx_open_fac=0; idx_open_fac<vet_aux_fac_viz[close_fac]; idx_open_fac++){
            index_open_fac = (rand() % (num_fac - idx_open_fac)) + idx_open_fac;
            open_fac = vet_idx_fac_viz[index_open_fac];
            vet_idx_fac_viz[index_open_fac] = vet_idx_fac_viz[idx_open_fac];
            vet_idx_fac_viz[idx_open_fac] = open_fac;

            //Verifica se ambas estão abertas ou ambas fechadas
            if((temp.vet_fac_qnt[close_fac] != 0 && temp.vet_fac_qnt[open_fac] != 0)
            || (temp.vet_fac_qnt[close_fac] == 0 && temp.vet_fac_qnt[open_fac] == 0)) continue;

            memcpy(&temp, &s, sizeof(Solucao));

            //Caso close_fac esteja fechada e open_fac aberta: swap
            if(temp.vet_fac_qnt[close_fac] == 0){
                int aux = close_fac;
                close_fac = open_fac;
                open_fac = aux;
            }

            for(int idx_cli = temp.vet_aux_fac_cli[close_fac] - 1; idx_cli >= 0; idx_cli--){ //Remove cada cliente que é abastecido por fac
                cli = temp.mat_fac_cli[close_fac][idx_cli];

                if(temp.mat_cli[cli][1] == close_fac){ //Caso fac esteja em f_2
                    pos = 1;
                }else{ // Caso fac esteja em f_1
                    pos = 0;
                }

                //Remove cliente da matriz de lista de clientes por facilidade
                temp.vet_aux_fac_cli[close_fac]--;

                //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                for(int y=0; y<vet_aux_incompat[cli]; y++) temp.mat_num_incompat[mat_incompat[cli][y]][close_fac] -= 1;

                temp.vet_fac_qnt[close_fac] -= temp.mat_cli[cli][pos + 2];

                //Talvez não precise
                for(int x=0; x<vet_aux_pref[cli]; x++) vet_idx_nova_fac[x] = mat_fac_pref[cli][x]; //Copia facilidades preferenciais

                //Escolhe melhor fac
                for(int idx_fac_pref=0; idx_fac_pref < vet_aux_pref[cli]; idx_fac_pref++){
                    nova_fac = mat_fac_pref[cli][idx_fac_pref];
                    if(vet_capacidade[nova_fac] - temp.vet_fac_qnt[nova_fac] < temp.mat_cli[cli][pos + 2]
                    || temp.mat_num_incompat[cli][nova_fac] > 0
                    || nova_fac == close_fac
                    ) continue;

                    //int custoParcial = mat_custo[cli][nova_fac] * temp.mat_cli[cli][pos + 2];
                    int custoTotal = mat_custo[cli][nova_fac] * temp.mat_cli[cli][pos + 2] + ((temp.vet_fac_qnt[nova_fac] > 0) ? 0 : vet_custo_fixo[nova_fac]);

                    //Salva o valor e a facilidade com melhor Custo Total
                    if(melhorFacCustoTotal == -1 || melhorCustoTotal > custoTotal){
                        melhorFacCustoTotal = nova_fac;
                        melhorCustoTotal = custoTotal;

                        //Salva o valor e a facilidade com melhor Custo Parcial
                        /*if(melhorFacCustoParcial == -1 || melhorCustoParcial > custoParcial){
                            melhorFacCustoParcial = nova_fac;
                            melhorCustoParcial = custoParcial;
                        }*/
                    } 
                }

                //Caso não encontre nenhuma
                if(/*melhorFacCustoParcial == -1 && */melhorFacCustoTotal == -1) continue;

                //Adiciona a melhorFacCustoTotal
                if(temp.mat_cli[cli][(pos + 1) % 2] != melhorFacCustoTotal){ //Verifica se cli já possui melhorFacCustoTotal
                    temp.mat_cli[cli][pos] = melhorFacCustoTotal;
                    temp.vet_fac_qnt[melhorFacCustoTotal] += temp.mat_cli[cli][pos + 2];
                    temp.mat_fac_cli[melhorFacCustoTotal][temp.vet_aux_fac_cli[melhorFacCustoTotal]] = cli;
                    temp.vet_aux_fac_cli[melhorFacCustoTotal]++;
        
                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        temp.mat_num_incompat[mat_incompat[cli][y]][melhorFacCustoTotal]++;
                    }
                }else{
                    temp.mat_cli[cli][0] = melhorFacCustoTotal;
                    temp.vet_fac_qnt[melhorFacCustoTotal] += temp.mat_cli[cli][pos + 2];
                    temp.mat_cli[cli][2] += temp.mat_cli[cli][3];
                    temp.mat_cli[cli][3] = 0;
                    temp.mat_cli[cli][1] = -1;
                }

                //Escolhe nova_fac
                /*int escolheu = 0;
                for(int idx=0; idx < vet_aux_pref[cli]; idx++){
                    index_nova_fac = (rand() % (vet_aux_pref[cli] - idx)) + idx;
                    nova_fac = vet_idx_nova_fac[index_nova_fac];
                    vet_idx_nova_fac[index_nova_fac] = vet_idx_nova_fac[idx];
                    vet_idx_nova_fac[idx] = nova_fac;

                    if(vet_capacidade[nova_fac] - s.vet_fac_qnt[nova_fac] < s.mat_cli[cli][pos + 2]
                    || s.mat_num_incompat[cli][nova_fac] > 0
                    || nova_fac == close_fac
                    ) continue;

                    //Adiciona a nova_fac
                    if(s.mat_cli[cli][(pos + 1) % 2] != nova_fac){ //Verifica se cli já possui nova_fac
                        s.mat_cli[cli][pos] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                        s.vet_aux_fac_cli[nova_fac]++;
            
                        for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                            s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                        }
                    }else{
                        s.mat_cli[cli][0] = nova_fac;
                        s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][pos + 2];
                        s.mat_cli[cli][2] += s.mat_cli[cli][3];
                        s.mat_cli[cli][3] = 0;
                        s.mat_cli[cli][1] = -1;
                    }

                    escolheu = 1;
                    break;
                }*/
            }
        }
    }
}

void swapSuppliers(Solucao& s){
    int parar_1 = 0, parar_2 = 0;
    int idx_cli_1, idx_cli_2;
    int cli_1, cli_2;
    int vet_aux_1[num_cli], vet_aux_2[num_cli]; //Vetores de índices (para aleatoriedade)
    int pos_1, pos_2; //Facilidade escolhida (1 ou 2)
    int fac_1, fac_2; //índice da Facilidade
    int outra_fac;
    int aux; //Variável auxiliar mesmo...
    int cap_resid_1, cap_resid_2; //Capaciade residual
    int delta, carga_1, carga_2;

    for(int x=0; x<num_cli; x++){ //Preenche vetores auxiliares
        vet_aux_1[x] = x;
    }

    for(int index_1=0; index_1 < num_cli && parar_1 == 0; index_1++){ //Escolhe cli_1 (de forma aleatória)
        idx_cli_1 = (rand() % (num_cli - index_1)) + index_1;
        cli_1 = vet_aux_1[idx_cli_1];
        vet_aux_1[idx_cli_1] = vet_aux_1[index_1];
        vet_aux_1[index_1] = cli_1;

        //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
        pos_1 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, s.mat_cli[cli_1][1] + 1); //Confia que funciona ;D
        fac_1 = s.mat_cli[cli_1][pos_1];
        
        #ifdef DBG_SWAPSUPPLIERS
        printf("De Cliente %i:\n",cli_1);
        printf("<%i, %i, %i, %i>\n",s.mat_cli[cli_1][0],s.mat_cli[cli_1][1],s.mat_cli[cli_1][2],s.mat_cli[cli_1][3]);
        printf("Facilidade: %i\n\n",fac_1);
        #endif

        //Calcula capacidade residual da facilidade escolhida
        cap_resid_1 = vet_capacidade[fac_1] - s.vet_fac_qnt[fac_1];
        
        //Copia os clientes preferenciais
        for(int x=0; x<vet_aux_cli_pref[fac_1]; x++){
            vet_aux_2[x] = mat_cli_pref[fac_1][x];
        }
        for(int index_2=0; index_2 < vet_aux_cli_pref[fac_1] && parar_2 == 0; index_2++){ //Escolhe cli_2 (de forma aleatória)
            idx_cli_2 = (rand() % (vet_aux_cli_pref[fac_1] - index_2)) + index_2;
            cli_2 = vet_aux_2[idx_cli_2];
            vet_aux_2[idx_cli_2] = vet_aux_2[index_2];
            vet_aux_2[index_2] = cli_2;

            //Escolhe facilidade com BIAS para f_2, se f_2 for -1 ele sempre resultará em 0
            pos_2 = (1 - (rand() % 100) / (100 - (SS_bias * 100))) * MIN(1, s.mat_cli[cli_2][1] + 1);
            fac_2 = s.mat_cli[cli_2][pos_2];

            #ifdef DBG_SWAPSUPPLIERS
            printf("Para Cliente %i:\n",cli_2);
            printf("<%i, %i, %i, %i>\n",s.mat_cli[cli_2][0],s.mat_cli[cli_2][1],s.mat_cli[cli_2][2],s.mat_cli[cli_2][3]);
            printf("Facilidade: %i\n\n",fac_2);
            #endif

            //Calcula capacidade residual da facilidade escolhida
            cap_resid_2 = vet_capacidade[fac_2] - s.vet_fac_qnt[fac_2];

            if(cli_2 == cli_1
            || fac_1 == fac_2
            || cap_resid_2 < s.mat_cli[cli_1][pos_1 + 2] - s.mat_cli[cli_2][pos_2 + 2]
            || cap_resid_1 < s.mat_cli[cli_2][pos_2 + 2] - s.mat_cli[cli_1][pos_1 + 2]
            ) continue;

            if(mat_bin_incompat[cli_1][cli_2]){ 
                if(s.mat_num_incompat[cli_1][fac_2] > 1
                || s.mat_num_incompat[cli_2][fac_1] > 1
            ) continue;
            }else{
                if(s.mat_num_incompat[cli_1][fac_2] > 0
                || s.mat_num_incompat[cli_2][fac_1] > 0
            ) continue;
            }

            //Calcula delta
            carga_1 = s.mat_cli[cli_1][pos_1 + 2];
            carga_2 = s.mat_cli[cli_2][pos_2 + 2];
            delta = carga_1 * (mat_custo[cli_1][fac_2] - mat_custo[cli_1][fac_1])
            + carga_2 * (mat_custo[cli_2][fac_1] - mat_custo[cli_2][fac_2]);

            //Realiza a troca
            //No Cliente 1:
            s.vet_fac_qnt[fac_1] -= s.mat_cli[cli_1][pos_1 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<s.vet_aux_fac_cli[fac_1]; y++){
                if(s.mat_fac_cli[fac_1][y] == cli_1){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    s.mat_fac_cli[fac_1][y] = s.mat_fac_cli[fac_1][s.vet_aux_fac_cli[fac_1] - 1];
                    s.vet_aux_fac_cli[fac_1]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_1]; y++) s.mat_num_incompat[mat_incompat[cli_1][y]][fac_1] -= 1;
            
            //Adiciona a fac_2
            if(s.mat_cli[cli_1][(pos_1 + 1) % 2] != fac_2){ //Verifica se cli_1 já possui fac_2
                s.mat_cli[cli_1][pos_1] = fac_2;
                s.vet_fac_qnt[fac_2] += s.mat_cli[cli_1][pos_1 + 2];
                s.mat_fac_cli[fac_2][s.vet_aux_fac_cli[fac_2]] = cli_1;
                s.vet_aux_fac_cli[fac_2]++;
    
                for(int y=0; y<vet_aux_incompat[cli_1]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    s.mat_num_incompat[mat_incompat[cli_1][y]][fac_2]++;
                }
            }else{
                s.mat_cli[cli_1][0] = fac_2;
                s.vet_fac_qnt[fac_2] += s.mat_cli[cli_1][pos_1 + 2];
                s.mat_cli[cli_1][2] += s.mat_cli[cli_1][3];
                s.mat_cli[cli_1][3] = 0;
                s.mat_cli[cli_1][1] = -1;
            }


            //No Cliente 2:
            s.vet_fac_qnt[fac_2] -= s.mat_cli[cli_2][pos_2 + 2];

            //Remove cliente da matriz de lista de clientes por facilidade
            for(int y=0; y<s.vet_aux_fac_cli[fac_2]; y++){
                if(s.mat_fac_cli[fac_2][y] == cli_2){
                    //Troca o elemento a ser removido com o último e decrementa o contadora
                    s.mat_fac_cli[fac_2][y] = s.mat_fac_cli[fac_2][s.vet_aux_fac_cli[fac_2] - 1];
                    s.vet_aux_fac_cli[fac_2]--;
                    break;
                }
            }

            //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
            for(int y=0; y<vet_aux_incompat[cli_2]; y++) s.mat_num_incompat[mat_incompat[cli_2][y]][fac_2] -= 1;
            
            //Adiciona a fac_1
            if(s.mat_cli[cli_2][(pos_2 + 1) % 2] != fac_1){ //Verifica se cli_2 já possui fac_1
                s.mat_cli[cli_2][pos_2] = fac_1;
                s.vet_fac_qnt[fac_1] += s.mat_cli[cli_2][pos_2 + 2];
                s.mat_fac_cli[fac_1][s.vet_aux_fac_cli[fac_1]] = cli_2;
                s.vet_aux_fac_cli[fac_1]++;
    
                for(int y=0; y<vet_aux_incompat[cli_2]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    s.mat_num_incompat[mat_incompat[cli_2][y]][fac_1]++;
                }
            }else{
                s.mat_cli[cli_2][0] = fac_1;
                s.vet_fac_qnt[fac_1] += s.mat_cli[cli_2][pos_2 + 2];
                s.mat_cli[cli_2][2] += s.mat_cli[cli_2][3];
                s.mat_cli[cli_2][3] = 0;
                s.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_1 e nova f_2 sejam iguais (junta ambas) 
            if(s.mat_cli[cli_1][0] == s.mat_cli[cli_1][1]){
                s.mat_cli[cli_1][2] += s.mat_cli[cli_1][3];
                s.mat_cli[cli_1][3] = 0;
                s.mat_cli[cli_1][1] = -1;
            }
            if(s.mat_cli[cli_2][0] == s.mat_cli[cli_2][1]){
                s.mat_cli[cli_2][2] += s.mat_cli[cli_2][3];
                s.mat_cli[cli_2][3] = 0;
                s.mat_cli[cli_2][1] = -1;
            }

            //Verificação caso f_2 seja mais barata que f_1
            if(s.mat_cli[cli_1][1] != -1
            && mat_custo[cli_1][s.mat_cli[cli_1][1]] < mat_custo[cli_1][s.mat_cli[cli_1][0]]){
                aux = s.mat_cli[cli_1][0];
                s.mat_cli[cli_1][0] = s.mat_cli[cli_1][1];
                s.mat_cli[cli_1][1] = aux;
                aux = s.mat_cli[cli_1][2];
                s.mat_cli[cli_1][2] = s.mat_cli[cli_1][3];
                s.mat_cli[cli_1][3] = aux;
            }
            if(s.mat_cli[cli_2][1] != -1
            && mat_custo[cli_2][s.mat_cli[cli_2][1]] < mat_custo[cli_2][s.mat_cli[cli_2][0]]){
                aux = s.mat_cli[cli_2][0];
                s.mat_cli[cli_2][0] = s.mat_cli[cli_2][1];
                s.mat_cli[cli_2][1] = aux;
                aux = s.mat_cli[cli_2][2];
                s.mat_cli[cli_2][2] = s.mat_cli[cli_2][3];
                s.mat_cli[cli_2][3] = aux;
            }

            parar_1 = parar_2 = 1; //Para os for's
            //calcular_fo(s);
            s.fo += delta;

        }
    }
}

void changeSupplier(Solucao& s){
    int cli, index_cli, index_fac;
    int antiga_fac,nova_fac; //Facilidade antiga e nova
    int vet_idx_cli[num_cli]; //Vetor de index de clientes (aleatoriedade sem repetição)
    for(int x=0; x<num_cli; x++) vet_idx_cli[x] = x;
    int vet_idx_fac[num_fac]; //Vetor de index de facilidades (aleatoriedade sem repetição)
    int cap_resid_a, cap_resid_b, n_carga;
    int parar = 0;

    for(int x=0; x<num_cli  && parar == 0; x++){
        index_cli =  (rand() % (num_cli - x)) + x;
        cli = vet_idx_cli[index_cli];
        vet_idx_cli[index_cli] = vet_idx_cli[x];
        vet_idx_cli[x] = cli;

        //Copia as facilidades preferenciais para o vetor auxiliar de index
        for(int y=0; y<vet_aux_pref[cli]; y++) vet_idx_fac[y] = mat_fac_pref[cli][y];
        
        if(s.mat_cli[cli][1] != -1){ //Possui 2 facilidades
            antiga_fac = s.mat_cli[cli][1]; /*---LEMBRAR DE TROCAR TODOS---*/
            int base = 0;

            for(int base=0; base<vet_aux_pref[cli]; base++){
                index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                nova_fac = vet_idx_fac[index_fac];
                vet_idx_fac[index_fac] = vet_idx_fac[base];
                vet_idx_fac[base] = nova_fac;
                cap_resid_a = MAX(0,vet_capacidade[s.mat_cli[cli][0]] - s.vet_fac_qnt[s.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                cap_resid_b = MAX(0,vet_capacidade[nova_fac] - s.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade

                if(nova_fac == s.mat_cli[cli][0] 
                || nova_fac == s.mat_cli[cli][1]
                //|| mat_custo[cli][nova_fac] >= mat_custo[cli][s.mat_cli[cli][1]] //Talvez tenha que tirar*
                || s.mat_num_incompat[cli][nova_fac] != 0
                || cap_resid_b <= 0 
                || (cap_resid_a + cap_resid_b) < s.mat_cli[cli][3]  //Verifica se a carga atual do 2º fornecedor
                ) continue;                                         //cabe na capacidade residual de f_1 e nova_fac

                if(mat_custo[cli][s.mat_cli[cli][0]] <= mat_custo[cli][nova_fac]){ //Caso f_1 seja mais barato que nova_fac
                    if(s.mat_cli[cli][3] <= cap_resid_a){ //Caso caiba tudo em f_1
                        s.mat_cli[cli][2] += s.mat_cli[cli][3]; //Passa tudo para f_1
                        s.vet_fac_qnt[s.mat_cli[cli][0]] += s.mat_cli[cli][3];
                        s.vet_fac_qnt[s.mat_cli[cli][1]] -= s.mat_cli[cli][3];

                        //Remove cliente da matriz de lista de clientes por facilidade

                        for(int y=0; y<s.vet_aux_fac_cli[s.mat_cli[cli][1]]; y++){
                            if(s.mat_fac_cli[s.mat_cli[cli][1]][y] == cli){
                                //Troca o elemento a ser removido com o último e decrementa o contadora
                                s.mat_fac_cli[s.mat_cli[cli][1]][y] = s.mat_fac_cli[s.mat_cli[cli][1]][s.vet_aux_fac_cli[s.mat_cli[cli][1]] - 1];
                                s.vet_aux_fac_cli[s.mat_cli[cli][1]]--;
                                break;
                            }
                        }

                        //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                        for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][s.mat_cli[cli][1]] -= 1;
                        
                        s.mat_cli[cli][1] = -1; //remove f_2
                        s.mat_cli[cli][3] = 0;

                    }else{//Caso NÃO caiba tudo em f_1

                        //Quantidade de carga para nova_fac
                        n_carga = s.mat_cli[cli][3] - cap_resid_a;
                        s.mat_cli[cli][2] += cap_resid_a;
                        s.vet_fac_qnt[s.mat_cli[cli][0]] += cap_resid_a;
                        s.vet_fac_qnt[s.mat_cli[cli][1]] -= s.mat_cli[cli][3];

                        //Remove cliente da matriz de lista de clientes por facilidade
                        for(int y=0; y<s.vet_aux_fac_cli[s.mat_cli[cli][1]]; y++){
                            if(s.mat_fac_cli[s.mat_cli[cli][1]][y] == cli){
                                //Troca o elemento a ser removido com o último e decrementa o contadora
                                s.mat_fac_cli[s.mat_cli[cli][1]][y] = s.mat_fac_cli[s.mat_cli[cli][1]][s.vet_aux_fac_cli[s.mat_cli[cli][1]] - 1];
                                s.vet_aux_fac_cli[s.mat_cli[cli][1]]--;
                                break;
                            }
                        }

                        //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                        for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][s.mat_cli[cli][1]] -= 1;
                        
                        s.mat_cli[cli][1] = nova_fac;
                        s.mat_cli[cli][3] = n_carga;
                        s.vet_fac_qnt[nova_fac] +=  n_carga;
                        s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                        s.vet_aux_fac_cli[nova_fac]++;

                        for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                            s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                        }
                    }

                }else{ //Caso nova_fac seja mais barato 

                    //Define a carga de nova_fac
                    n_carga = MIN(s.mat_cli[cli][3],cap_resid_b);
                    s.mat_cli[cli][2] += s.mat_cli[cli][3] - n_carga; //Atribui o restante à f_1
                    s.vet_fac_qnt[s.mat_cli[cli][0]] += s.mat_cli[cli][3] - n_carga;
                    s.vet_fac_qnt[s.mat_cli[cli][1]] -= s.mat_cli[cli][3];

                    //Remove cliente da matriz de lista de clientes por facilidade
                    for(int y=0; y<s.vet_aux_fac_cli[s.mat_cli[cli][1]]; y++){
                        if(s.mat_fac_cli[s.mat_cli[cli][1]][y] == cli){
                            //Troca o elemento a ser removido com o último e decrementa o contadora
                            s.mat_fac_cli[s.mat_cli[cli][1]][y] = s.mat_fac_cli[s.mat_cli[cli][1]][s.vet_aux_fac_cli[s.mat_cli[cli][1]] - 1];
                            s.vet_aux_fac_cli[s.mat_cli[cli][1]]--;
                            break;
                        }
                    }

                    //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                    for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][s.mat_cli[cli][1]] -= 1;
                    
                    s.mat_cli[cli][1] = nova_fac;
                    s.mat_cli[cli][3] = n_carga;
                    s.vet_fac_qnt[nova_fac] +=  n_carga;
                    s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                    s.vet_aux_fac_cli[nova_fac]++;

                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }

                    //Transforma nova_fac em f_1
                    int aux;
                    aux = s.mat_cli[cli][0];
                    s.mat_cli[cli][0] = s.mat_cli[cli][1];
                    s.mat_cli[cli][1] = aux;

                    aux = s.mat_cli[cli][2];
                    s.mat_cli[cli][2] = s.mat_cli[cli][3];
                    s.mat_cli[cli][3] = aux;
                }

                parar = 1;
                break;
            }

        }else{ //Possui 1 facilidade
            antiga_fac = rand() % 2;

            if(antiga_fac == 0){ //Caso seja f_a
                
                for(int base=0; base < vet_aux_pref[cli]; base++){
                    index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                    nova_fac = vet_idx_fac[index_fac];
                    vet_idx_fac[index_fac] = vet_idx_fac[base];
                    vet_idx_fac[base] = nova_fac;

                    cap_resid_a = MAX(0,vet_capacidade[s.mat_cli[cli][0]] - s.vet_fac_qnt[s.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                    cap_resid_b = MAX(0,vet_capacidade[nova_fac] - s.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade


                    if(nova_fac == s.mat_cli[cli][0] //Verifica qualquer tipo de incompatibilidade
                    || s.mat_cli[cli][2] > cap_resid_b
                    || s.mat_num_incompat[cli][nova_fac] != 0
                    ) continue;

                    s.vet_fac_qnt[s.mat_cli[cli][0]] -= s.mat_cli[cli][2];

                    //Remove cliente da matriz de lista de clientes por facilidade
                    for(int y=0; y<s.vet_aux_fac_cli[s.mat_cli[cli][0]]; y++){
                        if(s.mat_fac_cli[s.mat_cli[cli][0]][y] == cli){
                            //Troca o elemento a ser removido com o último e decrementa o contadora
                            s.mat_fac_cli[s.mat_cli[cli][0]][y] = s.mat_fac_cli[s.mat_cli[cli][0]][s.vet_aux_fac_cli[s.mat_cli[cli][0]] - 1];
                            s.vet_aux_fac_cli[s.mat_cli[cli][0]]--;
                            break;
                        }
                    }

                    //Subtrai 1 de cada facilidade incompatível com cli da matriz de incompatibilidade
                    for(int y=0; y<vet_aux_incompat[cli]; y++) s.mat_num_incompat[mat_incompat[cli][y]][s.mat_cli[cli][0]] -= 1;
                    
                    
                    s.mat_cli[cli][0] = nova_fac;
                    s.vet_fac_qnt[nova_fac] += s.mat_cli[cli][2];
                    s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                    s.vet_aux_fac_cli[nova_fac]++;
                    
                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }
                    
                    parar = 1;
                    break;
                }
            }else{ //Caso seja f_b (f_b == -1)

                for(int base=0; base < vet_aux_pref[cli]; base++){
                    index_fac = (rand() % (vet_aux_pref[cli] - base)) + base;
                    nova_fac = vet_idx_fac[index_fac];
                    vet_idx_fac[index_fac] = vet_idx_fac[base];
                    vet_idx_fac[base] = nova_fac;

                    cap_resid_a = MAX(0,vet_capacidade[s.mat_cli[cli][0]] - s.vet_fac_qnt[s.mat_cli[cli][0]]); //Capacidade residual da facilidade antiga
                    cap_resid_b = MAX(0,vet_capacidade[nova_fac] - s.vet_fac_qnt[nova_fac]); //Capacidade residual da nova facilidade
                    
                    if(nova_fac == s.mat_cli[cli][0]
                    || mat_custo[cli][nova_fac] > mat_custo[cli][s.mat_cli[cli][0]]
                    || s.mat_num_incompat[cli][nova_fac] != 0
                    || cap_resid_b == 0
                    ) continue;
                    
                    n_carga = MIN(cap_resid_b, s.mat_cli[cli][2] - 1); //Sempre deixa 1 carga restante pra f_a

                    if(n_carga == 0) continue; //Evita que entre um cliente que fique ocioso (carga 0)

                    //Retira carga de f_a
                    s.mat_cli[cli][2] -= n_carga;
                    s.vet_fac_qnt[s.mat_cli[cli][0]] -= n_carga;

                    //Insere nova_fac
                    s.mat_cli[cli][1] = nova_fac;
                    s.mat_cli[cli][3] = n_carga;
                    s.vet_fac_qnt[nova_fac] += n_carga;
                    s.mat_fac_cli[nova_fac][s.vet_aux_fac_cli[nova_fac]] = cli;
                    s.vet_aux_fac_cli[nova_fac]++;
                    
                    for(int y=0; y<vet_aux_incompat[cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        s.mat_num_incompat[mat_incompat[cli][y]][nova_fac]++;
                    }
                    
                    //Transforma nova_fac em f_1
                    int aux;
                    aux = s.mat_cli[cli][0];
                    s.mat_cli[cli][0] = s.mat_cli[cli][1];
                    s.mat_cli[cli][1] = aux;

                    aux = s.mat_cli[cli][2];
                    s.mat_cli[cli][2] = s.mat_cli[cli][3];
                    s.mat_cli[cli][3] = aux;

                    parar = 1;
                    break;
                }
            }
        }
    }

    calcular_fo(s);
}

void multi_start_greedy_const(Solucao& s){
    int melhor_fac, melhor_cli, index_melhor_cli; //Melhor (i,j)
    int index_cli, index_fac;
    int cli, fac;
    int cli_resid, fac_resid;
    double custo, melhor_custo; //Double pois pode receber a parte fracionada da amortização
    int n_carga;
    int equal_bests = 1;

    memset(&s.vet_aux_fac_cli,0,sizeof(s.vet_aux_fac_cli));
    memset(&s.vet_fac_qnt,0,sizeof(s.vet_fac_qnt));
    memset(&s.mat_num_incompat,0,sizeof(s.mat_num_incompat));
    for(int x=0; x<num_cli; x++){
        s.mat_cli[x][0] = -1;
        s.mat_cli[x][1] = -1;
        s.mat_cli[x][2] = 0;
        s.mat_cli[x][3] = 0;
    }

    int vet_aux_cli[num_cli], num_aux_cli;
    int vet_aux_fac[num_fac];

    for(int x=0; x<num_cli; x++) vet_aux_cli[x] = x; //Constrói vetor auxiliar de índex de cliente

    INICIO_HEU_CONST_GUL : ; //Ponto de início da heurística construtiva gulosa
    num_aux_cli = num_cli; //Inicia com todos os clientes (elimina os que já estiverem abastecidos)
    while(num_aux_cli > 0){//Enquanto houver cliente não abastecido
        melhor_custo = INT_MAX; //"infinito"
        melhor_cli = -1;

        for(int idx_cli=0; idx_cli < num_aux_cli; idx_cli++){
            index_cli = (rand() % (num_aux_cli - idx_cli)) + idx_cli;
            cli = vet_aux_cli[index_cli];
            vet_aux_cli[index_cli] = vet_aux_cli[idx_cli];
            vet_aux_cli[idx_cli] = cli;

            cli_resid = MAX(0,vet_demanda[cli] - (s.mat_cli[cli][2] + s.mat_cli[cli][3])); //Capacidade residual de cliente

            for(int x=0; x<vet_aux_pref[cli]; x++) vet_aux_fac[x] = mat_fac_pref[cli][x]; //Copia faciliades preferenciais
            
            for(int idx_fac=0; idx_fac < vet_aux_pref[cli]; idx_fac++){
                index_fac = (rand() % (vet_aux_pref[cli] - idx_fac)) + idx_fac;
                fac = vet_aux_fac[index_fac];
                vet_aux_fac[index_fac] = vet_aux_fac[idx_fac];
                vet_aux_fac[idx_fac] = fac;

                fac_resid = MAX(0,vet_capacidade[fac] - s.vet_fac_qnt[fac]); //Capacidade residual da facilidade

                if(fac_resid == 0
                || s.mat_num_incompat[cli][fac] != 0
                || (s.mat_cli[cli][0] != -1 && fac_resid < cli_resid)
                ) continue;

                //Calculo do custo (com a amortização)
                custo = mat_custo[cli][fac] * 1.0
                + 
                ((s.vet_fac_qnt[fac] > 0) ? 0 : greedy_amortize_factor * (1.0 * vet_custo_fixo[fac] * vet_demanda[cli]) / vet_capacidade[fac] * 1.0);

                if(custo < melhor_custo){ //Se for melhor, sempre atribui a melhor
                    index_melhor_cli = idx_cli;
                    melhor_cli = cli;
                    melhor_fac = fac;
                    melhor_custo = custo;
                    equal_bests = 1;
                } else if(custo < 1.0 * melhor_custo + equal_tolerance){ //Se estiver dentro do raio, aceita com uma probabilidade
                    equal_bests++;
                    if(rand() % equal_bests == 0){
                        index_melhor_cli = idx_cli;
                        melhor_cli = cli;
                        melhor_fac = fac;
                    }
                }
            }
        }

        if(melhor_cli == -1){ //Caso não tenha encontrado nenhum par reinicia tudo do 0
            printf("Não foi encontrado nenhum par, reiniciando construtiva gulosa\n");
            goto INICIO_HEU_CONST_GUL;
        }

        cli_resid = vet_demanda[melhor_cli] - (s.mat_cli[melhor_cli][2] + s.mat_cli[melhor_cli][3]);
        fac_resid = vet_capacidade[melhor_fac] - s.vet_fac_qnt[melhor_fac];
        int pos = MIN(1, s.mat_cli[melhor_cli][0] + 1);
        n_carga = MIN(cli_resid, fac_resid);

        //Insere nova_fac
       /*if(s.mat_cli[melhor_cli][(pos + 1) % 2] != melhor_fac){ //Verifica se melhor_cli já possui melhor_fac
            s.mat_cli[melhor_cli][pos] = melhor_fac;
            s.mat_cli[melhor_cli][pos + 2] = n_carga;
            s.vet_fac_qnt[melhor_fac] += n_carga;
            s.mat_fac_cli[melhor_fac][s.vet_aux_fac_cli[melhor_fac]] = melhor_cli;
            s.vet_aux_fac_cli[melhor_fac]++;

            for(int y=0; y<vet_aux_incompat[melhor_cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                s.mat_num_incompat[mat_incompat[melhor_cli][y]][melhor_fac]++;
            }
        }else{
            s.vet_fac_qnt[melhor_fac] += n_carga;
            s.mat_cli[melhor_cli][2] += n_carga;
            s.mat_cli[melhor_cli][3] = 0;
            s.mat_cli[melhor_cli][1] = -1;
        }
        Acho que não tem como aparecer a mesma facilidade para o mesmo cliente.
        Acredito que esse seja o jeito certo:*/
        s.mat_cli[melhor_cli][pos] = melhor_fac;
        s.mat_cli[melhor_cli][pos + 2] = n_carga;
        s.vet_fac_qnt[melhor_fac] += n_carga;
        s.mat_fac_cli[melhor_fac][s.vet_aux_fac_cli[melhor_fac]] = melhor_cli;
        s.vet_aux_fac_cli[melhor_fac]++;

        for(int y=0; y<vet_aux_incompat[melhor_cli]; y++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
            s.mat_num_incompat[mat_incompat[melhor_cli][y]][melhor_fac]++;
        }

        //Transforma nova_fac em f_1 caso tenha menor custo
        if(s.mat_cli[melhor_cli][1] != -1
            && mat_custo[melhor_cli][s.mat_cli[melhor_cli][1]] < mat_custo[melhor_cli][s.mat_cli[melhor_cli][0]]){
            int aux;
            aux = s.mat_cli[melhor_cli][0];
            s.mat_cli[melhor_cli][0] = s.mat_cli[melhor_cli][1];
            s.mat_cli[melhor_cli][1] = aux;
    
            aux = s.mat_cli[melhor_cli][2];
            s.mat_cli[melhor_cli][2] = s.mat_cli[melhor_cli][3];
            s.mat_cli[melhor_cli][3] = aux;
        }

        //Verifica se melhor_cli foi completamente atendido e o remove
        if(vet_demanda[melhor_cli] - (s.mat_cli[melhor_cli][2] + s.mat_cli[melhor_cli][3]) <= 0){
            vet_aux_cli[index_melhor_cli] = vet_aux_cli[num_aux_cli - 1];
            num_aux_cli--; //Diminui o tamanho do vetor 
        }
    }

    calcular_fo(s);
}

void heu_const_ale(Solucao& s){
    int n_fac; //1 ou 2 facilidades
    int n_dem; //Quantidade de demanda da 1ª facilidade
    int facilidade;
    int id,aux;
    //int incompat;
    int vet_aux[MAX_FAC];
    memset(&s.vet_aux_fac_cli,0,sizeof(s.vet_aux_fac_cli));
    memset(&s.vet_fac_qnt,0,sizeof(s.vet_fac_qnt));
    memset(&s.mat_num_incompat,0,sizeof(s.mat_num_incompat));

    for(int cli=0; cli<num_cli; cli++){
        n_fac = (rand() % 4); //se == 3 então 2 facilidades 
        s.mat_cli[cli][0] = -1;
        s.mat_cli[cli][1] = -1;
        s.mat_cli[cli][2] = 0;
        s.mat_cli[cli][3] = 0;

        //gera um vetor com os índices dos fac preferenciais para serem escolhidos aleatoriamente
        for(int x=0; x<vet_aux_pref[cli]; x++){
            vet_aux[x] = mat_fac_pref[cli][x];
        }

        //Determina o valor da demanda da 1ª facilidade
        if(n_fac < 3){
            n_dem = vet_demanda[cli];
        }else{
            n_dem = (rand() % (vet_demanda[cli] - 1)) + 1;
        }

        //Escolhe a 1ª facilidade
        for(int index=0; index<vet_aux_pref[cli]; index++){
            id = (rand() % (vet_aux_pref[cli] - index)) + index; //Escolhe aleatoriamente
            facilidade = vet_aux[id];
            aux = vet_aux[index];
            vet_aux[index] = vet_aux[id];
            vet_aux[id] = aux;

            #ifdef CONST_ALE_MODO_1
            if(s.vet_fac_qnt[facilidade] + n_dem > vet_capacidade[facilidade]){
                continue;
            }
            #endif
            #ifdef CONST_ALE_MODO_2
            if(n_fac < 3 && s.vet_fac_qnt[facilidade] + vet_demanda[cli] > vet_capacidade[facilidade]){
                continue;
            }
            #endif

           
            if(s.mat_num_incompat[cli][facilidade]==0){ //Caso não tenha incompatibilidade
                //Atribui facilidade ao cliente
                s.mat_cli[cli][0] = facilidade;
                s.mat_cli[cli][2] = n_dem;
                s.vet_fac_qnt[facilidade] +=  n_dem;
                s.mat_fac_cli[facilidade][s.vet_aux_fac_cli[facilidade]] = cli;
                s.vet_aux_fac_cli[facilidade]++;

                for(int x=0; x<vet_aux_incompat[cli]; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    s.mat_num_incompat[mat_incompat[cli][x]][facilidade]++;
                }

                /*Forma pior:
                for(int x=0; x<num_cli; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    if(mat_bin_incompat[cli][x]){
                        s.mat_num_incompat[x][facilidade]++;
                    }
                }*/
                break;
            }
        }

        //Verifica se não houve facilidade compatível (Tratamento de Exceção)
        if(s.mat_cli[cli][0] == -1){
            id = rand() % vet_aux_pref[cli]; //Escolhe aleatoriamente
            facilidade = vet_aux[id];
            //facilidade = rand() % vet_aux_pref[cli];
            s.mat_cli[cli][0] = facilidade; //Insere em uma facilidade aleatória (da lista de pref)
            s.mat_cli[cli][2] = n_dem;
            s.vet_fac_qnt[facilidade] +=  n_dem;
            s.mat_fac_cli[facilidade][s.vet_aux_fac_cli[facilidade]] = cli;
            s.vet_aux_fac_cli[facilidade]++;

            for(int x=0; x<vet_aux_incompat[cli]; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                s.mat_num_incompat[mat_incompat[cli][x]][facilidade]++;
            }

            /*Forma pior:
            for(int x=0; x<num_cli; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                if(mat_bin_incompat[cli][x]){
                    s.mat_num_incompat[x][facilidade]++;
                }
            }*/
        }


        //Adiciona 2ª facilidade se permitido
        if(n_fac == 3){
            for(int index=0; index<vet_aux_pref[cli]; index++){
                id = (rand() % (vet_aux_pref[cli] - index)) + index; //Escolhe aleatoriamente
                facilidade = vet_aux[id]; 
                aux = vet_aux[0];
                vet_aux[0] = vet_aux[id];
                vet_aux[id] = aux;
                
                //Verifica se extrapola a capacidade com o restante da demanda de cli
                if(s.vet_fac_qnt[facilidade] + vet_demanda[cli] - n_dem > vet_capacidade[facilidade] || facilidade == s.mat_cli[cli][0]){
                    continue;
                }

                
                if(s.mat_num_incompat[cli][facilidade] == 0){ //Caso não tenha incompatibilidade
                    //Atribui facilidade ao cliente
                    s.mat_cli[cli][1] = facilidade;
                    s.mat_cli[cli][3] = vet_demanda[cli] - n_dem;
                    s.vet_fac_qnt[facilidade] +=  vet_demanda[cli] - n_dem;
                    s.mat_fac_cli[facilidade][s.vet_aux_fac_cli[facilidade]] = cli;
                    s.vet_aux_fac_cli[facilidade]++;

                    for(int x=0; x<vet_aux_incompat[cli]; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        s.mat_num_incompat[mat_incompat[cli][x]][facilidade]++;
                    }

                    /*Forma pior:
                    for(int x=0; x<num_cli; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                        if(mat_bin_incompat[cli][x]){
                            s.mat_num_incompat[x][facilidade]++;
                        }
                    }*/
                    break;
                }
            }

            //Verifica se não houve facilidade compatível (Tratamento de Exceção)
            if(s.mat_cli[cli][1] == -1){
                do{
                    id = rand() % vet_aux_pref[cli]; //Escolhe aleatoriamente
                    facilidade = vet_aux[id];
                    //facilidade = rand() % vet_aux_pref[cli]; 
                }while(facilidade == s.mat_cli[cli][0]);

                s.mat_cli[cli][1] = facilidade; //Insere em uma facilidade aleatória (da lista de pref)
                s.mat_cli[cli][3] = vet_demanda[cli] - n_dem;
                s.vet_fac_qnt[facilidade] +=  vet_demanda[cli] - n_dem;
                s.mat_fac_cli[facilidade][s.vet_aux_fac_cli[facilidade]] = cli;
                s.vet_aux_fac_cli[facilidade]++;

                for(int x=0; x<vet_aux_incompat[cli]; x++){ //Atualiza a matriz do número de incompatibilidade de cada cliente para cada facilidade
                    s.mat_num_incompat[mat_incompat[cli][x]][facilidade]++;
                }
            }

            //Verifica se f_2 é mais barato que f_1 e troca
            if(mat_custo[cli][s.mat_cli[cli][1]] < mat_custo[cli][s.mat_cli[cli][0]]){
                aux = s.mat_cli[cli][0];
                s.mat_cli[cli][0] = s.mat_cli[cli][1];
                s.mat_cli[cli][1] = aux;
                aux = s.mat_cli[cli][2];
                s.mat_cli[cli][2] = s.mat_cli[cli][3];
                s.mat_cli[cli][3] = aux;
            }
        }

    }

    calcular_fo(s);
}

void cs_gsa(Solucao &s, const int& tempoLimite){
    Solucao melhor,temp;
    memcpy(&melhor, &s, sizeof(s));

    double TempAtual = Tzero;
    int eventos=0;
    int last_impr = 0;
    double p;
    int delta;
    double prob,x;
    time_t h = clock();

    while(TempAtual > Tc && (double)(clock() - h) / CLOCKS_PER_SEC < tempoLimite){
        for(int iter=0; iter < SAmax; iter++){
            last_impr++;
            p = (double)(iter+1) / SAmax;
            memcpy(&temp, &s, sizeof(Solucao));

            gerarVizinho_v2(temp,TempAtual,melhor.fo,p);

            delta = temp.fo - s.fo;
            if(delta < 0){
                memcpy(&s, &temp, sizeof(Solucao));

                if(temp.fo < melhor.fo){
                    memcpy(&melhor, &temp, sizeof(Solucao));
                    last_impr = 0;
                    //printf("MELHOROU: %i\n",melhor.fo);
                }
            }else{
                x = (rand() % 101) / 100;
                prob = pow(1 + (q_a - 1) * ((temp.fo - s.fo) / TempAtual),1/(q_a - 1));

                if(x < prob){
                    memcpy(&s, &temp, sizeof(s));
                }
            }

            //cluster_search(temp,melhor);
        }

        printf("TEMPERATURA: %f\nFO: %i\nMelhor FO:%i\n",TempAtual,temp.fo,melhor.fo);

        TempAtual *= a;
        /*if(last_impr >= r_i && eventos <= r_max){
            TempAtual = r_f;
            eventos++;
            last_impr = 0;
        }*/
       if(TempAtual <= TempChange && eventos < r_max){
        TempAtual = r_f;
        eventos++;
        last_impr = 0;
       }
    }

    memcpy(&s, &melhor, sizeof(s));
}

void simulated_annealing(Solucao &s, const int& tempoLimite){
    Solucao melhor;
    Solucao temp;
    Solucao nv;
    memcpy(&melhor, &s, sizeof(Solucao));
    memcpy(&temp, &s, sizeof(Solucao));
    float T = Tzero;
    int delta=0;
    time_t h;

    h = clock();
    while(T > Tc && (double)(clock() - h) / CLOCKS_PER_SEC < tempoLimite){
        printf("-----------------------------------------------\n");
        for(int x=0; x<SAmax; x++){
            memcpy(&temp, &s, sizeof(Solucao));
            gerarVizinho(temp);
            printf("nv FO: %i\n",temp.fo);
            
            delta = temp.fo - s.fo;
            if(delta < 0){
                memcpy(&s, &temp, sizeof(Solucao));

                if(temp.fo < melhor.fo){
                    //bl_pm_change(temp);
                    memcpy(&melhor, &temp, sizeof(Solucao));

                    #ifdef DBG_SA
                    //escrever_solucao(melhor);
                    printf("Melhorou FO: %i\n",melhor.fo);
                    #endif
                }
            }else{
                float prob = (rand() % 101) / 100.0f;
                if(prob < exp((-delta)/T)){
                    memcpy(&s, &temp, sizeof(Solucao));
                }
            }
        }

        /*bl_pm_change(temp);
        if(temp.fo < melhor.fo){
            memcpy(&melhor, &temp, sizeof(Solucao));
        }*/
        T = a * T;
    }
    memcpy(&s, &melhor, sizeof(Solucao));
}

void GSA(Solucao &s, const int& tempoLimite){
    clock_t h = clock();
    Solucao melhor;
    Solucao temp;
    double T = T0;
    double delta = 0;
    double p;
    int it = 1;

    memcpy(&melhor, &s, sizeof(Solucao));

    while((double)(clock() - h) / CLOCKS_PER_SEC < tempoLimite){
        T = atualizarTemperaturaGSA(T0,it,qT);
        gerarVizinhoGSA(s,temp);
        delta = temp.fo - s.fo;
        p = p_aceitacaoGSA(T,delta,qA);

        if((delta <= 0.0) || (rnd() < p)){
            memcpy(&s, &temp, sizeof(Solucao));

            if(s.fo < melhor.fo){
                bl_pm_change(s);
                memcpy(&melhor, &s, sizeof(Solucao));
            }
        }

        #ifdef DBG_GSA
        printf("TEMPERATURA: %f\nFO: %i\nMelhor FO:%i\n", T, s.fo, melhor.fo);
        #endif
        it++;
    }

    memcpy(&s, &melhor, sizeof(Solucao));
}

void gerarVizinho(Solucao& s){
    int probabilidade = rand() % 100;

    if(probabilidade <= p_CF * 100){ //Probabilidade ClopenFacilities
        clopenFacilities_autoral_v1(s);
    }else if(probabilidade <= p_SS * 100){ //Probabilidade SwapSuppliers
        swapSuppliers(s);
    }else{ //Proababilidade de ChangeSuppliers
        changeSupplier(s);
    }
}

void gerarVizinho_v2(Solucao &s,const double& TempAtual,const int& melhorFo,const double& p){
    double i = MAX(0,(s.fo - melhorFo/*Possível Inversão*/) / (MAX(fabs(melhorFo), 1e-10))); //Garantir: melhorFo > 0
    double lambda = (lambda_i * (TempAtual / Tzero) + lambda_f * (1 - (TempAtual / Tzero))) * MAX(1e-10/*0.0001 ou 1e-10*/,(1 - p) * (i * p)/*Pode ser muito pequeno*/);

    double u = 0.0;
    while(u==0.0||u==1.0){
        u = (rand() % 101) / 100.0;
    }
    double x = lambda * tan(M_PI * (u - 0.5));
    
    int num_passos = MAX(1,(int) fabs(x)); //Define o Valor-Piso de passos
    num_passos = MIN(lim_super, num_passos); //Define o Valor-Teto de passos
    
    //printf("NumPassos: %i\n",num_passos);

    for(int passo=0; passo<num_passos; passo++){
        int probabilidade = rand() % 100;

        if(probabilidade <= p_CF * 100){ //Probabilidade ClopenFacilities
            clopenFacilities_autoral_v1(s);
        }else if(probabilidade <= p_SS * 100){ //Probabilidade SwapSuppliers
            swapSuppliers(s);
        }else{ //Proababilidade de ChangeSuppliers
            changeSupplier(s);
        }
    }
}

void gerarVizinhoGSA(Solucao &s, Solucao& temp){

    memcpy(&temp, &s, sizeof(Solucao));

    int num_passos = calcular_num_passos();

    for(int passo=0; passo<num_passos; passo++){
        int probabilidade = rand() % 100;

        if(probabilidade <= p_CF * 100){ //Probabilidade ClopenFacilities
            clopenFacilities_autoral_v1(temp);
        }else if(probabilidade <= p_SS * 100){ //Probabilidade SwapSuppliers
            swapSuppliers(temp);
        }else{ //Proababilidade de ChangeSuppliers
            changeSupplier(temp);
        }
    }

    calcular_fo(temp);
}

double atualizarTemperaturaGSA(double T0, int t, double qT){
    if(fabs(qT - 1.0) < 1e-12){
        return T0 / log(2.0 + t);
    }

    double num = pow(2.0, qT - 1.0) - 1.0;
    double den = pow(1.0 + t, qT - 1.0) - 1.0;
    if(fabs(den) < 1e-12){
        den = 1e-12;
    }
    return T0 * num / den;
}

double p_aceitacaoGSA(double T, double delta, double qA){
    if(delta <= 0.0){
        return 1.0;
    }
    if(fabs(qA - 1.0) < 1e-12){
        return exp(-delta / T);
    }
    double base = 1.0 - (1.0 - qA) * delta / T;
    if(base <= 0.0){
        return 0.0;
    }
    return pow(base, 1.0 / (1.0 - qA));
}

int calcular_num_passos(){
    double r = rnd();

    //int num = 1 + (int)pow(r, -1.0 / (qV - 1.0));
    int num = (int)pow(r, -1.0 / (qV - 1.0));
    if(num > num_fac){
        return num_fac;
    }

    return num;
}

double rnd(){ //Retorna x = [0,1]
    return (double)rand() / RAND_MAX;
}

void calcular_fo(Solucao& s){
    s.fo = 0;

    for(int x=0; x<num_cli; x++){
        /*if(s.mat_cli[x][0] != -1){ //Talvez não precise dessa verificação
            s.fo += mat_custo[x][s.mat_cli[x][0]] * s.mat_cli[x][2];
        }*/
        s.fo += mat_custo[x][s.mat_cli[x][0]] * s.mat_cli[x][2];

        s.fo += s.mat_num_incompat[x][s.mat_cli[x][0]] * ALPHA; // Caso seja tratado na gerarVizinho e construtiva, pode comentar

        if(s.mat_cli[x][1] != -1){
            s.fo += mat_custo[x][s.mat_cli[x][1]] * s.mat_cli[x][3];
        }
    }

    for(int x=0; x<num_fac; x++){
        s.fo += MIN(1,s.vet_fac_qnt[x]) * vet_custo_fixo[x];
        s.fo += MAX(0,s.vet_fac_qnt[x] - vet_capacidade[x]) * BETA; // Caso seja tratado na gerarVizinho e construtiva, pode comentar
    }
}

void escrever_solucao(Solucao& s){
    for(int x=0; x<num_cli; x++){
        printf("<%i, %i, %i, %i>\n",s.mat_cli[x][0],s.mat_cli[x][1], s.mat_cli[x][2], s.mat_cli[x][3]);
    }

    printf("\nFO: %i\n", s.fo);
}

void testar_dados(char* arq){
    FILE* arquivo;
    if(arq == ""){
        arquivo = stdout;
    }else{
        arquivo = fopen(arq,"w");
    }

    fprintf(arquivo, "Facilities = %i;\n", num_fac);
    fprintf(arquivo,"Customers = %i;\n\n",num_cli);

    fprintf(arquivo, "Capacity = [");
    for (int x = 0; x < num_fac; x++) {
        fprintf(arquivo, "%d", vet_capacidade[x]);
        if (x < num_fac - 1) fprintf(arquivo, ", ");
    }
    fprintf(arquivo, "];\n");

    fprintf(arquivo, "FixedCost = [");
    for (int x = 0; x < num_fac; x++) {
        fprintf(arquivo, "%d", vet_custo_fixo[x]);
        if (x < num_fac - 1) fprintf(arquivo, ", ");
    }
    fprintf(arquivo, "];\n");

    fprintf(arquivo, "Demand = [");
    for (int x = 0; x < num_cli; x++) {
        fprintf(arquivo, "%d", vet_demanda[x]);
        if (x < num_cli - 1) fprintf(arquivo, ", ");
    }
    fprintf(arquivo, "];\n");

    fprintf(arquivo, "ShippingCost = [|");
    for (int x = 0; x < num_cli; x++) {
        for (int y = 0; y < num_fac; y++) {
            fprintf(arquivo, "%d", mat_custo[x][y]);
            if (y < num_fac - 1) fprintf(arquivo, ", ");
        }
        if (x < num_cli - 1) {
            fprintf(arquivo, "\n              |");
        } else {
            fprintf(arquivo, "|];\n");
        }
    }

    fprintf(arquivo, "\n");
    fprintf(arquivo, "Incompatibilities = %d;\n", num_incompat);
    fprintf(arquivo, "IncompatiblePairs = [|");
    int primeiro = 1;
    for (int i = 0; i < num_cli; i++) {
        for (int j = i + 1; j < num_cli; j++) {
            if (mat_bin_incompat[i][j]) {
                if (!primeiro) fprintf(arquivo, "|");
                fprintf(arquivo, "%d, %d", i + 1, j + 1);
                primeiro = 0;
            }
        }
    }
    fprintf(arquivo, "|];\n");

    if (arq != "") fclose(arquivo);

    fclose(arquivo);

}

void ler_dados(){
    FILE* arquivo = fopen(arq,"r");
    char lixo;

    fscanf(arquivo,"Facilities = %i;\n",&num_fac);
    fscanf(arquivo,"Customers = %i;\n\n",&num_cli);

    fscanf(arquivo,"Capacity = [");
    for(int x=0; x<num_fac; x++){
        fscanf(arquivo,"%i, ",&vet_capacidade[x]);
    }
    fscanf(arquivo,"];\n");
    
    fscanf(arquivo,"FixedCost = [");
    for(int x=0; x<num_fac; x++){
        fscanf(arquivo, "%i, ",&vet_custo_fixo[x]);
    }
    fscanf(arquivo,"];\n");

    fscanf(arquivo,"Demand = [");
    for(int x=0; x<num_cli; x++){
        fscanf(arquivo,"%i, ", &vet_demanda[x]);
    }
    fscanf(arquivo,"];\n");

    fscanf(arquivo,"ShippingCost = [|");
    for(int x=0; x<num_cli;x++){
        for(int y=0; y<num_fac; y++){
            fscanf(arquivo,"%i%*[, \t\n] ", &mat_custo[x][y]);
        }
        fscanf(arquivo,"|");
    }
    fscanf(arquivo,"];\n");

    fscanf(arquivo,"Incompatibilities = ");
    fscanf(arquivo,"%i;\n",&num_incompat);

    fscanf(arquivo,"IncompatiblePairs = [|");
    int cli_1,cli_2;
    memset(&vet_aux_incompat,0,sizeof(vet_aux_incompat));
    memset(&mat_bin_incompat,0,sizeof(mat_bin_incompat));
    for(int x=0; x<num_incompat; x++){
        fscanf(arquivo,"%i, %i |",&cli_1,&cli_2);
        cli_1--; //Possível erro
        cli_2--; //Possível erro
        mat_incompat[cli_1][vet_aux_incompat[cli_1]] = cli_2;
        vet_aux_incompat[cli_1]++;

        //A relação de incompatibilidade é simétrica (Equação 11)
        mat_incompat[cli_2][vet_aux_incompat[cli_2]] = cli_1;
        vet_aux_incompat[cli_2]++;

        mat_bin_incompat[cli_1][cli_2] = 1;
        mat_bin_incompat[cli_2][cli_1] = 1;
    }
    fclose(arquivo);

    num_pref = beta * ((int)sqrt(num_fac));
}

void criar_mat_pref(){
    int vet_temp[num_fac]; //Vetor cópia dos custos
    int vet_index[num_fac]; //Vetor de índices (referente ao vet_temp)

    //Inicia vetores auxiliares com valor 0
    memset(&vet_aux_pref,0,sizeof(vet_aux_pref));
    memset(&vet_aux_cli_pref,0,sizeof(vet_aux_cli_pref));


    for(int x=0; x<num_cli; x++){ //Copia valores de distância de cada facilidade
        for(int y=0; y<num_fac; y++){
            vet_temp[y] = mat_custo[x][y];
            vet_index[y] = y;
        }
        quicksortCentro(vet_temp,vet_index,0,num_fac-1);

        for(int y=0; y<num_fac; y++){
            if(y <= num_pref || vet_temp[y] <= vet_temp[0] + lambda){
                //Adiciona facilidade na lista de preferência (p/ cada cliente)
                mat_fac_pref[x][vet_aux_pref[x]] = vet_index[y];
                vet_aux_pref[x]++;
                
                //Adiciona cliente na lista de preferência (p/ cada facilidade)
                mat_cli_pref[vet_index[y]][vet_aux_cli_pref[vet_index[y]]] = x;
                vet_aux_cli_pref[vet_index[y]]++;
            }else{
                break;
            }
        }
        
    }
}

void criar_mat_viz(){
    int existe;
    int fac1,fac2;

    for(int cli=0; cli<num_cli; cli++){
        for(int idx_fac1=0; idx_fac1<vet_aux_pref[cli]; idx_fac1++){
            fac1 = mat_fac_pref[cli][idx_fac1];
            for(int idx_fac2=idx_fac1 + 1; idx_fac2<vet_aux_pref[cli]; idx_fac2++){
                fac2 = mat_fac_pref[cli][idx_fac2];
                existe = 0;
                for(int x=0; x<vet_aux_fac_viz[fac1]; x++){ //Verifica todas as fac vizinhas de fac1
                    if(mat_fac_viz[fac1][x] == fac2){
                        existe = 1;
                        break;
                    }
                }

                if(!existe){
                    mat_fac_viz[fac1][vet_aux_fac_viz[fac1]] = fac2;
                    vet_aux_fac_viz[fac1]++;

                    mat_fac_viz[fac2][vet_aux_fac_viz[fac2]] = fac1;
                    vet_aux_fac_viz[fac2]++;
                }
            }
        }
    }
}

void quicksortCentro(int* vetor,int* vetor2,int inicio, int fim){
    if(fim-inicio>=1){
        int l=inicio, r=fim, pivo=(inicio+fim)/2, aux=0;

        /*Posiciona os valores menores que o pivo � sua esquerda e,
        � sua direita, os valores maiores*/
        while(l <= r){

            //Encontra um valor menor que o piv� no vetor � direita
            while(vetor[l] < vetor[pivo]){
                l++;
            }

            //Encontra um valor maior que o piv� no vetor � esquerda
            while(vetor[r] > vetor[pivo]){
                r--;
            }

            //Permuta os valores encontrados
            if(l<=r){
                aux = vetor[r];
                vetor[r] = vetor[l];
                vetor[l] = aux;
                //Implementação para o problema atual
                aux = vetor2[r];
                vetor2[r] = vetor2[l];
                vetor2[l] = aux;
                //-----------------
                l++;
                r--;
            }
        }

        //Executa recursivamente o vetor do lado esquerdo ao pivo
        quicksortCentro(vetor,vetor2,inicio,l-1);

        //Executa recursivamente o vetor do lado direito ao pivo
        quicksortCentro(vetor,vetor2,l,fim);
    }
}
