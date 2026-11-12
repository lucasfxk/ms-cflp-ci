#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <climits>
#include "estruturas.h"

#define MAX(X,Y) ((X > Y) ? X : Y)
#define MIN(X,Y) ((X < Y) ? X : Y)


void construtivaAleatoriaGulosa(Solucao& s){
    int contMelhor;
    int melhorCli, melhorFac;
    int quantidade;
    double custo, melhorCusto;
    bool compativel = false, achou = false;

    int cargaTotal=0, cargaAtual=0;
    for(int x=0;x<num_cli;x++)cargaTotal += vet_demanda[x];
    

    int indexCli[num_cli], indexFac[num_fac];
    criarVetorDeIndex(indexCli,num_cli);
    criarVetorDeIndex(indexFac, num_fac);

    memset(s.vetorConsumoFac, 0, sizeof(s.vetorConsumoFac));
    montarMatrizSolucao(s);

    while(cargaAtual < cargaTotal){


        for(int i=0; i < num_cli; i++){
            int idxCli = (rand() % (num_cli - i)) + i;
            int cli = indexCli[idxCli];
            indexCli[idxCli] = indexCli[i];
            indexCli[i] = cli;

            if(vet_demanda[cli] - (s.matriz[cli][2] + s.matriz[cli][3]) <= 0) continue; //

            contMelhor = 0;
            achou = false;
            melhorCusto = INT_MAX;
            for(int j=0; j < num_fac; j++){
                int idxFac = (rand() % (num_fac - j)) + j;
                int fac = indexFac[idxFac];
                indexFac[idxFac] = indexFac[j];
                indexFac[j] = fac;

                
                if(vet_capacidade[fac] - s.vetorConsumoFac[fac] <= 0
                || (s.matriz[cli][0] != -1 && vet_demanda[cli] - s.matriz[cli][2] > vet_capacidade[fac] - s.vetorConsumoFac[fac])
                || s.matriz[cli][0] == fac) continue;

                eCompativel(s, cli, fac, compativel);
                if(!compativel) continue;


                custo = mat_custo[cli][fac];
                if(s.vetorConsumoFac[fac] == 0){
                    custo += fatorAmortizacao * vet_custo_fixo[fac] * MIN(vet_demanda[cli] - s.matriz[cli][2], vet_capacidade[fac] - s.vetorConsumoFac[fac]) / vet_capacidade[fac];
                }

                if(custo < melhorCusto){
                    melhorCusto = custo;
                    melhorCli = cli;
                    melhorFac = fac;

                    contMelhor = 1;
                    achou = true;
                }else if(custo < melhorCusto + janelaDeAceitacao){
                    contMelhor++;

                    if((rand() % 1000)/(double) 1000 <= (1 / (double) contMelhor)){
                        melhorCli = cli;
                        melhorFac = fac;

                        achou = true;
                    }
                }
            }


            if(achou){
                quantidade = MIN(
                                    vet_demanda[melhorCli] - (s.matriz[melhorCli][2] + s.matriz[melhorCli][3]),
                                    vet_capacidade[melhorFac] - s.vetorConsumoFac[melhorFac]
                                );
                
                inserirAbastecimento(s,melhorCli,melhorFac,quantidade);
                cargaAtual += quantidade;
            }else{
                //Escolhe um aleatório
            }
        }
    }

    calcular_fo(s);
}

void calcularCusto(int& cliente, int& facilidade, int& quantidade, int& custo){
    custo = mat_custo[cliente][facilidade] * quantidade;
}

void criarVetorDeIndex(int vetor[],int& tamanho){
    for(int i=0; i<tamanho; i++){
        vetor[i] = i;
    }
}

void eCompativel(Solucao& s, int& cliente, int& facilidade, bool& compativel){
    for(int i=0; i < vet_aux_incompat[cliente]; i++){
        int cliIncompativel = mat_incompat[cliente][i];

        if(s.matriz[cliIncompativel][0] == facilidade || s.matriz[cliIncompativel][1] == facilidade){
            compativel = false;
            return;
        }
    }

    compativel = true;
}

void montarMatrizSolucao(Solucao& s){
    for(int i=0; i < num_cli; i++){
        s.matriz[i][0] = s.matriz[i][1] = -1;
        s.matriz[i][2] = s.matriz[i][3] = 0;
    }
}

void inserirAbastecimento(Solucao& s, int& cliente, int& facilidade, int& q){
    if(s.matriz[cliente][0] == -1){
        s.matriz[cliente][0] = facilidade;
        s.matriz[cliente][2] = q;
    }else{
        s.matriz[cliente][1] = facilidade;
        s.matriz[cliente][3] = q;


        if(mat_custo[cliente][s.matriz[cliente][0]] > mat_custo[cliente][s.matriz[cliente][1]]){
            int temp;

            temp = s.matriz[cliente][0];
            s.matriz[cliente][0] = s.matriz[cliente][1];
            s.matriz[cliente][1] = temp;

            temp = s.matriz[cliente][2];
            s.matriz[cliente][2] = s.matriz[cliente][3];
            s.matriz[cliente][3] = temp;
        }
    }

    s.vetorConsumoFac[facilidade] += q;

}
