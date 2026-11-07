#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "estruturas.h"

#define BETA 1000

#define MAX(X,Y) ((X > Y) ? X : Y)
#define MIN(X,Y) ((X < Y) ? X : Y)

const int lambda = /*8*/100 ;
const int beta = 1 ;


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

void ler_dados(char* arq, int& num_pref){
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

    num_pref = BETA * ((int)sqrt(num_fac));
}

void criar_mat_pref(int& num_pref){
    int vet_temp[num_fac]; //Vetor cópia dos custos
    int vet_index[num_fac]; //Vetor de índices (referente ao vet_temp)

    //Inicia vetores auxiliares com valor 0
    memset(&vet_aux_fac_pref,0,sizeof(vet_aux_fac_pref));
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
                mat_fac_pref[x][vet_aux_fac_pref[x]] = vet_index[y];
                vet_aux_fac_pref[x]++;
                
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
        for(int idx_fac1=0; idx_fac1<vet_aux_fac_pref[cli]; idx_fac1++){
            fac1 = mat_fac_pref[cli][idx_fac1];
            for(int idx_fac2=idx_fac1 + 1; idx_fac2<vet_aux_fac_pref[cli]; idx_fac2++){
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

void calcular_fo(Solucao &s){
    s.fo = 0;
    int vetorCargaFac[num_fac];
    memset(&vetorCargaFac, 0, sizeof(vetorCargaFac));


    //Somatório do custo por carga
    for(int cli=0; cli < num_cli; cli++){
        s.fo += s.matriz[cli][2] * mat_custo[cli][s.matriz[cli][0]]; //Suponho que sempre haverá pelo menos uma facilidade
        vetorCargaFac[s.matriz[cli][0]] += s.matriz[cli][2];

        if(s.matriz[cli][1] != -1){
            s.fo += s.matriz[cli][3] * mat_custo[cli][s.matriz[cli][1]];
            vetorCargaFac[s.matriz[cli][1]] += s.matriz[cli][3];
        }
    }

    //Somatório de custo fixo
    for(int fac=0; fac < num_fac; fac++){
        s.fo += MIN(vetorCargaFac[fac], 1) * vet_custo_fixo[fac];

        //Penalidades
        s.fo += MAX(0, vetorCargaFac[fac] - vet_capacidade[fac]) * BETA;
    }
}
