#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <climits>
#include "header.h"

#define MAX(X,Y) ((X > Y) ? X : Y)
#define MIN(X,Y) ((X < Y) ? X : Y)

#define BETA 1000

//const char arq[] = "../Instancias/cflp-ci_toy.txt";
const char arq[] = "../Instancias/wlp01.txt";

int main(){
    srand(time(NULL));
    ler_dados();
    //testar_dados("");

    Solucao s,temp;
    int tamanho = num_cli * 25 / 100;
    construtivaGulosa(s);
    printf("FO: %i\n",s.fo);
    //imprimirSolucao(s);

    while(true){
        memcpy(&temp, &s, sizeof(Solucao));
        while(true){
            gerarVizinhoMineReduced(temp, tamanho);

            if(temp.fo < s.fo){
                printf("FO: %i\n", temp.fo);
                memcpy(&s, &temp, sizeof(Solucao));
                break;
            }
        }
    }    
    
    printf("FO: %i\n", s.fo);
    //imprimirSolucao(s);
    return 0;
}

void gerarVizinhoMineReduced(Solucao& s, int& tamanho){ 
    Solucao temp;
    memcpy(&temp, &s, sizeof(Solucao));

    int vetorIndex[tamanho], vetorIndexGlobal[num_cli];
    int cliente, facilidade, fac1, fac2;
    int index,indexFac;
    int vetorFac[num_fac],vetorBinarioFac[num_fac];
    int auxVetorFac=0;
    
    //Fase de Seleção de Clientes
    for(int x=0; x < num_cli; x++) vetorIndexGlobal[x] = x;
    for(int idx=0; idx < tamanho; idx++){
        index = rand() % (num_cli - idx) + idx;
        cliente = vetorIndexGlobal[index];
        vetorIndexGlobal[index] = vetorIndexGlobal[idx];
        vetorIndexGlobal[idx] = cliente;

        vetorIndex[idx] = cliente;
    }

    // Fase de Desconexão
    memset(&vetorBinarioFac, 0, sizeof(vetorBinarioFac));
    for(int cli=0; cli < tamanho; cli++){
        cliente = vetorIndex[cli];
        fac1 = temp.matriz[cliente][0];
        fac2 = temp.matriz[cliente][1];


        if(fac2 != -1){
            temp.vetorConsumoFac[fac2] -= temp.matriz[cliente][3]; //Retira carga da facilidade
            temp.fo -= mat_custo[cliente][fac2] * temp.matriz[cliente][3]; //retira o valor da FO
            temp.matriz[cliente][1] = -1; //Retira facilidade de cliente
            temp.matriz[cliente][3] = 0; //Retira carga da facilidade de cliente

            if(!vetorBinarioFac[fac2]){ //Caso não contenha
                vetorFac[auxVetorFac] = fac2;
                auxVetorFac++;
                vetorBinarioFac[fac2] = 1;
            }
        }

        if(fac1 != -1){ //Talvez não precise, considerar remoção
            temp.vetorConsumoFac[fac1] -= temp.matriz[cliente][2]; //Retira carga da facilidade
            temp.fo -= mat_custo[cliente][fac1] * temp.matriz[cliente][2]; //retira o valor da FO
            temp.matriz[cliente][0] = -1; //Retira facilidade de cliente
            temp.matriz[cliente][2] = 0; //Retira carga da facilidade de cliente

            if(!vetorBinarioFac[fac1]){ //Caso não contenha
                vetorFac[auxVetorFac] = fac1;
                auxVetorFac++;
                vetorBinarioFac[fac1] = 1;
            }
        }
    }
    for(int x=0; x < auxVetorFac; x++){
        if(temp.vetorConsumoFac[vetorFac[x]] == 0){
            temp.fo -= vet_custo_fixo[vetorFac[x]];
        }
    }

    // Fase de Reconexão
    bool achou = false;
    bool compativel=false;
    int quantidade;
    int custo, melhorCusto, melhorFac, facAleatoria;
    int cargaTotal=0, cargaAtual=0;
    int vetorIndexFac[num_fac];
    for(int x=0; x<num_fac; x++) vetorIndexFac[x] = x;
    for(int x=0; x<num_cli;x++)cargaTotal += vet_demanda[x];
    for(int x=0; x<num_fac; x++)cargaAtual += temp.vetorConsumoFac[x];

    while(cargaAtual < cargaTotal){

        for(int i=0; i < tamanho; i++){
            index = rand() % (tamanho - i) + i;
            cliente = vetorIndex[index];
            vetorIndex[index] = vetorIndex[i];
            vetorIndex[i] = cliente;

            if(vet_demanda[cliente] - (temp.matriz[cliente][2] + temp.matriz[cliente][3]) <= 0) continue;
    
            achou = false;
            melhorCusto = INT_MAX;
            for(int j=0; j < auxVetorFac; j++){
                indexFac = rand() % (auxVetorFac - j) + j;
                facilidade = vetorFac[indexFac];
                vetorFac[indexFac] = vetorFac[j];
                vetorFac[j] = facilidade;
    
                if(vet_capacidade[facilidade] - temp.vetorConsumoFac[facilidade] <= 0
                || (temp.matriz[cliente][0] != -1 && vet_demanda[cliente] - temp.matriz[cliente][2] > vet_capacidade[facilidade] - temp.vetorConsumoFac[facilidade])
                || temp.matriz[cliente][0] == facilidade) continue;

                eCompativel(temp, cliente, facilidade, compativel);
                if(!compativel) continue;

                custo = mat_custo[cliente][facilidade];

                if(custo < melhorCusto){
                    melhorCusto = custo;
                    melhorFac = facilidade;

                    achou = true;
                }
            }

            if(achou){
                quantidade = MIN(
                                    vet_demanda[cliente] - (temp.matriz[cliente][2] + temp.matriz[cliente][3]),
                                    vet_capacidade[melhorFac] - temp.vetorConsumoFac[melhorFac]
                                );
                
                inserirAbastecimento(temp,cliente,melhorFac,quantidade);
                cargaAtual += quantidade;
            }else{
                for(int j=0; j < num_fac; j++){
                    indexFac = rand() % (num_fac - j) + j;
                    facilidade = vetorIndexFac[indexFac];
                    vetorIndexFac[indexFac] = vetorIndexFac[j];
                    vetorIndexFac[j] = facilidade;

                    if(facilidade == temp.matriz[cliente][0]) continue;

                    eCompativel(temp, cliente, facilidade, compativel);
                    if(!compativel) continue;

                    quantidade = vet_demanda[cliente] - (temp.matriz[cliente][2] + temp.matriz[cliente][3]);
                    inserirAbastecimento(temp, cliente, facilidade, quantidade);

                    temp.fo += BETA * MAX(0, temp.vetorConsumoFac[facilidade] - vet_capacidade[facilidade]);

                    cargaAtual += quantidade;
                    achou = true;
                    break;
                }
                if(!achou){
                    return;
                }
            }
        }
    }

    // Facilidades Fechadas (remoção do custo de abertura)
    for(int x=0; x < auxVetorFac; x++){
        if(temp.vetorConsumoFac[vetorFac[x]] == 0){
            temp.fo -= vet_custo_fixo[vetorFac[x]];
        }
    }

    memcpy(&s, &temp, sizeof(Solucao));
}

void construtivaAleatoriaGulosa(Solucao& s){
    int contMelhor;
    int melhorCli, melhorFac;
    int quantidade;
    int facAleatoria;
    double custo, melhorCusto;
    bool compativel = false, achou = false;
    s.fo = 0;

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

            if(vet_demanda[cli] - (s.matriz[cli][2] + s.matriz[cli][3]) <= 0) continue;

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
                do{
                    facAleatoria = rand() % num_fac;
                    eCompativel(s, cli, facAleatoria, compativel);
                }while(!compativel || facAleatoria == s.matriz[cli][0]);
                
                quantidade = vet_demanda[cli] - (s.matriz[cli][2] + s.matriz[cli][3]);
                inserirAbastecimento(s, cli, facAleatoria, quantidade);
                cargaAtual += quantidade;
            }
        }
    }

    calcular_fo(s);
}

void construtivaGulosa(Solucao& s){
    int nFalha = 0;
    int melhorCli, melhorFac;
    int quantidade;
    double custo, melhorCusto;
    bool compativel = false, achou = false;
    s.fo = 0;

    int cargaTotal=0, cargaAtual=0;
    for(int x=0;x<num_cli;x++)cargaTotal += vet_demanda[x];
    

    int indexCli[num_cli], indexFac[num_fac];
    criarVetorDeIndex(indexCli,num_cli);
    criarVetorDeIndex(indexFac, num_fac);

    memset(s.vetorConsumoFac, 0, sizeof(s.vetorConsumoFac));
    montarMatrizSolucao(s);

    while(cargaAtual < cargaTotal){
        achou = false;
        melhorCusto = INT_MAX;

        for(int cli=0; cli < num_cli; cli++){

            if(vet_demanda[cli] - (s.matriz[cli][2] + s.matriz[cli][3]) <= 0) continue; 

            
            for(int fac=0; fac < num_fac; fac++){
                
                if((s.matriz[cli][0] != -1 && vet_demanda[cli] - s.matriz[cli][2] > vet_capacidade[fac] - s.vetorConsumoFac[fac])
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
            break;
        }
    }

    if(!achou){
        for(int cli=0; cli<num_cli; cli++){
            if(vet_demanda[cli] < s.matriz[cli][2] + s.matriz[cli][3]) continue;
            int facAleatoria;

            do{
                facAleatoria = rand() % num_fac;
                eCompativel(s, cli, facAleatoria, compativel);
            }while(!compativel || facAleatoria == s.matriz[cli][0]);
            
            quantidade = vet_demanda[cli] - (s.matriz[cli][2] + s.matriz[cli][3]);
            inserirAbastecimento(s, cli, facAleatoria, quantidade);
        }
    }
    
    calcular_fo(s);
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

    num_pref = BETA * ((int)sqrt(num_fac));
}

void criar_mat_pref(){
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
    memset(&s.vetorConsumoFac, 0, sizeof(s.vetorConsumoFac));


    //Somatório do custo por carga
    for(int cli=0; cli < num_cli; cli++){
        s.fo += s.matriz[cli][2] * mat_custo[cli][s.matriz[cli][0]]; //Suponho que sempre haverá pelo menos uma facilidade
        s.vetorConsumoFac[s.matriz[cli][0]] += s.matriz[cli][2];

        if(s.matriz[cli][1] != -1){
            s.fo += s.matriz[cli][3] * mat_custo[cli][s.matriz[cli][1]];
            s.vetorConsumoFac[s.matriz[cli][1]] += s.matriz[cli][3];
        }
    }

    //Somatório de custo fixo
    for(int fac=0; fac < num_fac; fac++){
        s.fo += MIN(s.vetorConsumoFac[fac], 1) * vet_custo_fixo[fac];

        //Penalidades
        s.fo += MAX(0, s.vetorConsumoFac[fac] - vet_capacidade[fac]) * BETA;
    }
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
    if(s.vetorConsumoFac[facilidade] == 0){
        s.fo += vet_custo_fixo[facilidade];
    }

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
    s.fo += mat_custo[cliente][facilidade] * q;

}

void imprimirSolucao(Solucao& s){
    printf("\nMATRIZ:\n");
    for(int x=0; x<num_cli; x++){
        printf("%i\t%i\t%i\t%i\n",s.matriz[x][0],s.matriz[x][1],s.matriz[x][2],s.matriz[x][3]);
    }

    printf("\nCONSUMO:\n");
    for(int x=0; x<num_fac; x++){
        printf("%i ",s.vetorConsumoFac[x]);
    }

    printf("\nFO: %i\n",s.fo);
}