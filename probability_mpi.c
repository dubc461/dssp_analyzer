#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define WORD_LENGTH 8
#define BASKET 10
#define COUNT 360/10

void cluster(int **sample_list, double **angle_list, int **natrual_distribution_index, int cluster_num, int line_max, int col_max, char *file_name){
  int init_row = 0, row = 0, num = 0;
  int i = 0, k = 0;
  char cluster_file[90] = {0};
  //printf("CLUSTER %d %d %d \n", cluster_num, line_max, col_max);
  for(i = 0; i < cluster_num; i++){
    sprintf(cluster_file, "%s.cluster_%d", file_name, i);
    FILE *fp;
    fp = fopen(cluster_file, "w");
    if(fp == NULL){
      printf("Cluster_Error %s %d %i\n", cluster_file, cluster_num, i);
      return;
    }
    num = 0;
    for(init_row = 0; init_row < line_max; init_row++){
      if(natrual_distribution_index[i][0] + 1 == sample_list[init_row][col_max]){
        for(k = 0; k < col_max; k++){
          fprintf(fp, "%lf ", angle_list[init_row][k]);
        }
        fprintf(fp, " \n");
        num = num + 1;
      }
/*
      else{
        printf("%d %d \n", natrual_distribution_index[i][0] + 1, sample_list[init_row][col_max]);
      }
*/
      if(num == natrual_distribution_index[i][1]){
        break;
      }
    }
    close(fp);
  }
  printf("Cluster Done! \n");
  return;
}

int array_compare(int **sample_list, int row_1, int row_2, int col, int col_max, int cluster){
  if(sample_list[row_1][col] == sample_list[row_2][col]) {
    if(col > 0){
      return array_compare(sample_list, row_1, row_2, col - 1, col_max, cluster);
    }
    else{
      sample_list[row_2][col_max] = cluster + 1;
      return 1;
    }
  }
  else{
    return 0;
  }
}

int natrual_distribution(int **sample_list, int **natrual_distribution_index, int line_max, int col_max){
  int init_row = 0, row = 0, num = 0;
  int i = 0;
  for(init_row = 0; init_row < line_max; init_row++){
    if(sample_list[init_row][col_max] != 0){
      continue;
    }
    for(row = init_row; row < line_max; row++){
      num = num + array_compare(sample_list, init_row, row, col_max - 1, col_max, init_row);
    }
    natrual_distribution_index[i][0] = init_row;
    natrual_distribution_index[i][1] = num;
    i = i + 1;
    num = 0;
  }
  return i;
}

double single_distribution(int **sample_list, int state_index, int col_A, int angle_A, int col_last, int line_max){
  int line = 0, sample_num = 0;
  for(line = 0; line < line_max; line++){
    if( state_index + 1 == sample_list[state_index][col_last] && sample_list[line][col_A] == angle_A ){
      sample_num = sample_num + 1;
    }
  }
  return sample_num/(double)line_max;
}

double joint_distrubution_2nd(int **sample_list, int state_index, int col_A, int col_B, int angle_A, int angle_B, int col_last, int line_max){
  int sample_num = 0;
  int line;
  for(line = 0; line < line_max; line++){
    if(state_index + 1 == sample_list[state_index][col_last] && sample_list[line][col_A] == angle_A && sample_list[line][col_B] == angle_B){
      sample_num = sample_num + 1;
    }
  }
  return sample_num/(double)line_max;
}

//P(x_n|x_n-1) = P(x_n,x_n-1)/P(x_n-1)
double conditional_probability_2nd(int **sample_list, int state_index, int col_A, int col_B, int angle_A, int angle_B, int col_last, int line_max){
  double probility = 0;
  probility = (double)joint_distrubution_2nd( sample_list, state_index, col_A, col_B, angle_A, angle_B, col_last, line_max) / single_distribution( sample_list, state_index, col_A, angle_A, col_last, line_max);
  return probility;
}

double joint_distrubution_3nd(int **sample_list, int state_index, int col_A, int col_B, int col_C, int angle_A, int angle_B, int angle_C, int col_last, int line_max){
  double sample_num = 0;
  int line;
  for(line = 0; line < line_max; line++){
    if(state_index + 1 == sample_list[state_index][col_last] && sample_list[line][col_A] == angle_A && sample_list[line][col_B] == angle_B && sample_list[line][col_C] == angle_C){
      sample_num = sample_num + 1;
    }
  }
  return sample_num/(double)line_max;
}

double conditional_probability_3nd(int **sample_list, int state_index, int col_A, int col_B, int col_C, int angle_A, int angle_B, int angle_C, int col_last, int line_max){
  double probility = 0, a = 0, b = 0;
  if(joint_distrubution_2nd(sample_list, state_index, col_A, col_B, angle_A, angle_B, col_last, line_max) == 0 ){
    probility = 0;
  }
  else{
    a = joint_distrubution_3nd(sample_list, state_index, col_A, col_B, col_C, angle_A, angle_B, angle_C, col_last, line_max);
    b = joint_distrubution_2nd(sample_list, state_index, col_A, col_B, angle_A, angle_B, col_last, line_max);
    probility = a / b;
  }
  return probility;
}

double joint_distrubution_4nd(int **sample_list, int state_index, int col_A, int col_B, int col_C, int col_D, int angle_A, int angle_B, int angle_C, int angle_D, int col_last, int line_max){
  double sample_num = 0;
  int line;
  for(line = 0; line < line_max; line++){
    if(state_index + 1 == sample_list[state_index][col_last] && sample_list[line][col_A] == angle_A && sample_list[line][col_B] == angle_B && sample_list[line][col_C] == angle_C && sample_list[line][col_D] == angle_D ){
      sample_num = sample_num + 1;
    }
  }
  return sample_num/(double)line_max;
}

double conditional_probability_4nd(int **sample_list, int state_index, int col_A, int col_B, int col_C, int col_D, int angle_A, int angle_B, int angle_C, int angle_D, int col_last, int line_max ){
  double probility = 0;
  if(joint_distrubution_3nd(sample_list, state_index, col_B, col_C, col_D, angle_B, angle_C, angle_D, col_last, line_max) == 0){
    probility = 0;
  }
  else{
    probility = joint_distrubution_4nd(sample_list, state_index, col_A, col_B, col_C, col_D, angle_A, angle_B, angle_C, angle_D, col_last, line_max ) / joint_distrubution_3nd(sample_list, state_index, col_A, col_B, col_C, angle_A, angle_B, angle_C, col_last, line_max);
  }
  return probility;
}


int main(int argc,char *argv[]){
  
  int comm_sz;
  int my_rank = 0;
  MPI_Status mystatus;
  
  int line_length = 0;
  int file_length = 0;
  int line = 0;
  int cols = 0;
  int i = 0, j = 0, k = 0;
  int position = 0;
  
  int single_argument_buffer[4] = {0};
  int size_argument_buffer = sizeof(int) * 4;
  int size_buffer_double = 0;
  int package_length = 0;


  MPI_Init(&argc, &argv);
  //MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  if(my_rank == 0){
    FILE *fp ,*wp;
    int **sample_list;
    double **angle_list;
    static int **natrual_distribution_index;
    
    fp = fopen(argv[1], "r");
    //fp = fopen("/home/dubc/SS_01/H_2-0.dat", "r");
    if (fp == NULL) {
      printf(" Do not exist\n");
      exit(1);
    }
    
    while(fgetc(fp) != '\n'){
      line_length = (int) ftell(fp) + 1;
    }
    fseek(fp, SEEK_SET, SEEK_END);
    file_length = (int) ftell(fp);
    rewind(fp);
    
    line = file_length/line_length;
    cols = line_length/WORD_LENGTH;
    
    sample_list =  (int **) malloc(sizeof(int *) * file_length/line_length);
    angle_list =  (double **) malloc(sizeof(double *) * file_length/line_length);
    
    for (j = 0; j < file_length/line_length; j++){
      sample_list[j] = (int *) malloc(sizeof(int) * (line_length/WORD_LENGTH + 1));
      for (i = 0; i < line_length/WORD_LENGTH + 1; i++){
        sample_list[j][i] = 0;
      }
    }
    
    for (j = 0; j < file_length/line_length; j++){
      angle_list[j] = (double *) malloc(sizeof(double) * (line_length/WORD_LENGTH + 1));
      for (i = 0; i < line_length/WORD_LENGTH + 1; i++){
        angle_list[j][i] = 0;
      }
    }
    
    natrual_distribution_index =  (int **) malloc(sizeof(int *) * file_length/line_length);
    for (j = 0; j < file_length/line_length; j++){
      natrual_distribution_index[j] = (int *) malloc(sizeof(int) * line_length/WORD_LENGTH);
      for (i = 0; i < line_length/WORD_LENGTH; i++){
        natrual_distribution_index[j][i] = 0;
      }
    }
    
    i = 0;
    int num_break = 0;
    double sample = 0;
    while( 1 ) {
      fscanf(fp, "%lf", &sample);
      if (ftell(fp) == file_length) {
        break;
      }
      if(sample == 360 || sample == -360){
        sample = 0;
      }
      sample_list[i][num_break] = (int) ((sample + 180) / (COUNT) + 0.5);
      angle_list[i][num_break] = sample;
      num_break = num_break + 1;
      if(num_break % (line_length / WORD_LENGTH) == 0){
        i = i + 1;
        num_break = 0;
      }
    }
    fclose(fp);
    printf("Reading Done! %d \n", i);
    
    //PACKAGE
    int package_length = line / (comm_sz - 1);
    size_buffer_double = sizeof(double) * cols * package_length ;
    double *line_buffer = NULL;
    
    line_buffer =  (double *) malloc(sizeof(double) * size_buffer_double);
    printf("%d %d %d %d \n", line, cols, package_length, comm_sz);
    MPI_Pack( &line, 1, MPI_INT, single_argument_buffer, size_argument_buffer, &position, MPI_COMM_WORLD );
    MPI_Pack( &cols, 1, MPI_INT, single_argument_buffer, size_argument_buffer, &position, MPI_COMM_WORLD );
    MPI_Pack( &package_length, 1, MPI_INT, single_argument_buffer, size_argument_buffer, &position, MPI_COMM_WORLD );
    MPI_Pack( &size_buffer_double , 1, MPI_INT, single_argument_buffer, size_argument_buffer, &position, MPI_COMM_WORLD );
    
    MPI_Bcast( single_argument_buffer, size_argument_buffer, MPI_PACKED, 0, MPI_COMM_WORLD);
    position = 0;
    
    for(i = 1; i < comm_sz; i++){
      for(j = ( i - 1 ) * package_length; j <  i * package_length; j++) {
        MPI_Pack( angle_list[j], cols, MPI_DOUBLE, line_buffer, size_buffer_double, &position, MPI_COMM_WORLD );
      }
      MPI_Send( line_buffer, size_buffer_double, MPI_PACKED, i, 0, MPI_COMM_WORLD );
      position = 0;
    }


//
    int line_left = line - (comm_sz - 1) * package_length;
    //printf("LEFT:%d \n", line_left);

    int **sample_list_left =  (int **) malloc(sizeof(int *) * line_left);
    for (j = 0; j < line_left; j++){
      sample_list_left[j] = (int *) malloc(sizeof(int) * (cols + 1));
      for (i = 0; i < cols + 1; i++){
        if( i != cols ){
          sample_list_left[j][i] = sample_list[j + (line - line_left)][i];
          //printf("K:%lf ", angle_list[j + (line - line_left)][i]);
        }
        else{
          sample_list_left[j][i] = 0;
        }
      }
      //printf("\n");
    }

    int **natrual_distribution_indexs_left = NULL;
    natrual_distribution_indexs_left =  (int **) malloc(sizeof(int *) * (line - line_left));
    for (j = 0; j < (line - line_left); j++){
      natrual_distribution_indexs_left[j] = (int *) malloc(sizeof(int) * 2);
      for (i = 0; i < 2; i++){
        natrual_distribution_indexs_left[j][i] = 0;
      }
    }
    if(line_left > 1 ){ 
      natrual_distribution(sample_list_left, natrual_distribution_indexs_left, 1, 0);
    }
    else{
      natrual_distribution_indexs_left[0][0] = line - line_left;
      natrual_distribution_indexs_left[0][1] = 0;
    }
//
    int ND_index = 0;
    int ND_index_sum = 0;
    int ND_I_size = 0;
    int *ND_I_line_buffer = NULL;
    int *ND_I_single_line = (int *) malloc(sizeof(int) * 2);
    int **ND_I_check = (int **) malloc(sizeof(int *) * line);
    for(i = 0; i < line; i++){
      ND_I_check[i] = (int *) malloc(sizeof(int) * 2);
      for(j = 0; j < 2; j++){
        ND_I_check[i][j] = 0;
      } 
    }

    for(i = 1; i < comm_sz; i++){
      MPI_Recv( &ND_index, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &mystatus );
    //  printf("%d %d \n", i, ND_index);
      ND_I_size = sizeof(int) * 2 * ND_index;
      ND_I_line_buffer = (int *) malloc(sizeof(int) * ND_I_size);

      MPI_Recv( ND_I_line_buffer, ND_I_size, MPI_PACKED, i, 0, MPI_COMM_WORLD, &mystatus );
      position = 0;
      for(j = 0; j < ND_index; j++) {
        MPI_Unpack( ND_I_line_buffer, ND_I_size, &position, ND_I_single_line, 2, MPI_INT, MPI_COMM_WORLD );
        ND_I_check[j + ND_index_sum][0] = ND_I_single_line[0] + ( i - 1 ) * package_length;
        ND_I_check[j + ND_index_sum][1] = ND_I_single_line[1];
        //printf("%d %d \n", ND_I_check[j + ND_index_sum][0], ND_I_check[j + ND_index_sum][1]);
      }
      position = 0;
      ND_index_sum = ND_index_sum + ND_index;
      //printf("SUM:%d \n", ND_index_sum);
    }

    int *cluster_array = (int *) malloc(sizeof(int) * line);
    int *cluster_array_part = (int *) malloc(sizeof(int) * package_length);
    int *cluster_array_buffer = (int *) malloc(sizeof(int) * package_length);
    int cluster_array_buffer_size = sizeof(int) * package_length;
    for (i = 0; i < package_length; i++){
      cluster_array_buffer[i] = 0;
    }
    for (i = 0; i < line; i++){
      cluster_array[i] = 0;
    }

    for(i = 1; i < comm_sz; i++){
      MPI_Recv( cluster_array_buffer, cluster_array_buffer_size, MPI_PACKED, i, 0, MPI_COMM_WORLD, &mystatus );
      position = 0;
      MPI_Unpack( cluster_array_buffer, cluster_array_buffer_size, &position, cluster_array_part, package_length, MPI_INT, MPI_COMM_WORLD );
      for(j = 0; j < package_length; j++) {
        cluster_array[ ( i - 1 ) * package_length + j ] = cluster_array_part[j] + ( i - 1 ) * package_length;
//        printf("%d %d %d\n", ( i - 1 ) * package_length + j, cluster_array_part[j], cluster_array[ ( i - 1 ) * package_length + j ]);
      }
      position = 0;
    }    

    for(i = 0; i < ND_index_sum; i++) {
      for(j = 0; j < ND_index_sum; j++) {
        if(i == j || ND_I_check[j][0] == -1 || ND_I_check[i][0] == -1) {
          continue;
        }
        if( array_compare(sample_list, ND_I_check[i][0], ND_I_check[j][0], cols - 1, cols, ND_I_check[i][0]) == 1 ) {
          for( k = 0; k < line; k++){
            if(cluster_array[k] - 1 == ND_I_check[j][0]) {
              cluster_array[k] = ND_I_check[i][0] + 1;
            }
          }
          ND_I_check[j][0] = -1;
          ND_I_check[i][1] = ND_I_check[i][1] + ND_I_check[j][1];
        }
      }
    }

   for(i = 0; i < line; i++){
     sample_list[i][cols] = cluster_array[i];
//     printf("%d %d \n", i, cluster_array[i]);
    }

    char outfile[50] = {0};
    sprintf(outfile, "%s.prob", argv[1]);
    wp = fopen(outfile, "w");

    int ND_I = 0;
    for(i = 0; i < ND_index_sum; i++) {
      if(ND_I_check[i][0] != -1) {
        //printf("%d %d %d \n", i, ND_I_check[i][0], ND_I_check[i][1]); 
        natrual_distribution_index[ND_I][0] = ND_I_check[i][0];
        natrual_distribution_index[ND_I][1] = ND_I_check[i][1];
        ND_I = ND_I + 1;
      }
    }
    cluster(sample_list, angle_list, natrual_distribution_index, ND_I, line, cols, outfile);
  //  printf("%s Seems Done! \n", argv[1]);
   
   int state_index = 0, state_count = 0;
   double eq1 = 0, eq2 = 1, eq3 = 1, eq4 = 1;
   for(i=0;i<ND_I;i++){
   state_index = natrual_distribution_index[i][0];
   state_count = natrual_distribution_index[i][1];
   
   //eq1
   for(k = 0; k < line_length/WORD_LENGTH; k++){
   if(k == 0){
   eq1 = single_distribution(sample_list, state_index, k, sample_list[state_index][k], line_length/WORD_LENGTH, file_length/line_length);
   }
   else{
   eq1 = eq1 * single_distribution(sample_list, state_index, k, sample_list[state_index][k], line_length/WORD_LENGTH, file_length/line_length);
   }
   }
   
   //eq2
   for(k = 0; k < line_length/WORD_LENGTH - 1; k++){
   eq2 = eq2 * conditional_probability_2nd(sample_list, state_index,  k, k+1, sample_list[state_index][k], sample_list[state_index][k+1], line_length/WORD_LENGTH, file_length/line_length);
   }
   
   //eq3
   if(line_length/WORD_LENGTH-2 > 0){
   for(k=0;k<line_length/WORD_LENGTH-2;k++){
   eq3 = eq3 * conditional_probability_3nd(sample_list, state_index, k, k+1, k+2, sample_list[state_index][k], sample_list[state_index][k+1], sample_list[state_index][k+2], line_length/WORD_LENGTH, file_length/line_length);
   }
   }
   else{
   eq3 = 0;
   }
   
   //eq4
   if(line_length/WORD_LENGTH-3 > 0){
   for(k=0;k<line_length/WORD_LENGTH-3;k++){
   eq4 = eq4 * conditional_probability_4nd(sample_list, state_index, k, k+1, k+2, k+3, sample_list[state_index][k], sample_list[state_index][k+1], sample_list[state_index][k+2], sample_list[state_index][k+3], line_length/WORD_LENGTH, file_length/line_length);
   }
   }
   else{
   eq4 = 0;
   }
   
   fprintf(wp, "%lf %lf %lf %lf %lf\n", (double)state_count/(double)(file_length/line_length), eq1, eq2, eq3, eq4);
   }
   fclose(wp);
   
   }

  else{
    MPI_Bcast( single_argument_buffer, size_argument_buffer, MPI_PACKED, 0, MPI_COMM_WORLD);
    position = 0;
    MPI_Unpack(single_argument_buffer, size_argument_buffer, &position, &line, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(single_argument_buffer, size_argument_buffer, &position, &cols, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(single_argument_buffer, size_argument_buffer, &position, &package_length, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(single_argument_buffer, size_argument_buffer, &position, &size_buffer_double, 1, MPI_INT, MPI_COMM_WORLD);
//    printf("%d %d %d %d \n", line, cols, package_length, size_buffer_double);
    
    double *lines = NULL;
    lines =  (double *) malloc(sizeof(double) * size_buffer_double);
    
    double *single_line_buffer = NULL;
    single_line_buffer =  (double *) malloc(sizeof(double) * cols);
    
    double **lines_array = NULL;
    lines_array =  (double **) malloc(sizeof(double *) * package_length);
    for (j = 0; j < package_length; j++){
      lines_array[j] = (double *) malloc(sizeof(double) * cols);
      for (i = 0; i < cols; i++){
        lines_array[j][i] = 0;
      }
    }
    
    int **lines_array_sample = NULL;
    lines_array_sample = (int **) malloc(sizeof(int *) * package_length);
    
    MPI_Recv( lines, size_buffer_double, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &mystatus);
    
    position = 0;
    for(j = 0; j < package_length; j++) {
      MPI_Unpack(lines, size_buffer_double, &position, single_line_buffer, cols, MPI_DOUBLE, MPI_COMM_WORLD);
      //      memcpy ( lines_array[j], single_line_buffer, sizeof(single_line_buffer) );
      for(i = 0; i < cols; i++){
        lines_array[j][i] = single_line_buffer[i];
      }
      //printf("%lf %lf %lf %lf\n", single_line_buffer[0], single_line_buffer[1], lines_array[j][0], lines_array[j][1]);
    }
    for (j = 0; j < package_length; j++){
      lines_array_sample[j] = (int *) malloc(sizeof(int) * (cols + 1) );
      for (i = 0; i < cols; i++){
        lines_array_sample[j][i] = (int) ((lines_array[j][i] + 180) / (COUNT) + 0.5);
//        printf("%d ", lines_array_sample[j][i]);
      }
//      printf("\n");
    }
    
    int ND_I = 0;
    int **natrual_distribution_indexs = NULL;
    natrual_distribution_indexs =  (int **) malloc(sizeof(int *) * package_length);
    
    for (j = 0; j < package_length; j++){
      natrual_distribution_indexs[j] = (int *) malloc(sizeof(int) * 2);
      for (i = 0; i < 2; i++){
        natrual_distribution_indexs[j][i] = 0;
      }
    }
    
    ND_I = natrual_distribution(lines_array_sample, natrual_distribution_indexs, package_length, cols);
    MPI_Send( &ND_I, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );

    //cluster(lines_array_sample, lines_array, natrual_distribution_index, ND_I, line, cols, outfile);
//    printf("%d \n", ND_I);

    int ND_I_size = sizeof(int) * 2 * ND_I;
    int *ND_I_line_buffer = (int *) malloc(sizeof(int) * ND_I_size);
 
    position = 0;
    for(j = 0; j < package_length; j++) {
      MPI_Pack( natrual_distribution_indexs[j], 2, MPI_INT, ND_I_line_buffer, size_buffer_double, &position, MPI_COMM_WORLD );
    }
    MPI_Send( ND_I_line_buffer, ND_I_size, MPI_PACKED, 0, 0, MPI_COMM_WORLD );

    int *cluster_array = NULL;
    int *cluster_array_buffer = NULL;
    int cluster_array_buffer_size = sizeof(int) * package_length;

    cluster_array = (int *) malloc(sizeof(int) * package_length);
    cluster_array_buffer = (int *) malloc(sizeof(int) * package_length);
    for (i = 0; i < package_length; i++){
      cluster_array[i] = lines_array_sample[i][cols];
      cluster_array_buffer[i] = 0;
    }

    position = 0;
    MPI_Pack( cluster_array, package_length, MPI_INT, cluster_array_buffer, cluster_array_buffer_size, &position, MPI_COMM_WORLD );
    position = 0;
    MPI_Send( cluster_array_buffer, cluster_array_buffer_size, MPI_PACKED, 0, 0, MPI_COMM_WORLD );       
  }

  MPI_Finalize();
    printf("%s Seems Done! \n", argv[1]);

  return 0;
}
