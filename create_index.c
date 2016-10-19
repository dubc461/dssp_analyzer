//
//  create_index.c
//  dssp_mainchain_dihedral
//
//  Created by DuBochuan on 16/6/30.
//  Copyright © 2016年 DuBochuan. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define ROW 10000
#define PDB_INDEX 108682

int main(int argc, char *argv[]) {
  //PDB_INDEX
  FILE *pdb_index_file;
  if((pdb_index_file = fopen(argv[1],"r")) == NULL) {
    printf("No file\n");
    exit(1);
  }

  char pdb_id[50],pdb_file[50], pdb_index[PDB_INDEX][60];
  int i_pdb=0;

  while (!feof(pdb_index_file)) {
    fscanf(pdb_index_file, "%s", pdb_id);
    snprintf(pdb_index[i_pdb], sizeof(pdb_index[i_pdb]), "/data17/DSSP/%s.dssp", pdb_id);
    //printf("%s %s %d\n", pdb_index[i_pdb], pdb_id, i_pdb);
    i_pdb = i_pdb + 1;
  }
  fclose(pdb_index_file);

  int pdb_num=0;
  FILE *dssp_file;
  for (pdb_num=0; pdb_num < i_pdb; pdb_num++) {
    if((dssp_file = fopen(pdb_index[pdb_num], "r")) == NULL) {
      continue;
      printf("%s:No file\n",pdb_index[pdb_num]);
    }
    else{
      //跳过摘要
      int i=128;
      char line_break = '\n';
      while (line_break == '\n') {
        fseek(dssp_file,i,SEEK_SET);
        line_break = getc(dssp_file);
        i=i+(128+1);
      }

      i=i-(128+1)-128;
      while (line_break != '\n') {
        fseek(dssp_file,i,SEEK_SET);
        line_break = getc(dssp_file);
        i=i+1;
      }

      //DSSP_Record
      fseek(dssp_file,i,SEEK_SET);
      char line[ROW][150];
      int rows=0;
      while (!feof(dssp_file)) {
        fgets(line[rows], 138, dssp_file);
        rows=rows+1;
      }
      fclose(dssp_file);

      //Create_Index
      int line_max=0;
      char ss_type=0,chain=0;
      int residue_start=0,residue_end=0,residue_length=0;
      line_max=rows;
      rows=0;
      for (rows=0; rows < line_max; rows++) {
        if (rows==0) {
          ss_type = line[rows][16];
          chain = line[rows][11];
          residue_start = rows+1;
        }
        else if (ss_type != line[rows][16]) {
          residue_end = rows + 1;
          residue_length = residue_end - residue_start;
          if (ss_type == ' ') {
            ss_type='Z';
          }
          strncpy(pdb_file, pdb_index[pdb_num] + 13, 4);
          pdb_file[4] = '\0';
          printf("%s,%4d,%3d,%1c,%1c\n",pdb_file,residue_start,residue_length,ss_type,chain);
          ss_type = line[rows][16];
          chain = line[rows][11];
          residue_start = residue_end;
        }
      }
    }
  }
  return 0;
}

