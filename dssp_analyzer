//
//  dssp_analyzer.c
//  dssp_analyzer
//
//  Created by DuBochuan on 16/7/16.
//  Copyright © 2016年 DuBochuan. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PDB_INDEX 2000
#define Single_SS_Max 7506570

typedef struct ss {
  char pdb_id[5];
  int pdb_residue;
  int ss_length;
  char ss_type;
  char ss_chain;
}ss;

int main(int argc, char** argv) {
  char search_type = '0';
  int search_length = 0;
  int search_angle_type = 0;
  if (argc == 1) {
    printf("%s ss_type length fist_angle_type(0,1)\n",argv[0]);
    exit(1);
  }
  else{
    search_type = *argv[1];
    search_length = atoi(argv[2]);
    search_angle_type = atoi(argv[3]);
    //printf("%c,%d\n",search_type,search_length);
  }

  FILE *dssp_index;
  if((dssp_index = fopen("pdb_ss_index.txt","r")) == NULL) {
    printf("No file\n");
    exit(1);
  }

  char line[30];
  char pdb_id[5],pdb_residue_c[5],ss_length_c[5];
  char ss_type,ss_chain;
  int num=0;

  static ss ss_list[Single_SS_Max];
  while (fgets(line, 19, dssp_index)) {
    strncpy(pdb_id, line, 4);
    pdb_id[4]='\0';
    strncpy(pdb_residue_c, line+5, 4);
    pdb_residue_c[4]='\0';
    strncpy(ss_length_c, line+10, 3);
    ss_length_c[4]='\0';
    ss_type=line[14];
    ss_chain=line[16];

    strcpy(ss_list[num].pdb_id, pdb_id);
    ss_list[num].pdb_residue = atoi(pdb_residue_c);
    ss_list[num].ss_length = atoi(ss_length_c);
    ss_list[num].ss_type = ss_type;
    ss_list[num].ss_chain = ss_chain;
    num=num+1;
  }
  fclose(dssp_index);
  int is;
  for (is=0; is<num; is++) {
    if (ss_list[is].ss_type == search_type && ss_list[is].ss_length >= search_length) {

      FILE *dssp_file;
      char url[50];
      sprintf(url, "/data17/DSSP/%s.dssp",ss_list[is].pdb_id);
      if((dssp_file = fopen(url,"r")) == NULL) {
        //printf("No file\n");
        continue;
      }

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

      int line_position = 0, count=2*search_length-2;
      line_position=ss_list[is].pdb_residue - 1;
      fseek(dssp_file,i+line_position*137,SEEK_SET);

      char lines[150],phi_c[8],psi_c[8];
      int angle_i;
      double phi,psi;
      double dihedral_angle[1000]={0};
      if(ss_list[is].ss_length < search_length+2+search_angle_type){
        continue;
      }

      for (angle_i=0; angle_i<ss_list[is].ss_length; angle_i++) {
        fgets(lines, 140, dssp_file);
        strncpy(phi_c, lines + 103, 6);
        phi_c[6] = '\0';
        strncpy(psi_c, lines + 109, 6);
        psi_c[6] = '\0';
        phi = atof(phi_c);
        psi = atof(psi_c);
        dihedral_angle[2*angle_i] = phi;
        dihedral_angle[2*angle_i+1] = psi;
      }
      fclose(dssp_file);

      int angle_index = 0, line = 0;
      for(angle_index = 0; angle_index<(2*angle_i-3-2*search_angle_type); angle_index++){
        printf("%+7.2lf ",dihedral_angle[angle_index+2+search_angle_type]);
        if((angle_index+1)%count==0){
          printf("\n");
          line = line + 1;
        }
        if((2*angle_i-3-2*search_angle_type)-line*count<count){
          break;
        }
      }
    }
  }
  return 0;
}
