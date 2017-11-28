


int * extended_oIndex;
int * extended_IdStruct;
int * extended_IdHost;
int * extended_IdIGM;



void init_stfOutput (struct stfOutput * tmpstf);
void free_stfOutput (struct stfOutput * tmpstf);
void fill_SubIDS    (struct stfOutput * tmpstf);
void fill_ProgIDs (struct stfOutput * tmpstf, char * tffile);
int load_stf_extended_output (char * prefix, int filenum);
void read_properties_file (struct stfOutput * tmpstf);
void free_extended_arrays (void);
int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);
int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts);
