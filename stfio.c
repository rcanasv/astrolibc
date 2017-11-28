
//---------- Read Propertes File ----------//

void read_properties_file (struct stfOutput * tmpstf)
{
  int    i, j, k;

  int    dummyi;
  long   dummyl;
  float  dummyf;
  double dummyd;

  FILE * f;
  char   propts_fname [NAME_LENGTH];
  int    mystructs;

  //
  // Open properties file to read total number of structures
  // and number of processors if stf was run with MPI
  //
  sprintf (propts_fname, "%s.properties", tmpstf->prefix);
  if ((f = fopen (propts_fname, "r")) == NULL)
  {
    sprintf (propts_fname, "%s.properties.0", tmpstf->prefix);
    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s\n", propts_fname);
      exit (0);
    }
  }
  fgets  (longbuffer, NAME_LENGTH, f);  sscanf (longbuffer, "%d  %d", &dummyi, &(tmpstf->nprocs));
  fgets  (longbuffer, NAME_LENGTH, f);  sscanf (longbuffer, "%d  %d", &dummyi, &(tmpstf->nstruct));
  fclose (f);

  //
  // Allocate memory for structure properties
  //
  if ( (tmpstf->strctProps = (struct objProps *) malloc ((tmpstf->nstruct+1) * sizeof(struct objProps))) == NULL)
  {
    printf ("Couldn't allocate memory.\n");
    printf ("Exiting...\n");
    exit (0);
  }
  else
    tmpstf->iprops = 1;

  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    tmpstf->strctProps[i].SubIDs    = NULL;
    tmpstf->strctProps[i].ProgIDs   = NULL;
    tmpstf->strctProps[i].ProgMrrts = NULL;
  }

  //
  // Read file(s) and store desired structure properties
  //
  int offst = 1;

  for (i = 0; i < tmpstf->nprocs; i++)
  {
    if (tmpstf->nprocs == 1)
      sprintf (propts_fname, "%s.properties", tmpstf->prefix);
    else
      sprintf (propts_fname, "%s.properties.%d", tmpstf->prefix, i);

//     printf ("Openning file  %s  \n", propts_fname);

    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s", propts_fname);
      exit (0);
    }

    fgets  (longbuffer, NAME_LENGTH, f);
    fgets  (longbuffer, NAME_LENGTH, f);
    sscanf (longbuffer, "%d  %d", &mystructs, &dummyi);
    fgets  (longbuffer, 3000, f);

    for (j = 0; j < mystructs; j++)
    {
      fgets (longbuffer, 3000, f);
      sscanf (longbuffer, "%d  %d  %d  %d  %d  %d  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf                \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf            \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf            \
              %lf  %lf",                                                                                \
              &(tmpstf->strctProps[j+offst].ID), &dummyi, &(tmpstf->strctProps[j+offst].DirectHostID),  \
              &(tmpstf->strctProps[j+offst].HostID), &(tmpstf->strctProps[j+offst].NumSubs),            \
              &(tmpstf->strctProps[j+offst].Type), &(tmpstf->strctProps[j+offst].NumPart), &dummyd,     \
              &(tmpstf->strctProps[j+offst].Pos[0]), &(tmpstf->strctProps[j+offst].Pos[1]),             \
              &(tmpstf->strctProps[j+offst].Pos[2]), &dummyd, &dummyd, &dummyd,                         \
              &(tmpstf->strctProps[j+offst].Vel[0]), &(tmpstf->strctProps[j+offst].Vel[1]),             \
              &(tmpstf->strctProps[j+offst].Vel[2]), &dummyd, &dummyd, &dummyd,                         \
              &(tmpstf->strctProps[j+offst].TotMass), &dummyd, &dummyd, &dummyd, &dummyd,               \
              &(tmpstf->strctProps[j+offst].Efrac), &dummyd, &(tmpstf->strctProps[j+offst].Rsize),      \
              &dummyd, &dummyd, &dummyd, &(tmpstf->strctProps[j+offst].RHalfMass),                      \
              &(tmpstf->strctProps[j+offst].Rvmax), &(tmpstf->strctProps[j+offst].Vmax),                \
              &(tmpstf->strctProps[j+offst].Vdisp), &dummyd, &dummyd, &dummyd, &dummyd, &dummyd,        \
              &dummyd, &dummyd, &dummyd, &dummyd, &(tmpstf->strctProps[j+offst].Lambda),                \
              &(tmpstf->strctProps[j+offst].L[0]), &(tmpstf->strctProps[j+offst].L[1]),                 \
              &(tmpstf->strctProps[j+offst].L[2])                                                       \
            );
    }
    offst += mystructs;
    fclose(f);
  }
}




int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct)
{
  int i, j, k;
  FILE * f;

  char buffer [NAME_LENGTH];

  sprintf (buffer, "%s.filesofgroup", prefix);

  int tmpid;
  int nfiles;

  f = fopen (buffer, "r");
  do
  {
    fgets (buffer, NAME_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nfiles);
    fgets (buffer, NAME_LENGTH, f);
    get_n_num_from_string (buffer, nfiles, files_of_strct);
  }
  while (tmpid != strct_id);

  fclose (f);

  return nfiles;
}




int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts)
{
  int i, j, k;
  FILE * f;

  int tmpid;
  int nprogs;

  char buffer [LONG_LENGTH];

  f = fopen (tffile, "r");

  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);

  do
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);

    if (tmpid != strct_id)
      for (i = 0; i < nprogs; i++)
        fgets (buffer, LONG_LENGTH, f);
  }
  while (tmpid != strct_id);


  *prog_ids   = (int *)   malloc (nprogs * sizeof(int));
  *prog_mrrts = (float *) malloc (nprogs * sizeof(float));

  for (i = 0; i < nprogs; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %f", &prog_ids[0][i], &prog_mrrts[0][i]);
  }
  fclose (f);

  return nprogs;
}




int load_stf_extended_output (char * prefix, int filenum)
{
    FILE * f;
    char buffer[NAME_LENGTH];
    sprintf(buffer, "%s.extended.%d", prefix, filenum);
    int nparts = 0;
    int i;

    if ((f = fopen(buffer, "r")) == NULL)
      return 0;

    while (fgets(buffer, NAME_LENGTH, f) != NULL)
      nparts++;
    rewind(f);

//     printf("nparts %d\n", nparts);

    extended_oIndex   = (int *) malloc (nparts * sizeof(int));
    extended_IdStruct = (int *) malloc (nparts * sizeof(int));
    extended_IdHost   = (int *) malloc (nparts * sizeof(int));
    extended_IdIGM    = (int *) malloc (nparts * sizeof(int));

    for (i = 0; i < nparts; i++)
    {
      fgets(buffer, NAME_LENGTH, f);
      sscanf(buffer, "%d  %d  %d  %d  ", &extended_oIndex[i], &extended_IdStruct[i], &extended_IdHost[i], &extended_IdIGM[i]);
    }

    fclose (f);

    return nparts;
}

void free_extended_arrays (void)
{
  free (extended_oIndex);
  free (extended_IdStruct);
  free (extended_IdHost);
  free (extended_IdIGM);
}




void init_stfOutput (struct stfOutput * tmpstf)
{
  tmpstf->nprocs     = 0;
  tmpstf->nstruct    = 0;
  tmpstf->iprops     = 0;
  tmpstf->iparts     = 0;
  tmpstf->strctProps = NULL;
  tmpstf->strctParts = NULL;
}

void free_stfOutput (struct stfOutput * tmpstf)
{
  int i;

  if (tmpstf->iprops)
  {
    for (i = 1; i <= tmpstf->nstruct; i++)
    {
      if (tmpstf->strctProps[i].SubIDs != NULL)
        free (tmpstf->strctProps[i].SubIDs);

      if (tmpstf->strctProps[i].ProgIDs != NULL)
        free (tmpstf->strctProps[i].ProgIDs);

      if (tmpstf->strctProps[i].ProgMrrts != NULL)
        free (tmpstf->strctProps[i].ProgMrrts);
    }
  }

  if (tmpstf->iparts)
  {
    for (i = 1; i <= tmpstf->nstruct; i++)
      free (tmpstf->strctParts[i]);
    free (tmpstf->strctParts);
  }
}



void fill_SubIDS (struct stfOutput * tmpstf)
{
  int i, j;
  int bob;
  int tmp;


  for (i = 1; i <= tmpstf->nstruct; i++)
    if (tmpstf->strctProps[i].HostID == -1)
    {
      bob = tmpstf->strctProps[i].NumSubs;
      tmpstf->strctProps[i].SubIDs = (int *) malloc (bob * sizeof(int));
      tmpstf->strctProps[i].dummy = 0;
      for (j = 0; j < bob; j++)
        tmpstf->strctProps[i].SubIDs[j] = 0;
    }

  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    if (tmpstf->strctProps[i].HostID != -1)
    {
      bob = i;
      while (tmpstf->strctProps[bob].HostID != -1)
      {
        bob = tmpstf->strctProps[bob].DirectHostID;
      }
      tmp = tmpstf->strctProps[bob].dummy++;
      tmpstf->strctProps[bob].SubIDs[tmp] = i;
    }
  }
}


void fill_ProgIDs (struct stfOutput * tmpstf, char * tffile)
{
  int i, j;
  int bob;
  int tmpid ;
  int nprogs;

  FILE * f = fopen (tffile, "r");

  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);

  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);

    tmpstf->strctProps[i].NumProg = nprogs;
    if (nprogs > 0)
    {
      tmpstf->strctProps[i].ProgIDs   = (int *)    malloc (nprogs * sizeof(int));
      tmpstf->strctProps[i].ProgMrrts = (double *) malloc (nprogs * sizeof(double));
      for (j = 0; j < nprogs; j++)
      {
        fgets  (buffer, LONG_LENGTH, f);
        sscanf (buffer, "%d  %lf", &(tmpstf->strctProps[i].ProgIDs[j]), &(tmpstf->strctProps[i].ProgMrrts[j]));
      }
    }
  }

  fclose (f);
}
