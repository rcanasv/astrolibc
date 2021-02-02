

/*
  int * files;
  int nfiles;
  nfiles = read_stf_filesofgroup (stfprefix, ID, &files);
*/



/*

  if ( (strcmp(format, "gadget") != 0) && (strcmp(format,"ramses") != 0) )
  {
    printf ("Incorrect specified format  %s, exiting...\n", format);
    exit (0);
  }


//     printf ("%d\n", snap1.strctProps[ID].NumPart);

//     if (!(p = malloc (snap1.strctProps[ID].NumPart * sizeof(struct pdata_s))))
//     {
//       printf ("Cannot allocate memory for particle information for structure %d\n", i);
//       printf ("Exiting\n");
//       exit (0);
//     }
//
//     int ninextended = 0;
//     int dummystruct = 0;
//     int dummyindex = 0;
//     int acum = 0;
//     for (i = 0; i < nfiles; i++)
//     {
//       printf ("file %d  %d\n", i, files[i]);
//       ninextended = load_stf_extended_output (stfprefix, files[i]);
//       if (ninextended)
//       {
//         read_ramses_snapshot (directory, galprefix, files[i]);
//       printf ("HERE\n");
//         for (j = 0; j < ninextended; j++)
//         {
//           dummystruct = extended_IdStruct[j];
//           dummyindex  = extended_oIndex[j];
//           if (dummystruct == ID)
//           {
//             p[acum].Pos[0] = ramses_pos[0][dummyindex];
//             p[acum].Pos[1] = ramses_pos[1][dummyindex];
//             p[acum].Pos[2] = ramses_pos[2][dummyindex];
//             p[acum].Vel[0] = ramses_vel[0][dummyindex];
//             p[acum].Vel[1] = ramses_vel[1][dummyindex];
//             p[acum].Vel[2] = ramses_vel[2][dummyindex];
//             p[acum].Mass   = ramses_mass  [dummyindex];
//             p[acum].Id     = ramses_id    [dummyindex];
//             p[acum].Age    = ramses_age   [dummyindex];
//             p[acum].Type = 4;
//             acum++;
//           }
//         }
//       printf ("HERE\n");
//         free_extended_arrays ();
//         free_ramses_arrays ();
//       }
//     }
//
//     write_snapshot(p, acum, header1, outprefix);


//     exit(0);


//     int    * prog_ids;
//     float  * prog_mrrts;
//     int      nprogs;
//     char tfname [NAME_LENGTH];
//     sprintf (tfname, "tftest");
//     nprogs = load_treefrog (tfname, 1000001, &prog_ids, &prog_mrrts);
//     exit(0);




//     struct stfOutput snap2;

//     init_stfOutput (&snap2);

//     sprintf (snap1.prefix, "782_gals");
//     sprintf (snap2.prefix, "772_gals");



//     char fname [NAME_LENGTH];
//     sprintf (fname, "MassSize_veloci_%03d.dat", atoi(stfprefix));
//     f = fopen(fname, "w");
//     for (i = 1; i <= snap1.nstruct; i++)
//       if (snap1.strctProps[i].Type > 7)
//       {
//         fprintf (f, "%e  ", snap1.strctProps[i].RHalfMass);
//         fprintf (f, "%e  ", snap1.strctProps[i].TotMass);
//         fprintf (f, "\n");
//       }
//     fclose (f);

//     printf ("Properties file read\n");


//     read_properties_file (&snap2);

//     fill_SubIDS (&snap1);

//     fill_SubIDS (&snap2);

//     fill_ProgIDs (&snap1, "tftest");
//
//     int    id1, id2, tmpid;
//     double m1, m2;
//     double ratio;
//     int    nmajor, nminor, nsmooth;
//     double maxratio;
//     double tmmajor, tmminor, tmsmooth;
//
//
//     FILE * fmergers = fopen ("mergers.dat", "w");
//
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (snap1.strctProps[i].Type != 7)
//       {
//         nmajor  = 0;
//         nminor  = 0;
//         nsmooth = 0;
//
//         m1       = 0.0;
//         ratio    = 0.0;
//         maxratio = 0.0;
//
//         tmmajor  = 0.0;
//         tmminor  = 0.0;
//         tmsmooth = 0.0;
//
//         if (snap1.strctProps[i].NumProg > 1)
//         {
//           for (j = 0; j < snap1.strctProps[i].NumProg; j++)
//           {
//             tmpid = snap1.strctProps[i].ProgIDs[j];
//             if (snap2.strctProps[tmpid].TotMass > m1)
//             {
//               m1  = snap2.strctProps[tmpid].TotMass;
//               id1 = tmpid;
//             }
//           }
//
//           for (j = 0; j < snap1.strctProps[i].NumProg; j++)
//           {
//             id2 = snap1.strctProps[i].ProgIDs[j];
//
//             if (id1 != id2)
//             {
//               m2 = snap2.strctProps[id2].TotMass;
//
//               ratio = m2/m1;
//
//               if (ratio > 0.2)
//               {
//                 nmajor++;
//                 tmmajor += m2;
//               }
//               else
//               {
//                 if (ratio > 0.001)
//                 {
//                   nminor++;
//                   tmminor += m2;
//                 }
//                 else
//                 {
//                   nsmooth++;
//                   tmsmooth += m2;
//                 }
//               }
//             }
//
//             if (ratio > maxratio)
//               maxratio = ratio;
//           }
//           fprintf (fmergers, "%e   %e   %e    %d   %e   %d   %e   %d  %e\n", snap1.strctProps[i].TotMass, m1, maxratio, \
//                    nmajor, tmmajor, nminor, tmminor, nsmooth, tmsmooth);
//         }
//       }
//     }
//     fclose (fmergers);

//     int gal1, gal2;
//     double * pos1, * pos2;
//     double rmax1, rmax2;
//     double mass1, mass2;
//     double rsize1, rsize2;
//     double dx;
//     double ratio;
//     double larger;
//
//
//     FILE * fpairs  = fopen ("pairs.dat","w");
//
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (snap1.strctProps[i].Type == 7)
//       {
//         for (j = 0; j < snap1.strctProps[i].NumSubs; j++)
//         {
//           for (k = j+1; k < snap1.strctProps[i].NumSubs; k++)
//           {
//             gal1 = snap1.strctProps[i].SubIDs[j];
//             gal2 = snap1.strctProps[i].SubIDs[k];
//
//             mass1 = snap1.strctProps[gal1].TotMass;
//             mass2 = snap1.strctProps[gal2].TotMass;
//
//             pos1 = snap1.strctProps[gal1].Pos;
//             pos2 = snap1.strctProps[gal2].Pos;
//
//             rmax1 = snap1.strctProps[gal1].Rvmax;
//             rmax2 = snap1.strctProps[gal2].Rvmax;
//
//             rsize1 = snap1.strctProps[gal1].Rsize;
//             rsize2 = snap1.strctProps[gal2].Rsize;
//
//             dx = (pos1[0] - pos2[0]) * (pos1[0] - pos2[0]) + \
//                  (pos1[1] - pos2[1]) * (pos1[1] - pos2[1]) + \
//                  (pos1[2] - pos2[2]) * (pos1[2] - pos2[2]);
//
//             if (mass1 > mass2)
//             {
//               larger = mass1;
//               ratio  = mass2 / mass1;
//             }
//             else
//             {
//               larger = mass2;
//               ratio  = mass1 / mass2;
//             }
//
//             if (dx < rmax1*rmax1 || dx < rmax2*rmax2)
//               fprintf (fpairs, "%e  %e  %e  %d\n", larger, sqrt(dx), ratio, 0);
//             else
//             {
//               if (dx < rsize1*rsize1 || dx < rsize2*rsize2)
//                 fprintf (fpairs, "%e  %e  %e  %d\n", larger, sqrt(dx), ratio, 1);
//             }
//           }
//         }
//       }
//     }
//
//     fclose (fpairs);



//     free_stfOutput (&snap2);




//     load_particles (format, nstruct, strctProps)
//
//
//
    //
    //  Load extended output to extract particles of structures
    //
//     int dummystruct;
//     int dummyindex;


//     P = malloc ((snap1.nstruct+1) * sizeof(struct pdata_s *));
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (!(P[i] = malloc (snap1.strctProps[i].NumPart * sizeof(struct pdata_s))))
//       {
//         printf ("Cannot allocate memory for particle information for structure %d\n", i);
//         printf ("Exiting\n");
//         exit (0);
//       }
//     }
//
//     int gadget = 0;
//     int ramses = 0;
//
//     if ( (strcmp(format, "gadget") == 0) )
//       gadget = 1;
//     if ( (strcmp(format, "ramses") == 0) )
//       ramses = 1;
//
//     int * acum = (int *) malloc ((snap1.nstruct+1) * sizeof(int));
//     for (i = 1; i <= snap1.nstruct; i++)
//       acum[i] = 0;

//     int ninextended;
//     if (ramses)
//     {
//       for (i = 0; i < NUMFILES; i++)
//       {
//         ninextended = load_stf_extended_output (stfprefix, i);
//         if (ninextended)
//         {
//           read_ramses_snapshot (directory, galprefix, i);
//           for (j = 0; j < ninextended; j++)
//           {
//             dummystruct = extended_IdStruct[j];
//             dummyindex  = extended_oIndex[j];
//
//             P[dummystruct][acum[dummystruct]].Pos[0] = ramses_pos[0][dummyindex];
//             P[dummystruct][acum[dummystruct]].Pos[1] = ramses_pos[1][dummyindex];
//             P[dummystruct][acum[dummystruct]].Pos[2] = ramses_pos[2][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[0] = ramses_vel[0][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[1] = ramses_vel[1][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[2] = ramses_vel[2][dummyindex];
//             P[dummystruct][acum[dummystruct]].Mass   = ramses_mass  [dummyindex];
//             P[dummystruct][acum[dummystruct]].Id     = ramses_id    [dummyindex];
//             P[dummystruct][acum[dummystruct]].Age    = ramses_age   [dummyindex];
//
//             P[dummystruct][acum[dummystruct]].Type = 4;
//
//             acum[dummystruct]++;
//           }
//
//           free_extended_arrays ();
//           free_ramses_arrays ();
//         }
//       }

//       double boxlen;
//       double time;
//       double aexp;
//       double H0;
//       double Omega_m;
//       double Omega_l;
//       double Omega_b;
//       char name [NAME_LENGTH];
//       char dummys [NAME_LENGTH];
//       double t;
//
//       sprintf (name, "%s/info_%s.txt", directory, galprefix);
//       FILE * ff = fopen(name,"r");
//       for (i = 0; i < 7; i++)
//         fgets(buffer, 100, ff);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &boxlen);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &time);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &aexp);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &H0);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_m);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_l);
//       fclose (ff);
//
//       double * axp_frw;
//       double * hexp_frw;
//       double * tau_frw;
//       double * t_frw;
//       int      n_frw = 1000;
//       double   time_tot = friedman(Omega_m, Omega_l, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);
//
//       // Find neighbouring conformal time
//       i = 1;
//       while (tau_frw[i] > time  && i < n_frw)
//         i = i+1;
//
//       // Interpolate time
//       double time_simu = t_frw[i]   * (time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
//                         t_frw[i-1] * (time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);
//
//       printf ("Time simu    %lf\n", (time_tot + time_simu) / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//       printf ("Hubble time  %lf\n", time_tot / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//
//       printf ("i               %d\n", i);
//       printf ("time            %e\n", time);
//       printf ("time_tot        %e\n", time_tot);
//       printf ("time_simu       %e\n", time_simu);
//       printf ("t_tot + t_simu  %e\n", time_tot + time_simu);
//
//       for (i = 1; i <= snap1.nstruct; i++)
//       {
//         for (j = 0; j < snap1.strctProps[i].NumPart; j++)
//         {
//           k = 1;
//           while (tau_frw[k] > P[i][j].Age  && k < n_frw)
//             k++;
//
//           t = t_frw[k]   * (P[i][j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
//               t_frw[k-1] * (P[i][j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);
//           P[i][j].Age = (time_simu - t) / (H0*1e5/3.08e24) / (365*24*3600.0);
//         }
//       }
//       free (axp_frw);
//       free (hexp_frw);
//       free (tau_frw);
//       free (t_frw);
//     }
//     char tmpbuffer[100];


//     if (gadget)
//     {
//       for (i = 0; i < NUMFILES; i++)
//       {
//         if (NUMFILES == 1)
//           sprintf (tmpbuffer, "%s", galprefix);
//         else
//           sprintf (tmpbuffer, "%s.%d", galprefix, i);
//
//         ninextended = 0;
//         ninextended = load_stf_extended_output (stfprefix, i);
//
//         if (ninextended)
//         {
//           read_gadget_snapshot (tmpbuffer);
//           for (j = 0; j < ninextended; j++)
//           {
//             dummystruct = extended_IdStruct[j];
//             dummyindex  = extended_oIndex[j];
//
//             P[dummystruct][acum[dummystruct]].Pos[0] = p[dummyindex].Pos[0];
//             P[dummystruct][acum[dummystruct]].Pos[1] = p[dummyindex].Pos[1];
//             P[dummystruct][acum[dummystruct]].Pos[2] = p[dummyindex].Pos[2];
//
//             P[dummystruct][acum[dummystruct]].Vel[0] = p[dummyindex].Vel[0];
//             P[dummystruct][acum[dummystruct]].Vel[1] = p[dummyindex].Vel[1];
//             P[dummystruct][acum[dummystruct]].Vel[2] = p[dummyindex].Vel[2];
//
//             P[dummystruct][acum[dummystruct]].Mass   = p[dummyindex].Mass;
//             P[dummystruct][acum[dummystruct]].Id     = p[dummyindex].Id;
//
//             P[dummystruct][acum[dummystruct]].Type = 2;
//
//             acum[dummystruct]++;
//           }
//           free (p);
//           free (ID);
//           free_extended_arrays ();
//         }
//       }
//     }


//     char bfr [NAME_LENGTH];
//     sprintf (bfr, "sSFR_%s.dat", galprefix);
//     FILE * fgal = fopen (bfr, "w");
//
//     double Mnew;
//     int nage;
//     double Agesum;
//     double mtot;
//     double Dt = 100e6;
//     for (i = 1; i <= snap1.nstruct; i++)
//       if (snap1.strctProps[i].Type != 7)
//       {
//         Mnew = 0;
//         mtot = 0;
//         nage = 0;
//         Agesum = 0;
//         for (j = 0; j < snap1.strctProps[i].NumPart; j++)
//         {
//           if (P[i][j].Age < Dt && P[i][j].Age > 0.0)
// //           if (P[i][j].Age < Dt)
//           {
//             Mnew += P[i][j].Mass;
//             nage++;
//             Agesum += P[i][j].Age;
//           }
//           mtot += P[i][j].Mass;
//         }
//         if (nage)
//           fprintf (fgal, "%d   %e   %e   %e   %e\n", i, snap1.strctProps[i].TotMass, mtot, Mnew/Dt, Agesum/(double)nage);
//       }
//
//     fclose (fgal);



    for (i = 1; i <= snap1.nstruct; i++)
      if (snap1.strctProps[i].Type != 7)
      {
        Mnew = 0;
        for (j = 0; j < snap1.strctProps[i].NumPart; j++)
          Mnew += P[i][j].Mass;
        printf ("%d  %e   %e\n", i, snap1.strctProps[i].TotMass, Mnew);
      }   */

//     for (i = 1; i <= snap1.nstruct; i++)
//       free (P[i]);
//     free (P);

//    free_stfOutput (&snap1);



//     FILE * ff;
//     int PID, STID, HSTID, IGMID;
//
//     printf("Distributing particles\n");
//     for (i = 0; i < NUMFILES; i++)
//     {
//       sprintf (buffer, "%s.%d", galprefix, i);
//       read_gadget_snapshot (buffer);
//       sprintf (buffer, "%s.extended.%d", stfprefix, i);
//       ff = fopen (buffer, "r");
//       for (j = 0; j < NumPart; j++)
//       {
// 	fgets  (longbuffer, NAME_LENGTH, ff);
// 	sscanf (longbuffer, "%d %d %d %d", &PID, &STID, &HSTID, &IGMID);
// 	if (STID > 0)
// 	{
// 	  P[STID-1][acum[STID-1]].Id = p[j].Id;
// 	  P[STID-1][acum[STID-1]].Pos[0] = p[j].Pos[0];
// 	  P[STID-1][acum[STID-1]].Pos[1] = p[j].Pos[1];
// 	  P[STID-1][acum[STID-1]].Pos[2] = p[j].Pos[2];
// 	  P[STID-1][acum[STID-1]].Vel[0] = p[j].Vel[0];
// 	  P[STID-1][acum[STID-1]].Vel[1] = p[j].Vel[1];
// 	  P[STID-1][acum[STID-1]].Vel[2] = p[j].Vel[2];
// 	  P[STID-1][acum[STID-1]].Mass = p[j].Mass;
// 	  P[STID-1][acum[STID-1]].Type = p[j].Type;
// 	  acum[STID-1]++;
// 	}
//       }
//       fclose(ff);
//       free (p);
//     }
//     printf("Particles Distributed\n");

//     //
//     // Write Gadget-formatted files
//     //
//     for (i = 0; i < nstruct; i++)
//     {
//       printf ("struct %d nparts %d acum %d \n", i, nparts[i], acum[i]);
//       if (outType == 0)
//         sprintf (output_fname, "%s_%03d.gdt_%03d", outprefix, hostid[i], idins[i]-1);
//       if (outType == 1)
//         sprintf (output_fname, "%s_%03d", outprefix, i+1);
//       write_snapshot(P[i], acum[i], header1, output_fname);
//     }
