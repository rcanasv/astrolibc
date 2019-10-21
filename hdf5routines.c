/*
 *  \file   hdf5.c
 *  \brief  This file contains functions for HDF5.
 *
 */

#include "base.h"
#include "hdf5routines.h"



void hdf5_get_attribute (hid_t obj_id, char * attr_name, void * buffer, size_t vsize)
{
 int i, j, k;

 int           dummyi;
 long          dummyl;
 float         dummyf;
 double        dummyd;
 char          dummys[NAME_LENGTH];

 int         * dummyiarr = NULL;
 long        * dummylarr = NULL;
 float       * dummyfarr = NULL;
 double      * dummydarr = NULL;

 int           sizearr = 1;

 size_t        size;

 hid_t         attr_id;
 hid_t         type_id;
 hid_t         space_id;

 hsize_t       ndims;
 hsize_t     * dims    = NULL;
 hsize_t     * maxdims = NULL;

 H5T_class_t   type_class;  /* H5T_INTEGER, H5T_FLOAT, H5T_STRING */
 H5T_order_t   type_order;  /* H5T_ORDER_LE, H5T_ORDER_BE */


 if (H5Aexists (obj_id, attr_name) > 0)
 {
   attr_id     = H5Aopen      (obj_id, attr_name, H5P_DEFAULT);
   type_id     = H5Aget_type  (attr_id);
   type_class  = H5Tget_class (type_id);
   size        = H5Tget_size  (type_id);

   space_id    = H5Aget_space (attr_id);
   ndims       = H5Sget_simple_extent_ndims (space_id);
   dims        = (hsize_t *) malloc (ndims * sizeof(hsize_t));
   maxdims     = (hsize_t *) malloc (ndims * sizeof(hsize_t));
   ndims       = H5Sget_simple_extent_dims (space_id, dims, maxdims);

   switch (type_class)
   {
     case H5T_INTEGER:
       if (size == sizeof(int))
       {
         if (ndims == 0)
         {
           H5Aread (attr_id, type_id, &dummyi);
           if (vsize == sizeof(int))  *(int  *) buffer = (int)  dummyi;
           if (vsize == sizeof(long)) *(long *) buffer = (long) dummyi;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummyiarr = (int *) malloc (sizearr * sizeof(int));
           H5Aread (attr_id, type_id, dummyiarr);

           if (vsize == sizeof(int))
             for (i = 0; i < sizearr; i++)
             {
               *(int*) buffer = (int) dummyiarr[i];
               buffer = (char *) buffer + sizeof(int);
             }

           if (vsize == sizeof(long))
             for (i = 0; i < sizearr; i++)
             {
               *(long*) buffer = (long) dummyiarr[i];
               buffer = (char *) buffer + sizeof(long);
             }
           free (dummyiarr);
         }
       }
       else
       {
         if (ndims == 0)
         {
           H5Aread (attr_id, type_id, &dummyl);
           if (vsize == sizeof(int))  *(int * ) buffer = (int)  dummyl;
           if (vsize == sizeof(long)) *(long *) buffer = (long) dummyl;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummylarr = (long *) malloc (sizearr * sizeof(long));
           H5Aread (attr_id, type_id, dummylarr);

           if (vsize == sizeof(int))
             for (i = 0; i < sizearr; i++)
             {
               *(int*) buffer = (int) dummylarr[i];
               buffer = (char *) buffer + sizeof(int);
             }

           if (vsize == sizeof(long))
             for (i = 0; i < sizearr; i++)
             {
               *(long*) buffer = (long) dummylarr[i];
               buffer = (char *) buffer + sizeof(long);
             }
           free (dummylarr);
         }
       }
       break;


     case H5T_FLOAT:
       if (size == sizeof(float))
       {
         if (ndims == 0)
         {
           H5Aread (attr_id, type_id, &dummyf);
           if (vsize == sizeof(float))   *(float  *) buffer = (float)  dummyf;
           if (vsize == sizeof(double))  *(double *) buffer = (double) dummyf;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummyfarr = (float *) malloc (sizearr * sizeof(float));
           H5Aread (attr_id, type_id, dummyfarr);

           if (vsize == sizeof(float))
             for (i = 0; i < sizearr; i++)
             {
               *(float*) buffer = (float) dummyfarr[i];
               buffer = (char *) buffer + sizeof(float);
             }

           if (vsize == sizeof(double))
             for (i = 0; i < sizearr; i++)
             {
               *(double*) buffer = (double) dummyfarr[i];
               buffer = (char *) buffer + sizeof(double);
             }
           free (dummyfarr);
         }
       }
       else
       {
         if (ndims == 0)
         {
           H5Aread (attr_id, type_id, &dummyd);
           if (vsize == sizeof(float))  *(float  *) buffer = (float)   dummyd;
           if (vsize == sizeof(double)) *(double *) buffer = (double)  dummyd;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummydarr = (double *) malloc (sizearr * sizeof(double));
           H5Aread (attr_id, type_id, dummydarr);

           if (vsize == sizeof(float))
           {
             for (i = 0; i < sizearr; i++)
             {
               *(float*) buffer = (float) dummydarr[i];
               buffer = (char *) buffer + sizeof(float);
             }
           }

           if (vsize == sizeof(double))
           {
             for (i = 0; i < sizearr; i++)
             {
               *(double*) buffer = (double) dummydarr[i];
               buffer = (char *) buffer + sizeof(double);
             }
           }
           free (dummydarr);
         }
       }
       break;

     case H5T_STRING:
       H5Aread (attr_id, type_id, buffer);
       break;

     default:
       printf ("No class %d\n", (int) type_class);

   }

   H5Tclose (type_id);
   H5Sclose (space_id);
   H5Aclose (attr_id);
 }
 else
 {
   printf ("Cannot find attribute named  %s\n", attr_name);
   printf ("Exiting\n");
   exit (0);
 }

 if (dims    != NULL) free (dims);
 if (maxdims != NULL) free (maxdims);
}



void hdf5_get_data (hid_t obj_id, char * dset_name, void * buffer, size_t vsize)
{
 int i, j, k;

 int           dummyi;
 long          dummyl;
 float         dummyf;
 double        dummyd;
 char          dummys[NAME_LENGTH];

 int         * dummyiarr = NULL;
 long        * dummylarr = NULL;
 float       * dummyfarr = NULL;
 double      * dummydarr = NULL;

 int           sizearr = 1;

 size_t        size;

 hid_t         attr_id;
 hid_t         dset_id;
 hid_t         type_id;
 hid_t         space_id;

 hsize_t       ndims;
 hsize_t     * dims    = NULL;
 hsize_t     * maxdims = NULL;

 H5T_class_t   type_class;  /* H5T_INTEGER, H5T_FLOAT, H5T_STRING */
 H5T_order_t   type_order;  /* H5T_ORDER_LE, H5T_ORDER_BE */


 if (H5Oexists_by_name (obj_id, dset_name, H5P_DEFAULT) > 0)
 {
   dset_id     = H5Dopen      (obj_id, dset_name, H5P_DEFAULT);

   type_id     = H5Dget_type  (dset_id);
   type_class  = H5Tget_class (type_id);
   size        = H5Tget_size  (type_id);

   space_id    = H5Dget_space (dset_id);
   ndims       = H5Sget_simple_extent_ndims (space_id);
   dims        = (hsize_t *) malloc (ndims * sizeof(hsize_t));
   maxdims     = (hsize_t *) malloc (ndims * sizeof(hsize_t));
   ndims       = H5Sget_simple_extent_dims (space_id, dims, maxdims);

   switch (type_class)
   {
     case H5T_INTEGER:
       if (size == sizeof(int))
       {
         if (ndims == 0)
         {
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummyi);
           if (vsize == sizeof(int))  *(int  *) buffer = (int)  dummyi;
           if (vsize == sizeof(long)) *(long *) buffer = (long) dummyi;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummyiarr = (int *) malloc (sizearr * sizeof(int));
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummyiarr);

           if (vsize == sizeof(int))
             for (i = 0; i < sizearr; i++)
             {
               *(int*) buffer = (int) dummyiarr[i];
               buffer = (char *) buffer + sizeof(int);
             }

           if (vsize == sizeof(long))
             for (i = 0; i < sizearr; i++)
             {
               *(long*) buffer = (long) dummyiarr[i];
               buffer = (char *) buffer + sizeof(long);
             }
           free (dummyiarr);
         }
       }
       else
       {
         if (ndims == 0)
         {
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummyl);
           if (vsize == sizeof(int))  *(int * ) buffer = (int)  dummyl;
           if (vsize == sizeof(long)) *(long *) buffer = (long) dummyl;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummylarr = (long *) malloc (sizearr * sizeof(long));
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummylarr);

           if (vsize == sizeof(int))
             for (i = 0; i < sizearr; i++)
             {
               *(int*) buffer = (int) dummylarr[i];
               buffer = (char *) buffer + sizeof(int);
             }

           if (vsize == sizeof(long))
             for (i = 0; i < sizearr; i++)
             {
               *(long*) buffer = (long) dummylarr[i];
               buffer = (char *) buffer + sizeof(long);
             }
           free (dummylarr);
         }
       }
       break;


     case H5T_FLOAT:
       if (size == sizeof(float))
       {
         if (ndims == 0)
         {
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummyf);
           if (vsize == sizeof(float))   *(float  *) buffer = (float)  dummyf;
           if (vsize == sizeof(double))  *(double *) buffer = (double) dummyf;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummyfarr = (float *) malloc (sizearr * sizeof(float));
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummyfarr);

           if (vsize == sizeof(float))
             for (i = 0; i < sizearr; i++)
             {
               *(float*) buffer = (float) dummyfarr[i];
               buffer = (char *) buffer + sizeof(float);
             }

           if (vsize == sizeof(double))
             for (i = 0; i < sizearr; i++)
             {
               *(double*) buffer = (double) dummyfarr[i];
               buffer = (char *) buffer + sizeof(double);
             }
           free (dummyfarr);
         }
       }
       else
       {
         if (ndims == 0)
         {
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dummyd);
           if (vsize == sizeof(float))  *(float  *) buffer = (float)   dummyd;
           if (vsize == sizeof(double)) *(double *) buffer = (double)  dummyd;
         }
         else
         {
           for (i = 0; i < ndims; i++)  sizearr *= dims[i];
           dummydarr = (double *) malloc (sizearr * sizeof(double));
           H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dummydarr);

           if (vsize == sizeof(float))
             for (i = 0; i < sizearr; i++)
             {
               *(float*) buffer = (float) dummydarr[i];
               buffer = (char *) buffer + sizeof(float);
             }

           if (vsize == sizeof(double))
             for (i = 0; i < sizearr; i++)
             {
               *(double*) buffer = (double) dummydarr[i];
               buffer = (char *) buffer + sizeof(double);
             }
           free (dummydarr);
         }
       }
       break;

     case H5T_STRING:
       H5Dread (dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
       break;

     default:
       printf ("No class %d\n", (int) type_class);
   }

   H5Tclose (type_id);
   H5Sclose (space_id);
   H5Dclose (dset_id);
 }
 else
 {
   printf ("Cannot find object named  %s\n", dset_name);
   printf ("Exiting\n");
   exit (0);
 }

 if (dims    != NULL) free (dims);
 if (maxdims != NULL) free (maxdims);
}



 herr_t hdf5_finfo (hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
 {
     hid_t group_id;
     //Open the group using its name via the c interface
     group_id = H5Gopen2(loc_id, name, H5P_DEFAULT);
     //Display group name.
     printf ("Group Name :  %s\n", name);
     //H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, file_data_info, NULL);
     //close the group via the c interface
     H5Gclose(group_id);
     return 0;
 }
