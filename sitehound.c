/*****************************************************************************
 *****************************************************************************
 ** sitehound.c by Dario Ghersi                                             **
 ** Version: 20101214                                                       **
 **                                                                         **
 ** Usage: please see manual for usage info                                 **
 **                                                                         **
 ** The following code contains modified parts from the Open Source         **
 ** Clustering Software in the functions create_distance_matrix,            **
 ** sl_cluster and average_cluster                                          **
 **                                                                         **
 ** Revisions:                                                              **
 ** 20101214  changed the atom type in the output (clusters.pdb)            **
 ** 20101201  added the option for using compressed EasyMIFs maps           **
 **           added the option to convert compressed maps back to dx files  **
 *****************************************************************************
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "prototypes.h"

/*---------------------------------------------------------------------------*\
 *  GLOBAL VARIABLES                                                         *
\*---------------------------------------------------------------------------*/

unsigned int npoints[3]; /* number of points in the grid in x, y, z */
double resolution = 1.0; /* the resolution of the grid */
double energy_cutoff = -0.4; /* the energy cutoff to filter out points */
double spatial_cutoff = 7.8; /* the spatial cutoff needed to cut
                                the partition */
unsigned int cl_num_contacts = 10; /* the number of clusters whose contact 
                                     residues will be computed */
double lower_corner[3]; /* coordinates of the lower corner of the box */
char filename[MAX_STR_LEN] = "-1"; /* the name of the map */
char maptype[MAX_STR_LEN] = ""; /* the type of the map (easymifs, autogrid) */
char linkage[MAX_STR_LEN] = "average" ; /* clustering parameter 
                                           ('average' or 'single') */
unsigned int cluster_atom_num = 0; /* to append clusters in PDB */
char clusters_chain_id; /* the chain identifier for the clusters in PDB */
char used_chains[TOTAL_CHAIN_IDS] = {'\0'}; /* array for storing the chains in the PDB */
unsigned int decompress = 0; /* flag to decompress a .cmp map back to .dx */

long seed = 1; /* the initial seed for the random jiggling of the points */

/*---------------------------------------------------------------------------*\
 *  FUNCTIONS DEFINITIONS                                                    * 
\*---------------------------------------------------------------------------*/

char *basename(char *s)
{
  /* return */
  char *b;
  char p[1000];
  char *temp;

  b = malloc(sizeof(char) * MAX_STR_LEN);
  temp = strtok(s, ".");
  strncpy(b, temp, MAX_STR_LEN);
  p[0] = '\0';
  while(temp = strtok(NULL, ".")) {
    strcat(b, p);
    strcpy(p, ".");
    strcat(p, temp);
  }
  return b;
}

/* -------------------------------------------------------------------------- */

double ran(long *seed)
{
  long k;
  double ans;

  *seed ^= MASK;
  k = (*seed) / IQ;
  *seed = IA * (*seed -k * IQ) - IR * k;
  if (*seed < 0)
    *seed += IM;

  ans = AM * (*seed);
  *seed ^= MASK;
  return ans;
}

/* -------------------------------------------------------------------------- */

double euclidean_distance(const Data *p, const Atom *a)
{
  /* returned the squared euclidean distance between a point and an atom */

  double ed = 0.0; /* the squared euclidean distance */

  ed = (p->x - a->x) * (p->x - a->x) + (p->y - a->y) * (p->y - a->y) +
       (p->z - a->z) * (p->z - a->z);
  return(ed);
}

/* -------------------------------------------------------------------------- */

void free_memory_restable(Restable *residue_list)
{
  Restable *next_residue;
  while (residue_list != NULL) {
    next_residue = residue_list->next;
    free(residue_list);
    residue_list = next_residue;
  }
  free(residue_list);
}

/* -------------------------------------------------------------------------- */

void set_parameters(int argc, char **argv)
{
  /* get the command line parameters and set them in the global variables */
  char *program_name;
  int return_value;

  program_name = argv[0];
  while ((argc > 1) && (argv[1][0] == '-')) {
    switch (argv[1][1]) {
      case 't': /* type of map */
        strtok(argv[1], "=");
        strncpy(maptype, strtok(NULL, " "), MAX_STR_LEN);
        break;
      case 'f': /* filename */
        strtok(argv[1], "=");
        strncpy(filename, strtok(NULL, " "), MAX_STR_LEN);
        break;
      case 'l': /* linkage option */
        if (argv[1][3] == 's') {
          strcpy(linkage, "single");
	}
        else if (argv[1][3] == 'a') {
          strcpy(linkage, "average");
        }
        else {
          fprintf(stderr, "Linkage option not recognized\nAborting...\n");
          exit(1);
        }
        break;

      case 'r': /* resolution */
        strtok(argv[1], "=");
	return_value = sscanf(strtok(NULL, " "), "%lf", &resolution);
        if (return_value == 0) {
          fprintf(stderr, "Resolution option not recognized\nAborting...\n");
          exit(1);
        }
        break;

      case 'e': /* energy cutoff */
        strtok(argv[1], "=");
	return_value = sscanf(strtok(NULL, " "), "%lf", &energy_cutoff);
        if (return_value == 0) {
          fprintf(stderr, "Energy option not recognized\nAborting...\n");
          exit(1);
        }
        if (energy_cutoff > 0.0) {
          fprintf(stderr, "Energy should be a negative number\nAborting...\n");
          exit(1);
        } 
        break;

      case 's': /* spatial cutoff */
        strtok(argv[1], "=");
	return_value = sscanf(strtok(NULL, " "), "%lf", &spatial_cutoff);
        if (return_value == 0) {
          fprintf(stderr, "Spatial option not recognized\nAborting...\n");
          exit(1);
        }
        break;

      case 'p': /* number of clusters for contact residues */
        strtok(argv[1], "=");
        return_value = sscanf(strtok(NULL, " "), "%d", &cl_num_contacts);
        if (return_value == 0) {
	  fprintf(stderr, "Contact residues option not recognized\nAborting...\n");
	  exit(1);
	}
        break;

      case 'z': /* decompress a .cmp map back to .dx file and exit */
        decompress = 1;
	break;
        
      default:
        fprintf(stderr, "Option not recognized\nAborting...\n");
        exit(1);
    }
    ++argv;
    --argc;
  }
}

/*---------------------------------------------------------------------------*/

void filter_autogrid_map(const char *filename, const double energy_cutoff)
{
  /* read in an autogrid map and writes a new file with only the filtered
     points (whose energy is below the user-defined threshold */

  FILE *infile, *outfile;
  char *cp, outfilename[MAX_STR_LEN];
  char line[MAX_STR_LEN];
  char *temp = NULL;
  bool end_of_header = FALSE;
  double center[3]; /* coordinates of the center */
  double coords[3]; /* coordinates of a point */
  double energy;
  unsigned int cx, cy, cz; /* counter variables for computing the coordinates */
  unsigned int index_point = 1; /* keep track of the points */

  /* open/check the file to be filtered */
  if (!(infile = fopen(filename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", filename);
    exit(1);
  }

  /* read the header of the autogrid file */
  while ( (!end_of_header) && fgets(line, MAX_STR_LEN, infile)) {
    if (strstr(line, "SPACING")) { /* grab the resolution of the grid */
      strtok(line, " ");
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &resolution);
    }

    if (strstr(line, "NELEMENTS")) { /* grab the number of points in x, y, z */
      strtok(line, " ");
      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[X]);
      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[Y]);
      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[Z]);
    }
    
    if (strstr(line, "CENTER")) { /* grab the coordinates of the center */
      strtok(line, " ");
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &center[X]);
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &center[Y]);
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &center[Z]);
      end_of_header = TRUE;
    }
  }
  if (! end_of_header) {
    fprintf(stderr, "Problem with autogrid file\nAborting...\n");
    exit(1);
  }

  /* compute the coordinates of the lower corner of the box */
  lower_corner[X] = center[X] - resolution * npoints[X] / 2;
  lower_corner[Y] = center[Y] - resolution * npoints[Y] / 2;
  lower_corner[Z] = center[Z] - resolution * npoints[Z] / 2;

  /* open the output file */
  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, ".tmp");
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s for writing\nAborting...\n", outfilename);
  }

  /* read the data from the autogrid file and write the filtered points to 
     the output file */
  for (cz = 0; cz < npoints[Z] + 1; cz++) {
    coords[Z] = lower_corner[Z] + resolution * cz;
    for (cy = 0; cy < npoints[Y] + 1; cy++) {
      coords[Y] = lower_corner[Y] + resolution * cy;
      for (cx = 0; cx < npoints[X] + 1; cx++) {
        coords[X] = lower_corner[X] + resolution * cx;
        if (fgets(line, MAX_STR_LEN, infile) == NULL) {
          fprintf(stderr, "Problems with autogrid file\n...Aborting\n");
          exit(1);
        }
        sscanf(line, "%lf", &energy); /* read the energy value */
        if (energy <= energy_cutoff) { /* if energy < cutoff write point */
	  fprintf(outfile, "%d %.3f %.3f %.3f %.3f\n", index_point,
                  coords[X], coords[Y], coords[Z], energy);
        }
        index_point++;
      }
    }
  }

  /* set npoints to the actual number of points  */
  npoints[X]++; npoints[Y]++; npoints[Z]++;

  /* close input and output files */
  fclose(infile);
  fclose(outfile);
}

/*---------------------------------------------------------------------------*/

double *expand(unsigned int *bytes, unsigned int num_lines, unsigned int num_values)
{
  /* expansion part of the LZW compression algorithm                          *
     Consult http://warp.povusers.org/EfficientLZW/index.html for an          *
     excellent discussion of the efficiency issues                            */

  FILE *temp_file;
  double *energy_values;
  unsigned int *bytes_ptr, index, old;
  unsigned int i, current_num_value = 0, string_length;
  const char dictionary[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.', '-', ' '};
  const unsigned int length_dictionary = 13;
  const unsigned int size_dict = 1 << BITS;
  const unsigned int last_position = 57; /* the last occupied position in the
                                          hash table */
  unsigned int dict_pos = last_position;
  char *dict[size_dict];
  char B;

  /* memory allocation */
  bytes_ptr = bytes;
  energy_values = malloc(sizeof(double) * num_values);

  /* initialize the hash table (with characters from 0 to 9, point and space) */
  for (i = 0; i < size_dict; i++) {
    dict[i] = malloc(2 * sizeof(char));
    dict[i][0] = '\0';
    dict[i][1] = '\0';
  }
  for (i = 0; i < length_dictionary; i++) {
    dict[(int) dictionary[i]] = malloc(2 * sizeof(char));
    dict[(int) dictionary[i]][0] = dictionary[i];
    dict[(int) dictionary[i]][1] = '\0';
  }

  /* expand the compressed data */
  if ((temp_file = tmpfile()) == NULL) { /* open a temporary file */
    fprintf(stderr, "Problem opening temporary file...Aborting\n");
    exit(1);
  }
  index = *bytes_ptr;
  bytes_ptr++;
  current_num_value++;
  fprintf(temp_file, "%s", dict[index]);
  old = index;
  while (current_num_value < num_lines) { /* main loop */
    index = *bytes_ptr;
    bytes_ptr++;
    if (dict[index][0] != '\0') { /* index is in the dictionary */
      fprintf(temp_file, "%s", dict[index]);
      B = dict[index][0];
      dict_pos++;

      /* add the string to the dictionary if there is room */
      if (dict_pos < size_dict) {
        free(dict[dict_pos]);
        string_length = strlen(dict[old]);
        dict[dict_pos] = malloc(sizeof(char) * (string_length + 2));
        strcpy(dict[dict_pos], dict[old]);
        dict[dict_pos][string_length] = B;
        dict[dict_pos][string_length + 1] = '\0';
      }
    }
    else { /* index is not in the dictionary */
      B = dict[old][0];
      dict_pos++;
      free(dict[dict_pos]);
      string_length = strlen(dict[old]);
      dict[dict_pos] = malloc(sizeof(char) * (string_length + 2));
      strcpy(dict[dict_pos], dict[old]);
      dict[dict_pos][string_length] = B;
      dict[dict_pos][string_length + 1] = '\0';
      fprintf(temp_file, "%s", dict[dict_pos]);
    }

    /* prepare to read the next index */
    current_num_value++;
    old = index;
  }

  /* fill in the energy values array */
  i = 0;
  rewind(temp_file);
  while (!feof(temp_file)) {
    fscanf(temp_file, "%lf\n", &energy_values[i]);
    i++;
    if (i == num_values - 1)
      break;
  }

  /* clean up the memory */
  for (i = 0; i < size_dict; i++)
    free(dict[i]);

  fclose(temp_file);
  return energy_values;
}

/*---------------------------------------------------------------------------*/

void filter_easymifs_map(const char *filename, const double energy_cutoff)
{
   /* read in an easymifs map and writes a new file with only the filtered
     points (whose energy is below the user-defined threshold */

  FILE *infile, *outfile;
  char *cp, outfilename[MAX_STR_LEN];
  char line[MAX_STR_LEN];
  char *temp = NULL;
  char open_mode[MAX_STR_LEN];
  bool end_of_header = FALSE;
  bool compressed = FALSE; /* whether we are using compressed maps */
  double coords[3]; /* coordinates of a point */
  double energy;
  double *energy_values = NULL;
  unsigned int cx, cy, cz; /* counter variables for computing the coordinates */
  unsigned int index_point = 1; /* keep track of the points */
  unsigned int counter = 0;
  unsigned int total_num_points;
  unsigned int *bytes = NULL; /* to store the compressed map */
  unsigned int *bytes_ptr;
  int size = 0;
  int byte = 0;
  const unsigned int num_bytes = BITS / 8;

  resolution = -1.000;
  /* if the extension of the file is .cmp then assume it is compressed */
  temp = strrchr(filename, '.');
  if (temp != NULL) {
    if (strcmp(temp, ".cmp") == 0) {
      compressed = TRUE;
      strncpy(open_mode, "rb", MAX_STR_LEN); /* set the open mode to binary */
    }
    else {
      strncpy(open_mode, "r", MAX_STR_LEN);
    }
  }

  /* open/check the file to be filtered */
  if (!(infile = fopen(filename, open_mode))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", filename);
    exit(1);
  }

  /* open the output file */
  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  if (decompress)
    strcat(outfilename, ".dx");
  else
    strcat(outfilename, ".tmp");
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s for writing\nAborting...\n", outfilename);
  }

  /* make sure the map is compressed before attempting the decompression */
  if (decompress && !compressed) {
    fprintf(stderr, "The map does not seem to be compressed\nAborting...\n");
    exit(1);
  }
  
  /* read the header of the easymifs file */
  while ( (!end_of_header) && fgets(line, MAX_STR_LEN, infile)) {
    if (decompress) { /* print the line and go on */
      fprintf(outfile, "%s", line);
      if (strstr(line, "data follows")) /* end of the header */
      end_of_header = TRUE;
    }

    if (strstr(line, "counts")) { /* grab the number of points in x, y, z */
      strtok(line, " ");
      strtok(NULL, " "); strtok(NULL, " "); strtok(NULL, " ");
      temp = strtok(NULL, " ");

      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[X]);
      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[Y]);
 
      temp = strtok(NULL, " "); sscanf(temp, "%d", &npoints[Z]);
    }

    if (strstr(line, "origin")) { /* grab the coordinates of the origin */
      strtok(line, " ");
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &lower_corner[X]);
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &lower_corner[Y]);
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &lower_corner[Z]);
    }

    if (strstr(line, "delta") && (resolution < 0.0)) { 
    /* grab the resolution of the grid */
      strtok(line, " ");
      temp = strtok(NULL, " "); sscanf(temp, "%lf", &resolution);
    }
   
    if (strstr(line, "data follows")) /* end of the header */
      end_of_header = TRUE;
  }
  
  if (! end_of_header) {
    fprintf(stderr, "Problem with easymifs file\nAborting...\n");
    exit(1);
  }

  /* calculate the total number of points */
  total_num_points = npoints[X] * npoints[Y] * npoints[Z];
  counter = 0;
  if (compressed) { /* read the bytes and expand the map */
    bytes = malloc(sizeof(unsigned int) * total_num_points * 2);
    bytes_ptr = bytes;
    while (!feof(infile)) {
      size = fread(&byte, num_bytes, 1, infile);
      if (size > 0) {
        counter++;
        *bytes_ptr = byte;
        bytes_ptr++;
      }
    }
    energy_values = expand(bytes, counter, total_num_points);
  }
   /* read the data from the file (or use the values stored in energy_values) 
      and write the filtered points to the output file */
    counter = 0;
    for (cx = 0; cx < npoints[X]; cx++) {
      coords[X] = lower_corner[X] + resolution * cx;
      for (cy = 0; cy < npoints[Y]; cy++) {
        coords[Y] = lower_corner[Y] + resolution * cy;
        for (cz = 0; cz < npoints[Z]; cz++) {
          coords[Z] = lower_corner[Z] + resolution * cz;
          if (compressed) {
            energy = energy_values[counter];
            counter++;
	  }
          else {
            if (fgets(line, MAX_STR_LEN, infile) == NULL) {
              fprintf(stderr, "Problems with EasyMIFs file\n...Aborting\n");
              exit(1);
            }
            sscanf(line, "%lf", &energy); /* read the energy value */
	  }
          if (decompress) { /* write the energy value */
            fprintf(outfile, "%.3f\n", energy);
	  }
          else if (energy <= energy_cutoff) { /* if energy < cutoff write point */
	    fprintf(outfile, "%d %.3f %.3f %.3f %.3f\n", index_point,
                    coords[X], coords[Y], coords[Z], energy);
          }
          index_point++;
        }
      }
    }

  /* clean up the memory */
  if (compressed) {
    free(bytes);
    free(energy_values);
  }

  /* close input and output files */
  fclose(infile);
  fclose(outfile);
}

/*---------------------------------------------------------------------------*/

void cluster_grid(const char *filename, const double cutoff, 
                  const char *linkage)
{
  /* pilot function that performs the clustering of the grid points */

  Data *data;
  double **distance_matrix;
  int *clusters;
  Node *results = NULL;
  Summary *summary;
  char *cp, infilename[MAX_STR_LEN];
  int n = 0; /* number of grid points */
  int nc; /* number of clusters */
  int i;
  
  /* get the filename */
  cp = strdup(filename);
  strcpy(infilename, basename(cp));
  strcat(infilename, ".tmp");

  /* load the data matrix */
  data = create_data_matrix((const char *) infilename, &n);

  /* compute the distance matrix */
  distance_matrix = create_distance_matrix(data, n);

  /* perform the clustering */
  fprintf(stdout, "   Performing clustering with %s linkage...\n", linkage);
  fprintf(stdout, "   Energy cutoff: %.3f\n", energy_cutoff); 
  fprintf(stdout, "   Spatial cutoff: %.3f\n", cutoff);
  if (strcmp(linkage, "single") == 0)
    results = sl_cluster(distance_matrix, n);
  else if(strcmp(linkage, "average") == 0)
    results = average_cluster(distance_matrix, n);

  /* free the allocated memory for the distance matrix */
  for (i = 0; i < n; i++)
    free(distance_matrix[i]);
  free(distance_matrix);

  /* cut the tree at the specified cutoff level */
  clusters = cuttree(results, n, cutoff, &nc);
  free(results);

  /* rank the points according to total cluster interaction energy */
  rank_data(data, clusters, nc, n);
  free(clusters);

  /* print the results */
  print_clusters(data, n, filename);

  /* create and print the cluster summary */
  summary = create_summary(data, n, nc);
  print_summary(summary, nc, filename);
  free(data);
  free(summary);

  fprintf(stdout, "   ...done\n");
}

/*---------------------------------------------------------------------------*/

Data *create_data_matrix(const char *filename, int *n)
{
  /* create a data matrix suitable for computing the distance matrix */
  FILE *infile;
  int i;
  char s[MAX_STR_LEN], *temp;
  Data *data;
  double sign; /* sign for the random jiggling */
  double random;

  if (!(infile = fopen(filename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", filename);
    exit(1);
  }

  /* count how many points the file contains */
  while (fgets(s, MAX_STR_LEN - 1, infile))
    (*n)++;

  /* allocate memory */
  data = (Data *) malloc(sizeof(Data) * (*n));

  /* fill in the data */
  rewind(infile);
  for (i = 0; i < *n; i++) {
    fgets(s, MAX_STR_LEN - 1, infile);
    temp = strtok(s, " ");
    sscanf(temp, "%lu", &data[i].index);
    data[i].x = strtod(strtok(NULL, " "), NULL);
    data[i].y = strtod(strtok(NULL, " "), NULL);
    data[i].z = strtod(strtok(NULL, " "), NULL);
    if (strcmp(linkage, "average") == 0) { /* jiggle the points to break ties */
      random = ran(&seed);
      sign = (random > 0.5) ? 1.0: -1.0;
      data[i].x += sign * ran(&seed) * 1.0E-6;
      data[i].y += sign * ran(&seed) * 1.0E-6;
      data[i].z += sign * ran(&seed) * 1.0E-6;
    }
    data[i].energy = strtod(strtok(NULL, " "), NULL);
  }

  fclose(infile);
  return data;
}

/*---------------------------------------------------------------------------*/

double **create_distance_matrix(Data * data, int n)
{
  /* create the distance matrix */

  int i, j;
  double **matrix;

  matrix = malloc(sizeof(double *) * n);
  matrix[0] = NULL;

  /* The zeroth row has zero columns. 
     We allocate it anyway for convenience. */
  for (i = 1; i < n; i++) {
    matrix[i] = malloc(sizeof(double) * i);
    if (matrix[i] == NULL)
      break; /* Not enough memory available */
  }
  if (i < n) { /* break condition encountered */
    j = i;
    for (i = 1; i < j; i++)
      free(matrix[i]);

    return NULL;
  }

  /* Calculate the distances and save them in the matrix */
  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++)
      matrix[i][j] = sqrt(pow(data[i].x - data[j].x, 2) +
                          pow(data[i].y - data[j].y, 2) +
			  pow(data[i].z - data[j].z, 2));

  return matrix;
}

/*---------------------------------------------------------------------------*/

double find_closest_pair(int n, double **distmatrix, int *ip, int *jp)
{

  int i, j;
  double distance = distmatrix[1][0];
  for (i = 0; i < n; i++)
  { for (j = 0; j < i; j++)
    { if (distmatrix[i][j]<distance)
      { distance = distmatrix[i][j];
        *ip = i;
        *jp = j;
      }
    }
  }
  return distance;
}

/*---------------------------------------------------------------------------*/

Node *sl_cluster(double **distmatrix, int n)
{
  /* perform single-linkage clustering using the SLINK algorithm */

  int i, j, k;
  const int nnodes = n - 1;
  int *vector;
  double *temp;
  int *index;
  Node *result;

  /* memory allocation */
  temp = malloc(sizeof(double) * nnodes);
  if (!temp)
    return NULL;
  index = malloc(sizeof(int) * n);
  if (!index) {
    free(temp);
    return NULL;
  }
  vector = malloc(sizeof(int) * nnodes);
  if (!vector) {
    free(index);
    free(temp);
    return NULL;
  }
  result = malloc(sizeof(Node) * nnodes);
  if (!result) {
    free(vector);
    free(index);
    free(temp);
    return NULL;
  }

  for (i = 0; i < nnodes; i++) {
    vector[i] = i;
    result[i].distance = DBL_MAX;
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++)
      temp[j] = distmatrix[i][j];
    for (j = 0; j < i; j++) {
      k = vector[j];
      if (result[j].distance >= temp[j]) {
        if (result[j].distance < temp[k])
	  temp[k] = result[j].distance;

	 result[j].distance = temp[j];
	 vector[j] = i;
      }
      else if (temp[j] < temp[k])
        temp[k] = temp[j];
    }
    for (j = 0; j < i; j++)
      if (result[j].distance >= result[vector[j]].distance)
        vector[j] = i;
  }
  free(temp);

  for (i = 0; i < nnodes; i++)
    result[i].left = i;

  qsort(result, nnodes, sizeof(Node), compare_nodes);

  for (i = 0; i < n; i++)
    index[i] = i;
  for (i = 0; i < nnodes; i++) {
    j = result[i].left;
    k = vector[j];
    result[i].left = index[j];
    result[i].right = index[k];
    index[k] = -i - 1;
  }
  free(vector);
  free(index);

  return result;
}

/*---------------------------------------------------------------------------*/

Node *average_cluster(double **distmatrix, int nelements)
{
  /* perform average-linkage clustering */
  int j;
  int n;
  int* clusterid;
  int* number;
  Node* result;

  clusterid = malloc(nelements*sizeof(int));
  if(!clusterid) return NULL;
  number = malloc(nelements*sizeof(int));
  if(!number)
  { free(clusterid);
    return NULL;
  }
  result = malloc((nelements-1)*sizeof(Node));
  if (!result)
  { free(clusterid);
    free(number);
    return NULL;
  }

  /* Setup a list specifying to which cluster a gene belongs, and keep track
   * of the number of elements in each cluster (needed to calculate the
   * average). */
  for (j = 0; j < nelements; j++) {
    number[j] = 1;
    clusterid[j] = j;
  }

  for (n = nelements; n > 1; n--) {
    int sum;
    int is = 1;
    int js = 0;
    result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

    /* Save result */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];

    /* Fix the distances */
    sum = number[is] + number[js];
    for (j = 0; j < js; j++) {
      distmatrix[js][j] = distmatrix[is][j]*number[is]
                        + distmatrix[js][j]*number[js];
      distmatrix[js][j] /= sum;
    }
    for (j = js+1; j < is; j++) {
      distmatrix[j][js] = distmatrix[is][j]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }
    for (j = is+1; j < n; j++) {
      distmatrix[j][js] = distmatrix[j][is]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update number of elements in the clusters */
    number[js] = sum;
    number[is] = number[n-1];

    /* Update clusterids */
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
  }
  free(clusterid);
  free(number);

  return result;
}

/*---------------------------------------------------------------------------*/

int *cuttree(Node * tree, int nelements, double cutoff, int *nclusters)
{
  int i, j, k;
  int icluster = 0;
  int n; /* number of nodes to join */
  int *nodeid;
  int *clusterid;

  *nclusters = 1;
  /* calculate how many clusters there will be */
  for (i = nelements - 2; i >= 0; i--)
    if (tree[i].distance > cutoff)
      (*nclusters)++;
    else
      break;
  n = nelements - *nclusters;
  clusterid = (int *) malloc(sizeof(int) * nelements);
  for (i = nelements - 2; i >= n; i--) {
    k = tree[i].left;
    if (k >= 0) {
      clusterid[k] = icluster;
      icluster++;
    }
    k = tree[i].right;
    if (k >= 0) {
      clusterid[k] = icluster;
      icluster++;
    }
  }
  nodeid = malloc(n * sizeof(int));
  if (!nodeid) {
    for (i = 0; i < nelements; i++)
      clusterid[i] = -1;
    return clusterid;
  }

  for (i = 0; i < n; i++)
    nodeid[i] = -1;
  for (i = n - 1; i >= 0; i--) {
    if (nodeid[i] < 0) {
      j = icluster;
      nodeid[i] = j;
      icluster++;
    }
    else
      j = nodeid[i];

    k = tree[i].left;
    if (k < 0)
      nodeid[-k - 1] = j;
    else
      clusterid[k] = j;

    k = tree[i].right;
    if (k < 0)
      nodeid[-k - 1] = j;
    else
      clusterid[k] = j;
  }

  free(nodeid);
  return clusterid;
}

/*---------------------------------------------------------------------------*/

void rank_data(Data * data, const int *clusters, const int nc, const int n)
{
  double total_energies[nc];
  int i;

  for (i = 0; i < nc; i++)
    total_energies[i] = 0.0;

  /* compute the total clusters interaction energy */
  for (i = 0; i < n; i++)
    total_energies[clusters[i]] += data[i].energy;

  /* attach to each data point the total interaction energy of its cluster */
  for (i = 0; i < n; i++) {
    data[i].cl_energy = total_energies[clusters[i]];
    data[i].cl_number = clusters[i];
  }

  /* sort the data points */
  qsort(data, n, sizeof(Data), compare_data);
}

/*---------------------------------------------------------------------------*/

Summary *create_summary(Data * data, const int n, const int nc)
{
  Summary *summary;
  int i, s = 0, cn;
  double weight; /* weighting factor for the coordinates */

  summary = (Summary *) malloc(sizeof(Summary) * nc);

  /* first data point */
  summary[s].cl_energy = data[0].cl_energy;
  weight = data[0].energy / data[0].cl_energy;
  summary[s].x = data[0].x * weight;
  summary[s].y = data[0].y * weight;
  summary[s].z = data[0].z * weight;
  summary[s].nmembers = 1;
  cn = data[0].cl_number;

  for (i = 1; i < n; i++) { /* all data points but the first */
    if (cn != data[i].cl_number) { /* new cluster info */
      s++;
      cn = data[i].cl_number;
      summary[s].cl_energy = data[i].cl_energy;
      weight = data[i].energy / data[i].cl_energy;
      summary[s].x = data[i].x * weight;
      summary[s].y = data[i].y * weight;
      summary[s].z = data[i].z * weight;
      summary[s].nmembers = 1;
    }
    else { /* update the cluster "center of energy" */
      weight = data[i].energy / data[i].cl_energy;
      summary[s].x += data[i].x * weight;
      summary[s].y += data[i].y * weight;
      summary[s].z += data[i].z * weight;
      summary[s].nmembers++;
    }
  }
  return summary;
}

/*---------------------------------------------------------------------------*/

int compare_nodes(const void *a, const void *b)
{
  /* helper function for qsort */

  const Node *node1 = (const Node *) a;
  const Node *node2 = (const Node *) b;
  const double term1 = node1->distance;
  const double term2 = node2->distance;
  if (term1 < term2)
    return -1;

  if (term1 > term2)
    return +1;

  return 0;
}

/*---------------------------------------------------------------------------*/

int compare_data(const void *a, const void *b)
{
  /* helper function for qsort */
  double energy;
  const Data *da = (const Data *) a;
  const Data *db = (const Data *) b;


  energy = ((*da).cl_energy > (*db).cl_energy) - ((*da).cl_energy <
                                                  (*db).cl_energy);

  if (energy > 0 || energy < 0 || (*da).cl_number == (*db).cl_number)
    return energy;
  else
    return ((*da).cl_number > (*db).cl_number) - ((*da).cl_number <
                                                  (*db).cl_number);
}

/*---------------------------------------------------------------------------*/

void print_clusters(Data * data, const int n, const char *filename)
{
  char outfilename[MAX_STR_LEN], *cp;
  FILE *outfile;
  int i, index, cl_number;

  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, "_clusters.dat");

  /* write the results on file */
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", outfilename);
    exit(1);
  }
  fprintf(outfile, "%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%lu\n", 1,
          data[0].cl_energy, data[0].energy, data[0].x, data[0].y,
	  data[0].z, data[0].index); /* first value */
  index = data[0].cl_number;
  cl_number = 1;

  for (i = 1; i < n; i++) { /* remaining values */
    if (index != data[i].cl_number) {
      index = data[i].cl_number;
      cl_number++;
    }

    fprintf(outfile, "%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%lu\n",
            cl_number, data[i].cl_energy, data[i].energy, data[i].x,
            data[i].y, data[i].z, data[i].index);
  }
  fclose(outfile);
}

/*---------------------------------------------------------------------------*/

void print_summary(Summary * summary, const int nc, const char *filename)
{
  char outfilename[MAX_STR_LEN], *cp;
  FILE *outfile;
  int i;

  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, "_summary.dat");

  /* write the summary on file */
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", outfilename);
    exit(1);
  }

  for (i = 0; i < nc; i++)
    fprintf(outfile, "%d\t%8.3f\t%d\t%8.3f\t%8.3f\t%8.3f\n", i + 1,
      summary[i].cl_energy, summary[i].nmembers, \
      summary[i].x, summary[i].y, summary[i].z);

  fclose(outfile);
}

/*---------------------------------------------------------------------------*/

void create_pdb(const char *filename)
{
  FILE *infile, *outfile;
  char infilename[MAX_STR_LEN], outfilename[MAX_STR_LEN], *cp, 
       line[MAX_STR_LEN];
  int cn;
  char pdb_line[] = "ATOM  %5d  H   %s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
  char resname[4];
  int resnum; /* the number of the current cluster */
  char chain = ' '; /* the chain ID that will be used for the clusters */
  double tie, ie, x, y, z; /* total/point-wise interaction energy, coords */
  double max_tie; /* used for the normalization */
  unsigned int i, j;
  bool found = FALSE; /* flag for finding unused chain identifier */

  /* pick a suitable chain ID for the clusters */
  for (i = 0; i < TOTAL_CHAIN_IDS; i++) {
    found = FALSE;
    for (j = 0; j < TOTAL_CHAIN_IDS; j++) {
      if (used_chains[j] == '\0')
	break;
      if (used_chains[j] == chain_IDs[i]) {
	found = TRUE;
	break;
      }
    }
    if (!found) {
      chain = chain_IDs[i];
      break;
    }
  }

  /* open the clusters file */
  cp = strdup(filename);
  strcpy(infilename, basename(cp));
  strcat(infilename, "_clusters.dat");
  if (!(infile = fopen(infilename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", infilename);
    exit(1);
  }

  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, "_clusters.pdb");

  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", outfilename);
    fclose(infile);
    exit(1);
  }

  /* write the clusters on the pdb file, one residue per cluster
     with a normalized b-factor (the TIE of each cluster  is divided by the
     TIE of the first ranking cluster */

  fgets(line, MAX_STR_LEN, infile);	/* get the maximum TIE and store it */
  cp = strdup(line);
  cn = atoi(strtok(cp, "\t"));
  max_tie = fabs(strtod(strtok(NULL, "\t"), NULL));
  rewind(infile);

  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    cp = strdup(line);
    cn = atoi(strtok(cp, "\t"));
    cluster_atom_num++; /* increment the cluster atom num */
    resnum = cn; /* residue number is equal to cluster number */
    if (cn >= MAX_CL_NUM)
      break;

    /* normalize tie */
    tie = 100 * fabs(strtod(strtok(NULL, "\t"), NULL)) / max_tie;
    ie = fabs(strtod(strtok(NULL, "\t"), NULL));
    x = strtod(strtok(NULL, "\t"), NULL);
    y = strtod(strtok(NULL, "\t"), NULL);
    z = strtod(strtok(NULL, "\t"), NULL);

    if (cn < 10)
      sprintf(resname, "C0%d", cn);
    else
      sprintf(resname, "C%d", cn);

    fprintf(outfile, pdb_line, cluster_atom_num, resname, chain, resnum, x, y, z, ie, tie);
  }

  fclose(infile);
  fclose(outfile);

}

/*---------------------------------------------------------------------------*/

void create_dx(const char *filename)
{
  /* print a .dx file that contains all the clustering info */

  char infilename[MAX_STR_LEN], outfilename[MAX_STR_LEN], *cp;
  char line[MAX_STR_LEN];
  FILE *infile, *outfile;
  unsigned long int i = 0, index = 1; /* indices */
  unsigned long total_points; /* total number of points */
  unsigned int cx, cy, cz; /* counter for x, y, z */
  double tie; /* total interaction energy of the cluster */
  double *points_tie; /* array with the TIE values of the cluster the point
                         belongs to */

  /* allocate the memory for points_tie */
  total_points = npoints[X] * npoints[Y] * npoints[Z];
  points_tie = malloc(sizeof(double) * (total_points + 1));
  for (i = 0; i <= total_points; i++)
    points_tie[i] = 0.0;

  /* open the output file */
  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, "_clusters.dx");
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", outfilename);
    exit(1);
  }

  /* write the header */
  fprintf(outfile, "# Sitehound output\n#\n#\n#\n");
  fprintf(outfile, "object 1 class gridpositions counts %d %d %d\n", npoints[X],
          npoints[Y], npoints[Z]);
  fprintf(outfile, "origin %.3f\t%.3f\t%.3f\n", lower_corner[X], 
          lower_corner[Y], lower_corner[Z]);
  fprintf(outfile, "delta %.3f 0 0\ndelta 0 %.3f 0\ndelta 0 0 %.3f\n",
          resolution, resolution, resolution);
  fprintf(outfile, "object 2 class gridconnections counts %d %d %d\n", 
          npoints[X], npoints[Y], npoints[Z]);
  fprintf(outfile, "object 3 class array type double rank 0 items %lu ", 
          total_points);
  fprintf(outfile, "data follows\n");
 
  /* read the data from the cluster file */
  cp = strdup(filename);
  strcpy(infilename, basename(cp));
  strcat(infilename, "_clusters.dat");
  if (!(infile = fopen(infilename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", infilename);
    exit(1);
  }
  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    strtok(line, "\t"); /* first field (cluster number) */
    tie = strtod(strtok(NULL, "\t"), NULL); /* tie field */
    strtok(NULL, "\t"); /* ie field */
    strtok(NULL, "\t"); strtok(NULL, "\t"); strtok(NULL, "\t"); /* x, y, z */
    sscanf(strtok(NULL, "\t"), "%lu", &index); /* index field */
    points_tie[index] = tie; /* assign the tie of the cluster to the
                                corresponding point */
  }
  fclose(infile);
  
  /* write the output according to the DX convention (Z, Y, X) */
  if (strcmp(maptype, "autogrid") == 0) {
    for (cx = 0; cx < npoints[X]; cx++)
      for (cy = 0; cy < npoints[Y]; cy++)
        for (cz = 0; cz < npoints[Z]; cz++) {
          i++;
	  index = cz * (npoints[X] * npoints[Y]) + cy * npoints[X] + cx + 1;
          fprintf(outfile, "%.3f\n", points_tie[index]);
        }
  }
  else if (strcmp(maptype, "easymifs") == 0) {
    index = 1;
    for (cx = 0; cx < npoints[X]; cx++)
      for (cy = 0; cy < npoints[Y]; cy++)
        for (cz = 0; cz < npoints[Z]; cz++) {
          fprintf(outfile, "%.3f\n", points_tie[index]);
	  index++;
        }
  }
    
  free(points_tie);
  fclose(outfile);
  
}

/*---------------------------------------------------------------------------*/

void create_contacts(const char *filename)
{
  /* prints a file where each line contains the residues that are in contact
     with the corresponding cluster */

  char infilename[MAX_STR_LEN], outfilename[MAX_STR_LEN], stem[MAX_STR_LEN],*cp;
  char line[MAX_STR_LEN], temp[MAX_STR_LEN];
  char resname[MAX_STR_LEN];
  char *word;
  FILE *infile, *outfile;
  Atom *atom_list = NULL, *current_atom;
  Data p; /* current cluster point */
  unsigned int current_cl_num, old_cluster_number;
  unsigned int counter, i, n;
  Restable *first_residue = NULL, *current_residue;
  bool flag; /* flag to check if a residue has been written already */
  bool exit_flag = FALSE;

  /* store the coordinates and residue names of the PDB residues */
  current_atom = atom_list;
  cp = strdup(filename);
  strcpy(stem, basename(cp));
  strcpy(infilename, stem);
  strcat(infilename, ".pdb");
  if (!(infile = fopen(infilename, "r"))) {
    strcpy(infilename, stem);
    strcat(infilename, ".ent");
    if (!(infile = fopen(infilename, "r"))) {
      /* count how many underscore-separated tokens are present in the name */
      n = 0;
      cp = strdup(filename);
      word = strtok(cp, "_");
      if (word != NULL)
        n++;
      while (strtok(NULL, "_") != NULL)
        n++;
      cp = strdup(filename);
      strcpy(stem, strtok(cp, "_"));
      for (counter = 1; counter < n - 1; counter++) {
        word = strtok(NULL, "_");
        strcat(stem, "_");
        strcat(stem, word);
      }
      strcpy(infilename, stem);
      strcat(infilename, ".pdb");
      if (!(infile = fopen(infilename, "r"))) {
        strcpy(infilename, stem);
        strcat(infilename, ".ent");
      }
      if (!(infile = fopen(infilename, "r"))) {
        fprintf(stderr, "Cannot open %s\nAborting...\n", infilename);
        exit(1);
      }
    }
  }
  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {
      /* allocate memory for new atom */
      current_atom = malloc(sizeof(Atom));

      /* extract the atom serial number */
      strncpy(temp, line + PDB_ATOMID_SHIFT, 5);
      temp[5] = '\0';
      sscanf(temp, "%u", &current_atom->atomNum);

      /* extract the residue name */
      strncpy(current_atom->resID, line + PDB_RESID_SHIFT, 3);
      current_atom->resID[3] = '\0';

      /* extract the chain identifier */
      current_atom->chain = line[PDB_CHAIN_SHIFT];

      /* check if the chain has been already stored */
      for (i = 0; i < TOTAL_CHAIN_IDS; i++) {
        if (used_chains[i] == current_atom->chain) /* chain already stored */
	  break;
        if (used_chains[i] == '\0') { /* chain not yet stored */
	  used_chains[i] = current_atom->chain;
          break;
         }
      }

      /* extract the residue number */
      strncpy(temp, line + PDB_RESNUM_SHIFT, 4);
      temp[4] = '\0';
      sscanf(temp, "%u", &current_atom->resNum);

      /* extract the coordinates */
      strncpy(temp, line + PDB_X_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->x);
      strncpy(temp, line + PDB_Y_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->y);
      strncpy(temp, line + PDB_Z_SHIFT, 8);
      temp[8] = '\0';
      sscanf(temp, "%lf", &current_atom->z);

      current_atom->next = atom_list;
      atom_list = current_atom;
    }
    if (strncmp(line, "TER", 3) == 0) {
      strncpy(temp, line + PDB_ATOMID_SHIFT, 5);
      temp[5] = '\0';
      sscanf(temp, "%u", &cluster_atom_num);
    }
  }
  fclose(infile);

  /* store the last atom number */
  if (cluster_atom_num < current_atom->atomNum)
    cluster_atom_num = current_atom->atomNum;

  /* open the output file */
  cp = strdup(filename);
  strcpy(outfilename, basename(cp));
  strcat(outfilename, "_predicted.dat");
  if (!(outfile = fopen(outfilename, "w"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", outfilename);
    exit(1);
  }

  /* process the cluster file and print the contact residue list */
  cp = strdup(filename);
  strcpy(infilename, basename(cp));
  strcat(infilename, "_clusters.dat");
  if (!(infile = fopen(infilename, "r"))) {
    fprintf(stderr, "Cannot open %s\nAborting...\n", infilename);
    exit(1);
  }
  old_cluster_number = 0;
  while (fgets(line, MAX_STR_LEN, infile) != NULL) {
    /* read in the current cluster point */
    sscanf(line, "%d %lf %lf %lf %lf %lf", &current_cl_num, &p.cl_energy, 
           &p.energy, &p.x, &p.y, &p.z);
    if (current_cl_num > cl_num_contacts) { /* all relevant residues listed */
      print_residues(outfile, first_residue);
      free_memory_restable(first_residue);
      exit_flag = TRUE;
      break;
    }
    if (old_cluster_number < current_cl_num) { /* write new cluster index */
      if (current_cl_num == 1)
        fprintf(outfile, "%d", current_cl_num);
      else {
        print_residues(outfile, first_residue);
        free_memory_restable(first_residue);
	fprintf(outfile, "\n%d", current_cl_num);
      }
      first_residue = NULL;
    }
    old_cluster_number = current_cl_num;

    /* check which residues are in contact */
    current_atom = atom_list;
    counter = 0; /* initialize counter of residues */
    while (current_atom != NULL) {
      char residue_code[MAX_STR_LEN] = {'\0'};
      if (euclidean_distance(&p, current_atom) <= DIST_CONTACTS) {
        /* assemble the residue code */
	strcpy(residue_code, current_atom->resID);
        sprintf(resname, "%d", current_atom->resNum);
        strcat(residue_code, "_");
        strcat(residue_code, resname);
        strcat(residue_code, "_");
	residue_code[strlen(residue_code)] = current_atom->chain;
	residue_code[strlen(residue_code)] = '\0';
        flag = TRUE;

        /* compare it to the previous ones */
        current_residue = first_residue;
        while (current_residue != NULL) {
	  if (strcmp(current_residue->resname, residue_code) == 0) {
	    flag = FALSE;
	    break;
	  }
	  current_residue = current_residue->next;
	}

        if (flag) { /* the residue has never been written for this cluster */
          /* add the residue to the list */
	  current_residue = first_residue;
	  current_residue = malloc(sizeof(Restable));
	  strcpy(current_residue->resname, residue_code);
          current_residue->printed = FALSE;
	  current_residue->next = first_residue;
	  first_residue = current_residue;
	}
      }
      current_atom = current_atom->next;
    }
  }
  if (!exit_flag) {
    print_residues(outfile, first_residue);
    free_memory_restable(first_residue);
  }
  

  fclose(infile);
  fclose(outfile);
}

/*---------------------------------------------------------------------------*/

void print_residues(FILE *outfile, Restable *residue_list)
{
  /* print the residues sorted by chain and residue number */

  Restable *current_residue, *residue_to_print;
  unsigned int current_resnum, to_print_resnum;
  char current_chain[MAX_STR_LEN];
  char temp[MAX_STR_LEN];
  char to_print_chain[MAX_STR_LEN];
  bool printed_all; /* to check if all residues have been printed */
  

  while (1) {
    printed_all = TRUE;
    current_residue = residue_list;
    residue_to_print = NULL;

    while (current_residue != NULL) { /* screen all the residues */
      if (residue_to_print == NULL && !current_residue->printed) {
        /* assign the first non-printed residue to residue_to_print */
        residue_to_print = current_residue;
        current_residue = current_residue->next;
        printed_all = FALSE;
        continue; /* move to the next residue */
      }

      if (residue_to_print != NULL && !current_residue->printed) {
        /* check if current_residue comes before residue_to_print */
        strcpy(temp, current_residue->resname);
        strtok(temp, "_");
        current_resnum = atoi(strtok(NULL, "_"));
        strcpy(current_chain, strtok(NULL, "_"));

        strcpy(temp, residue_to_print->resname);
        strtok(temp, "_");
        to_print_resnum = atoi(strtok(NULL, "_"));
        strcpy(to_print_chain, strtok(NULL, "_"));

        if (strcmp(current_chain, to_print_chain) < 0) { /* compare chains */
          residue_to_print = current_residue;
          current_residue = current_residue->next;
          continue; /* move to the next residue */
        }
        else if (strcmp(current_chain, to_print_chain) == 0) { /* compare num */
          if (current_resnum < to_print_resnum) {
	    residue_to_print = current_residue;
            current_residue = current_residue->next;
            continue; /* move to the next residue */
	  }
        }
      }
      current_residue = current_residue->next;
    }
    /* check for printing the residue */
    if (residue_to_print != NULL) {
      fprintf(outfile, " %s", residue_to_print->resname);
      residue_to_print->printed = TRUE;
    }
    if (printed_all)
      break;
  }
}

/*---------------------------------------------------------------------------*\
 * MAIN PROGRAM                                                              *
\*---------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
  /* print a welcome message */
  fprintf(stdout, "\n***********************************************\n");
  fprintf(stdout, "* Welcome to sitehound, version: %s       *\n", VERSION);
  fprintf(stdout, "***********************************************\n\n");

  /* get the command-line parameters */
  set_parameters(argc, argv);

  if (!strcmp(filename, "-1")) {
    fprintf(stderr, "Please enter at least the name of the map\n");
    fprintf(stderr, "Usage: %s\n", USAGE);
    exit(1);
  }

  /* decompress a .cmp file back to .dx */
  if (decompress) {
    fprintf(stdout, "Decompressing the map\n");
    filter_easymifs_map(filename, energy_cutoff);
    return 0;
  }

  /* filter out the points above the energy threshold */
  if (strcmp(maptype, "autogrid") == 0)
    filter_autogrid_map(filename, energy_cutoff);
  else if (strcmp(maptype, "easymifs") == 0)
    filter_easymifs_map(filename, energy_cutoff);
  else {
    fprintf(stderr, "Please specify the maptype ('easymifs' or 'autogrid')\n");
    exit(1);
  }

  /* cluster the remaining points */
  cluster_grid(filename, spatial_cutoff, linkage);

  /* create the contact information file */
  create_contacts(filename);
    
  /* print the PDB file */
  create_pdb(filename);
  
  /* print the DX file */
  create_dx(filename);

  return 0;
}
