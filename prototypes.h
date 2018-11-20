/*****************************************************************************
 *****************************************************************************
 ** prototypes.h by Dario Ghersi                                            **
 ** Version: 20101214                                                       **
 ** (part of SITEHOUND)                                                     **
 *****************************************************************************
 *****************************************************************************/

/*---------------------------------------------------------------------------*\
 *  CONSTANTS                                                                *
\*---------------------------------------------------------------------------*/

#define VERSION "20101214"
#define USAGE "sitehound -f=MAP [-t=TYPE -l=LINKAGE -s=SPATIAL_CUTOFF -e=ENERGY_CUTOFF -p=CL_NUM -z]"
#define MAX_STR_LEN 100 /* maximum length for a string */
#define MAX_CL_NUM 100 /* the maximum number of clusters to plot to PDB */
#define DIST_CONTACTS 16.0 /* the squared distance that defines whether a 
                              residue is in contact with a cluster */
#define BITS 16 /* size of the code for compression (byte multiple) */

#define X 0
#define Y 1
#define Z 2

/* PDB FORMAT CONSTANTS */
#define PDB_ATOMID_SHIFT 6
#define PDB_RESID_SHIFT 17
#define PDB_CHAIN_SHIFT 21
#define PDB_RESNUM_SHIFT 22
#define PDB_X_SHIFT 30
#define PDB_Y_SHIFT 38
#define PDB_Z_SHIFT 46
#define TOTAL_CHAIN_IDS 37
/* all the possible PDB chain identifiers */
const char chain_IDs[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
                          'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
			  'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3',
                          '4', '5', '6', '7', '8', '9', ' '};

/* RANDOM NUMBER GENERATOR CONSTANTS */
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773 
#define IR 2836
#define MASK 123459876

/*---------------------------------------------------------------------------*\
 *  STRUCTURES                                                               *
\*---------------------------------------------------------------------------*/

typedef enum { FALSE, TRUE } bool;

typedef struct node {
  int left;
  int right;
  double distance;
} Node;

typedef struct data {
  unsigned long int index; /* index of the point */
  double x; /* cartesian coordinates of the point */
  double y;
  double z;
  double energy; /* interaction energy value */
  double cl_energy; /* total energy of the cluster of the point  */
  int cl_number; /* cluster index */
} Data;

typedef struct summary {
  double x; /* coordinates of the "center of energy" */
  double y;
  double z;
  double cl_energy;
  int nmembers;	/* how many elements the cluster has */
} Summary;

typedef struct atom {
  double x; /* coordinates of the atom */
  double y;
  double z;
  char chain; /* chain identifier of the residue */
  char resID[4]; /* residue name */
  unsigned int atomNum; /* atom number */
  unsigned int resNum; /* residue number */
  struct atom *next; /* pointer to the next atom in the linked list */
} Atom;

typedef struct restable {
  char resname[MAX_STR_LEN];
  bool printed; /* a boolean flag to mark already printed residues */
  struct restable *next;
} Restable;

/*---------------------------------------------------------------------------*\
 *  FUNCTION DECLARATIONS                                                    *
\*---------------------------------------------------------------------------*/

/* general purpose functions */
double ran(long *);
double euclidean_distance(const Data *, const Atom *);

/* input functions */
double *expand(unsigned int *, unsigned int, unsigned int);
void filter_autogrid_map(const char *, const double);
void set_parameters(int, char **);

/* clustering functions */
void cluster_grid(const char *, const double, const char *);
Data *create_data_matrix(const char *, int *);
double **create_distance_matrix(Data *, int);
double find_closest_pair(int, double **, int *, int *);
Node *sl_cluster(double **, int);
Node *average_cluster(double **, int);
int *cuttree(Node *, int, double, int *);
void rank_data(Data *, const int *, const int, const int);
Summary *create_summary(Data *, const int, const int);
int compare_nodes(const void *, const void *);
int compare_data(const void *, const void *);

/* output functions */
void print_clusters(Data *, const int, const char *);
void print_summary(Summary *, const int, const char *);
void create_pdb(const char *infilename);
void create_dx(const char *infilename);
void create_contacts(const char *infilename);
void print_residues(FILE *, Restable *);
void free_memory_restable(Restable *);
