/*********************************************************
 *  2-d Ising model, Metropolis Monte-Carlo simulation.  *
 *********************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>

# define error_reading -1
# define error_value -2

char filename[BUFSIZ + 1] = "";
int length = 20;
int step = 1000000;
double temperature = 2.0;

/*  the spin status  */
enum _spin {down = -1, up = 1};
typedef enum _spin spin;
/*  a spin matrix contains its content and size.  */
struct _matrix
{
    spin **content;
    int row;
    int col;
};
typedef struct _matrix matrix;
/*  below are for printing with colors using "escape code".  */
enum _fg_color {fg_Black = 30, fg_Red, fg_Green, fg_Yellow, fg_Blue, fg_Magenta, fg_Cyan, fg_White};
enum _bg_color {bg_Black = 40, bg_Red, bg_Green, bg_Yellow, bg_Blue, bg_Magenta, bg_Cyan, bg_White};
typedef enum _fg_color fg_color;
typedef enum _bg_color bg_color;

/*  allocate memory for a spin board.  */
matrix Init_matrix(int row, int col);

/*  release memory and set abandoned pointers to NULL.  */
void Delete_matrix(matrix *matp);

/*  read arguments from command line, generate initial board and show information about it.  */
matrix Read_data(int argc, char const *argv[]);

/*  initialize print with escape code  */
void Init_print_escape_code(void);

/*  print the matrix using escape code on terminal.  */
void Print_matrix(matrix mat);

/*  print the matrix in a formated format to a file.  */
void Print_matrix_simple(matrix mat, FILE *fp);

/*  calculates energy at one specific grid.  */
int Calc_one_energy(matrix mat, int posi, int posj);

/*  calculates energy of the whole board.  */
int Calc_energy(matrix mat);

/*  calculates total moment of magnet of a board.  */
int Calc_magnet(matrix mat);

/*  print a progress bar on terminal.  */
void Print_progress(double progress);

/*  do one operation on the board.  */
void Operate_matrix_one_time(matrix *matp);


int main(int argc, char const *argv[])
{
    int n = 0;
    matrix board = Read_data(argc, argv);
    FILE *ofp = NULL;
    int print_interval = (int)ceil(step / 100.0);
    int print_start = step % print_interval;
    time_t begin_time = time(NULL);

    printf("Running, please wait...\n");
    for (n = 0; n < step; ++ n)
    {
        if (! ((n - print_start) % print_interval))
            Print_progress((double)(n + 1) / step);
        Operate_matrix_one_time(& board);
    }
    printf("\n");
    printf("Time used for simulation: %3 s.\n", (int)(time(NULL) - begin_time));

    printf("Final   status: total magnetic moment is %5d m, total energy is %5d J.\n", 
        Calc_magnet(board), Calc_energy(board));

    Print_matrix(board);
    ofp = fopen("Metropolis_Monte_Carlo_result.txt", "wt");
    if (! ofp)
    {
        fprintf(stderr, "Cannot open \"Metropolis_Monte_Carlo_result.txt\"!\n");
        exit(error_reading);
    }
    Print_matrix_simple(board, ofp);
    printf("Minimum energy: total magnetic moment is ±%4d m, total energy is %5d J.\n", 
        board.row * board.col, - 2 * board.row * board.col);
    printf("A copy of the result has been set to \"Metropolis_Monte_Carlo_result.txt\".\n");

    fclose(ofp);
    ofp = NULL;
    Delete_matrix(& board);
    
    return 0;
}


matrix Init_matrix(int row, int col)
{
    int i = 0;
    matrix ret = {NULL, row, col};
    ret.content = row && col ? (spin **)malloc(row * sizeof(spin *)) : NULL;

    if (! ret.content)
        return ret;
    ret.content[0] = (spin *)malloc(row * col * sizeof(spin));
    for (i = 1; i < row; ++ i)
        ret.content[i] = ret.content[i - 1] + col;

    return ret;
}

void Delete_matrix(matrix *matp)
{
    int i = 0;

    if (matp->content[0])
        free(matp->content[0]);
    for (i = 0; i < matp->row; ++ i)
        matp->content[i] = NULL;
    if (matp->content)
        free(matp->content);
    matp->content = NULL;
    matp->row = 0;
    matp->col = 0;

    return;
}

matrix Read_data(int argc, const char *argv[])
{
    matrix ret = {NULL, 0, 0};
    int iarg = 1;
    FILE *ifp = NULL;
    int i = 0, j = 0;
    int half = 0;
    int tmp_value = 0;

    # ifdef __WIN32
    system("CHCP 65001 1> NUL");
    # endif
    iarg = 1;
    for (;;)
    {
        if (iarg >= argc)
            break;
        if (! strcmp(argv[iarg], "-h") || ! strcmp(argv[iarg], "--help"))
        {
            printf("This program performs a Metropolis Monte Carlo simulation of 2-dimensional Ising model.\n");
            printf("Usage: \n");
            printf("2d_Ising_Metropolis_Monte_Carlo_simulation.exe [OPTIONS]\n");
            printf("\n");
            printf("OPTIONS:\n");
            printf("-h, --help              Print this message and exit unnormally.\n");
            printf("-l, --length            Number of atoms on each side.\n");
            printf("-s, --step              Number of total steps.\n");
            printf("-t, --temperature       Temperature (in unit J/k_Boltzmann)\n");
            printf("-f, --filename          Name of file contains initial status, \n");
            printf("                        \"-\" means to read from stdin.\n");
            printf("\n");
            printf("The default argumets are \"-l 20 -s 1000000 -t 2.0\".\n");
            printf("And the default initial status is about the left half\n");
            printf("is up-spin, and about the right half is down-spin.\n");
            printf("\n");
            printf("Please pay attention that the Curie temperature \n");
            printf("of the analytic solution is 2/ln(1+sqrt(2)) = 2.269 J/k_B.\n");
            exit(EXIT_SUCCESS);
        }
        if (! strcmp(argv[iarg], "-l") || ! strcmp(argv[iarg], "--length"))
        {
            if (iarg == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            if (sscanf(argv[iarg], "%u", & length) != 1)
            {
                printf("Error! Cannot set value for \"length\" from \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            continue;
        }
        if (! strcmp(argv[iarg], "-s") || ! strcmp(argv[iarg], "--step"))
        {
            if (iarg == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            if (sscanf(argv[iarg], "%u", & step) != 1)
            {
                printf("Error! Cannot set value for \"step\" from \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            continue;
        }
        if (! strcmp(argv[iarg], "-t") || ! strcmp(argv[iarg], "--temperature"))
        {
            if (iarg == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            if (sscanf(argv[iarg], "%lf", & temperature) != 1)
            {
                printf("Error! Cannot set value for \"temperature\" from \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            continue;
        }
        if (! strcmp(argv[iarg], "-f") || ! strcmp(argv[iarg], "--filename"))
        {
            if (iarg == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[iarg]);
                exit(error_reading);
            }
            ++ iarg;
            strncpy(filename, argv[iarg], BUFSIZ + 1);
            ++ iarg;
            continue;
        }

        printf("Error! Unknown argument: \"%s\"\n", argv[iarg]);
        exit(error_reading);
    }

    if (length <= 1)
    {
        printf("Error! length must be at least 2, but it is %u.\n", length);
        exit(error_value);
    }
    if (step <= 0)
    {
        printf("Error! step must be at least 1, but it is %u.\n", step);
        exit(error_value);
    }  
    printf("length = %u, step = %u, temperature = %.1lf, filename = \"%s\"\n", 
        length, step, temperature, filename);

    ret = Init_matrix(length, length);
    if (strcmp(filename, ""))
    {
        if (strcmp(filename, "-"))
        {
            ifp = fopen(filename, "rt");
            if (! ifp)
            {
                printf("Error! Cannot read from \"%s\".\n", filename);
                exit(error_reading);
            }
            for (i = 0; i < ret.row; ++ i)
                for (j = 0; j < ret.col; ++ j)
                {
                    if (fscanf(ifp, "%d", & tmp_value) != 1)
                    {
                        fprintf(stderr, "Error setting initial status!\n");
                        exit(error_value);
                    }
                    ret.content[i][j] = tmp_value > 0 ? up : down;
                }
            fclose(ifp);
            ifp = NULL;
        }
        else
        {
            printf("Reading initial status from stdin, please input:\n");
            for (i = 0; i < ret.row; ++ i)
                for (j = 0; j < ret.col; ++ j)
                {
                    if (scanf("%d", & tmp_value) != 1)
                    {
                        fprintf(stderr, "Error setting initial status!\n");
                        exit(error_value);
                    }
                    ret.content[i][j] = tmp_value > 0 ? up : down;
                }
            while (getchar() != '\n')
                ;
        }
    }
    else
    {
        half = (ret.col - 1) / 2 + 1;
        if (ret.col % 2)
        {
            for (i = 0; i < half; ++ i)
            {
                for (j = 0; j < half; ++ j)
                    ret.content[i][j] = up;
                for (j = half; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
            for (i = half; i < ret.row; ++ i)
            {
                for (j = 0; j < half - 1; ++ j)
                    ret.content[i][j] = up;
                for (j = half - 1; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
        }
        else
        {
            for (i = 0; i < ret.row; ++ i)
            {
                for (j = 0; j < half; ++ j)
                    ret.content[i][j] = up;
                for (j = half; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
        }
    }

    printf("Initial status: total magnetic moment is %5d m, total energy is %5d J.\n", \
        Calc_magnet(ret), Calc_energy(ret));
    Init_print_escape_code();
    Print_matrix(ret);

    /*  set the seed of random here.  */
    srand((unsigned int)time(NULL));

    return ret;
}

inline void Init_print_escape_code(void)
{
    system("");

    return;
}

void Print_matrix(matrix mat)
{
    int i = 0, j = 0;

    for (i = 0; i < mat.row; ++ i)
    {
        for (j = 0; j < mat.col; ++ j)
        {
            if (j)
                printf(" ");
            if (mat.content[i][j] == up)
                printf("\033[%um↑\033[0m", fg_Red);
            else if (mat.content[i][j] == down)
                printf("\033[%um↓\033[0m", fg_Blue);
            else
                printf("\033[%um0\033[0m", fg_Green); /* should never happen */
        }
        printf("\n");
    }

    return;
}

void Print_matrix_simple(matrix mat, FILE *fp)
{
    int i = 0, j = 0;

    for (i = 0; i < mat.row; ++ i)
    {
        for (j = 0; j < mat.col; ++ j)
            fprintf(fp, " %2d" + ! j, mat.content[i][j]);
        fprintf(fp, "\n");
    }

    return;
}

int Calc_one_energy(matrix mat, int posi, int posj)
{
    int ene = 0;

    if (posi >= mat.row || posj >= mat.col)
    {
        printf("Error! pos (%u, %u) out of range (0, 0) to (%u, %u).", posi, posj, mat.row - 1, mat.col - 1);
        exit(error_value);
    }

    ene = 0;
    ene -= mat.content[posi][posj] * mat.content[posi > 0 ? posi - 1 : mat.row - 1][posj];
    ene -= mat.content[posi][posj] * mat.content[posi < mat.row - 1 ? posi + 1 : 0][posj];
    ene -= mat.content[posi][posj] * mat.content[posi][posj > 0 ? posj - 1 : mat.col - 1];
    ene -= mat.content[posi][posj] * mat.content[posi][posj < mat.col - 1 ? posj + 1 : 0];

    return ene;
}

int Calc_energy(matrix mat)
{
    int i = 0, j = 0;
    int ene = 0;

    ene = 0;
    for (i = 0; i < mat.row; ++ i)
        for (j = 0; j < mat.col; ++ j)
            ene += Calc_one_energy(mat, i, j);

    ene /= 2; /*  Each element has been calculated twice.  */

    return ene;
}

int Calc_magnet(matrix mat)
{
    int i = 0, j = 0;
    int mag = 0;

    mag = 0;
    for (i = 0; i < mat.row; ++ i)
        for (j = 0; j < mat.col; ++ j)
            mag += mat.content[i][j];

    return mag;
}

void Print_progress(double progress)
{
    const int nsymbol = 50;
    int ncomplete = (int)floor(progress * nsymbol);
    int i = 0;

    if (progress > 1.0 - (1.0 / nsymbol) / 2)
    {
        printf("\r[");
        for (i = 0; i < nsymbol; ++ i)
            printf("#");
        printf("]  100.00%%");
        return;
    }
    printf("\r[");
    for (i = 0; i < ncomplete; ++ i)
        printf("#");
    for (i = ncomplete; i < nsymbol; ++ i)
        printf("*");
    printf("]  %6.2lf%%", progress * 100);

    return;
}

void Operate_matrix_one_time(matrix *matp)
{
    /*  This algorithm is not very good as the probability of the last grid is a bit lower.  */
    int pos = (int)(matp->row * matp->col * ((double)rand() / (RAND_MAX + 1)));
    int posi = pos / matp->col, posj = pos % matp->col;
    int delta_ene = Calc_one_energy(* matp, posi, posj) * (-2);

    if (delta_ene <= 0 || exp(- delta_ene / temperature) >= (double)rand() / (RAND_MAX + 1))
        matp->content[posi][posj] = - matp->content[posi][posj];

    return;
}

