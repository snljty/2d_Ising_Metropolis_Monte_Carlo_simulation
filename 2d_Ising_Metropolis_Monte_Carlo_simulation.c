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

char filename[BUFSIZ + 1u] = "";
unsigned int length = 20u;
unsigned int step = 1E6;
double temperature = 2.0;

enum _bool {False, True};
typedef enum _bool Bool;
enum _spin {down = -1, up = 1};
typedef enum _spin spin;
struct _matrix
{
    spin **content;
    unsigned int row;
    unsigned int col;
};
typedef struct _matrix matrix;
enum _fg_color {fg_Black = 30, fg_Red, fg_Green, fg_Yellow, fg_Blue, fg_Magenta, fg_Cyan, fg_White};
enum _bg_color {bg_Black = 40, bg_Red, bg_Green, bg_Yellow, bg_Blue, bg_Magenta, bg_Cyan, bg_White};
typedef enum _fg_color fg_color;
typedef enum _bg_color bg_color;

matrix Init_matrix(unsigned int row, unsigned int col);

void Delete_matrix(matrix *matp);

matrix Read_data(int argc, const char *argv[]);

void Print_matrix(matrix mat);

void Print_matrix_simple(matrix mat, FILE *fp);

int Calc_one_energy(matrix mat, unsigned int posi, unsigned int posj);

int Calc_energy(matrix mat);

int Calc_magnet(matrix mat);

void Print_progress(double progress);

void Operate_matrix_one_time(matrix *matp);


int main(int argc, const char *argv[])
{
    unsigned int n = 0u;
    matrix board = Read_data(argc, argv);
    FILE *ofp = NULL;
    unsigned int print_interval = (unsigned int)ceil(step / 100.0);
    unsigned int print_start = step % print_interval;
    time_t begin_time = time(NULL);

    puts("Running, please wait...");
    for (n = 0u; n < step; ++ n)
    {
        if (! ((n - print_start) % print_interval))
            Print_progress((double)(n + 1u) / step);
        Operate_matrix_one_time(& board);
    }
    puts("");
    printf("Time used for simulation: %3u s.\n", (unsigned int)(time(NULL) - begin_time));

    printf("Final   status: total magnetic moment is %5d m, total energy is %5d J.\n", 
        Calc_magnet(board), Calc_energy(board));

    Print_matrix(board);
    ofp = fopen("Metropolis_Monte_Carlo_result.txt", "wt");
    if (! ofp)
    {
        puts("Cannot open \"Metropolis_Monte_Carlo_result.txt\"!");
        exit(error_reading);
    }
    Print_matrix_simple(board, ofp);
    printf("Minimum energy: total magnetic moment is ±%4d m, total energy is %5d J.\n", 
        board.row * board.col, -2 * board.row * board.col);
    puts("A copy of the result has been set to \"Metropolis_Monte_Carlo_result.txt\".");

    fclose(ofp);
    ofp = NULL;
    Delete_matrix(& board);
    
    return 0;
}


matrix Init_matrix(unsigned int row, unsigned int col)
{
    unsigned int i = 0u;
    matrix ret = {NULL, row, col};
    ret.content = row && col ? (spin **)malloc(row * sizeof(spin *)) : NULL;

    if (! ret.content)
        return ret;
    ret.content[0] = (spin *)malloc(row * col * sizeof(spin));
    for (i = 1u; i < row; ++ i)
        ret.content[i] = ret.content[i - 1] + col;

    return ret;
}

void Delete_matrix(matrix *matp)
{
    unsigned int i = 0u;

    if (matp->content[0])
        free(matp->content[0]);
    for (i = 0u; i < matp->row; ++ i)
        matp->content[i] = NULL;
    if (matp->content)
        free(matp->content);
    matp->content = NULL;
    matp->row = 0u;
    matp->col = 0u;

    return;
}

matrix Read_data(int argc, const char *argv[])
{
    matrix ret = {NULL, 0u, 0u};
    unsigned int t = 1u;
    FILE *ifp = NULL;
    unsigned int i = 0u, j = 0u;
    unsigned int half = 0u;
    int tmp_value = 0;

    t = 1u;
    while (True)
    {
        if (t >= argc)
            break;
        if (! strcmp(argv[t], "-h") || ! strcmp(argv[t], "--help"))
        {
            puts("This program performs a Metropolis Monte Carlo simulation of 2-dimensional Ising model.");
            puts("Usage: ");
            puts("2d_Ising_Metropolis_Monte_Carlo_simulation.exe [OPTIONS]");
            puts("");
            puts("OPTIONS:");
            puts("-h, --help              Print this message and exit unnormally.");
            puts("-l, --length            Number of atoms on each side.");
            puts("-s, --step              Number of total steps.");
            puts("-t, --temperature       Temperature (in unit J/k_Boltzmann)");
            puts("-f, --filename          Name of file contains initial status, ");
            puts("                        \"-\" means to read from stdin.");
            puts("");
            puts("The default argumets are \"-l 20 -s 1000000 -t 2.0\".");
            puts("And the default initial status is about the left half");
            puts("is up-spin, and about the right half is down-spin.");
            puts("");
            puts("Please pay attention that the Curie temperature ");
            puts("of the analytic solution is 2/ln(1+sqrt(2)) = 2.269 J/k_B.");
            exit(error_reading);
        }
        if (! strcmp(argv[t], "-l") || ! strcmp(argv[t], "--length"))
        {
            if (t == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            if (sscanf(argv[t], "%u", & length) != 1)
            {
                printf("Error! Cannot set value for \"length\" from \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            continue;
        }
        if (! strcmp(argv[t], "-s") || ! strcmp(argv[t], "--step"))
        {
            if (t == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            if (sscanf(argv[t], "%u", & step) != 1)
            {
                printf("Error! Cannot set value for \"step\" from \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            continue;
        }
        if (! strcmp(argv[t], "-t") || ! strcmp(argv[t], "--temperature"))
        {
            if (t == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            if (sscanf(argv[t], "%lf", & temperature) != 1)
            {
                printf("Error! Cannot set value for \"temperature\" from \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            continue;
        }
        if (! strcmp(argv[t], "-f") || ! strcmp(argv[t], "--filename"))
        {
            if (t == argc - 1)
            {
                printf("Error! Missing value for \"%s\".\n", argv[t]);
                exit(error_reading);
            }
            ++ t;
            strncpy(filename, argv[t], BUFSIZ + 1u);
            ++ t;
            continue;
        }

        printf("Error! Unknown argument: \"%s\"\n", argv[t]);
        exit(error_reading);
    }

    if (length <= 1)
    {
        printf("Error! length must be at least 2, but it is %u.\n", length);
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
            for (i = 0u; i < ret.row; ++ i)
                for (j = 0u; j < ret.col; ++ j)
                {
                    fscanf(ifp, "%d", & tmp_value);
                    ret.content[i][j] = tmp_value > 0 ? up : down;
                }
            fclose(ifp);
            ifp = NULL;
        }
        else
        {
            puts("Reading initial status from stdin, please input:");
            for (i = 0u; i < ret.row; ++ i)
                for (j = 0u; j < ret.col; ++ j)
                {
                    scanf("%d", & tmp_value);
                    ret.content[i][j] = tmp_value > 0 ? up : down;
                }
            while (getchar() != '\n')
                ;
        }
    }
    else
    {
        half = (ret.col - 1u) / 2u + 1u;
        if (ret.col % 2u)
        {
            for (i = 0u; i < half; ++ i)
            {
                for (j = 0u; j < half; ++ j)
                    ret.content[i][j] = up;
                for (j = half; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
            for (i = half; i < ret.row; ++ i)
            {
                for (j = 0u; j < half - 1u; ++ j)
                    ret.content[i][j] = up;
                for (j = half - 1u; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
        }
        else
        {
            for (i = 0u; i < ret.row; ++ i)
            {
                for (j = 0u; j < half; ++ j)
                    ret.content[i][j] = up;
                for (j = half; j < ret.col; ++ j)
                    ret.content[i][j] = down;
            }
        }
    }

    printf("Initial status: total magnetic moment is %5d m, total energy is %5d J.\n", 
        Calc_magnet(ret), Calc_energy(ret));
    Print_matrix(ret);

    /*  set the seed of random here.  */
    srand((unsigned int)time(NULL));

    return ret;
}

void Print_matrix(matrix mat)
{
    unsigned int i = 0u, j = 0u;

    for (i = 0u; i < mat.row; ++ i)
    {
        for (j = 0u; j < mat.col; ++ j)
        {
            if (j)
                printf(" ");
            if (mat.content[i][j] == up)
                printf("\033[%um↑\033[0m", fg_Red);
            else if (mat.content[i][j] == down)
                printf("\033[%um↓\033[0m", fg_Blue);
        }
        printf("\n");
    }

    return;
}

void Print_matrix_simple(matrix mat, FILE *fp)
{
    unsigned int i = 0u, j = 0u;

    for (i = 0u; i < mat.row; ++ i)
    {
        for (j = 0u; j < mat.col; ++ j)
            fprintf(fp, " %2d" + ! j, mat.content[i][j]);
        fprintf(fp, "\n");
    }

    return;
}

int Calc_one_energy(matrix mat, unsigned int posi, unsigned int posj)
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
    unsigned int i = 0u, j = 0u;
    int ene = 0;

    ene = 0;
    for (i = 0u; i < mat.row; ++ i)
        for (j = 0u; j < mat.col; ++ j)
            ene += Calc_one_energy(mat, i, j);

    ene /= 2; /*  Each element has been calculated twice.  */

    return ene;
}

int Calc_magnet(matrix mat)
{
    unsigned int i = 0u, j = 0u;
    int mag = 0;

    mag = 0;
    for (i = 0u; i < mat.row; ++ i)
        for (j = 0u; j < mat.col; ++ j)
            mag += mat.content[i][j];

    return mag;
}

void Print_progress(double progress)
{
    const unsigned int nsymbol = 50u;
    unsigned int ncomplete = (unsigned int)floor(progress * nsymbol);
    unsigned int i = 0u;

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
    unsigned int pos = (unsigned int)(matp->row * matp->col * ((double)rand() / (RAND_MAX + 1u)));
    unsigned int posi = pos / matp->col, posj = pos % matp->col;
    int delta_ene = Calc_one_energy(* matp, posi, posj) * (-2);

    if (delta_ene <= 0 || exp(- delta_ene / temperature) >= (double)rand() / (RAND_MAX + 1u))
        matp->content[posi][posj] = - matp->content[posi][posj];

    return;
}

