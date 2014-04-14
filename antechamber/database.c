/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    database                                                   *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# define MAXCHAR 256
# include "utility.c"
# define COLORTEXT "YES"
# define MAXCOMMAND 50
# define MAXLINE 50
# define debug 1
# define RANDOM_ARRAY_SIZE 1000
FILE *fpin;
FILE *fpout;
FILE *fperror;
FILE *fpprocess;
FILE *fptmp;

char ifilename[MAXCHAR];
char dfilename[MAXCHAR];
char tmp_file[MAXCHAR] = "TMP";
int i, j, k;

char rbegin[MAXCHAR];
char rend[MAXCHAR];
int rbegin_index = 0;
int rend_index = 0;

char separator[MAXCHAR];
int separator_index = 0;
int separator_type = 0;
int first_record_flag = 1;

char molname[MAXCHAR];
int molname_flag;

char def[MAXCHAR];
int def_index = 0;
long id = 1;

char field[MAXCHAR];
int molname_field_index = 0;
int molname_col = 0;
int molname_ld = 1;
int molname_fl = 0;
int molname_fr = 0;

char stop_line[MAXCHAR];
int stop_line_index = 0;

char omit_line[MAXLINE][MAXCHAR];
int omit_line_num = 0;

char add_line[MAXLINE][MAXCHAR];
int add_line_flag[MAXLINE];
int add_line_num = 0;

char command[MAXCOMMAND][MAXCHAR];
int ncom = 0;

long cmpd_num = 0;
long proc_begin = 0;
long proc_end = 0;
int proc_freq = 1;
int freq = 1;
double proc_prob = -1;
int process_flag;

int fileopen_flag;

long random_id;
long random_flag = 0;
long random_array_size = RANDOM_ARRAY_SIZE;
double *randomp;
double *randomr, randoms;
int ranr, ranp;
void qarn(double *r, double p[], int n);
void randomnum(int);

void calendar(int *year,
              int *month, int *day, int *hour, int *minute, int *second)
{
    time_t now = time(NULL);
    struct tm *timestruct = localtime(&now);
    *year = timestruct->tm_year + 1900;
    *month = timestruct->tm_mon + 1;
    *day = timestruct->tm_mday;
    *hour = timestruct->tm_hour;
    *minute = timestruct->tm_min;
    *second = timestruct->tm_sec;
}

int calendar_(int *year,
              int *month, int *day, int *hour, int *minute, int *second)
{
    calendar(year, month, day, hour, minute, second);
    return 0;
}

void qarn(double *r, double p[], int n)
{
    int i, m;
    double s, u, v;
    s = 65536.0;
    u = 2053.0;
    v = 13849.0;
    for (i = 0; i <= n - 1; i++) {
        *r = u * (*r) + v;
        m = (int) (*r / s);
        *r = *r - m * s;
        p[i] = *r / s;
    }
    return;
}

/* function to generat a set of random number between 0-1 */
void randomnum(int arraysize)
{
/* int *year, *month, *day, *hour, *minute, *second;
   calendar(year,month,day,hour,minute,second);
   printf("\n %5d%5d%5d", *hour,*minute,*second);
   randoms=hour+second;
 */
    time_t now = time(NULL);
    randoms = now + 100000 * randomp[random_flag];
    random_flag += 2;
    if (random_flag >= random_array_size)
        random_flag = 0;
    /*  printf("\n time is %5d",now);  */
    randomr = &randoms;
    qarn(randomr, randomp, arraysize);
}

double getrandomnum(void)
{
    if (random_id > random_array_size) {
        random_id = 0;
        randomnum(random_array_size);
        printf("\nregenerate a set of random numbers ...");
    }
    random_id++;
    return randomp[random_id];
}

void clean_addline(char *str)
{
    char newstr[MAXCHAR];
    int i;
    char tmpchar[2];
    strcpy(newstr, "");
    for (i = 0; i < strlen(str); i++) {
        if (i < strlen(str) - 6) {
            if (str[i] == 'T' && str[i + 1] == 'M'
                && str[i + 2] == 'P' && str[i + 3] == 'F'
                && str[i + 4] == 'I' && str[i + 5] == 'L'
                && str[i + 6] == 'E') {
                if (spaceline(tmp_file) == 1) {
                    fprintf(stdout, "ERROR: TMPFILE is null string, exit");
                    exit(0);
                }
                strcat(newstr, tmp_file);
                i += 7;
            }
            if (str[i] == 'M' && str[i + 1] == 'O'
                && str[i + 2] == 'L' && str[i + 3] == 'N'
                && str[i + 4] == 'A' && str[i + 5] == 'M'
                && str[i + 6] == 'E') {
                if (spaceline(molname) == 1) {
                    fprintf(stdout, "ERROR: MOLNAME is null string, exit");
                    exit(0);
                }
                strcat(newstr, molname);
                i += 7;
            }
        }
        tmpchar[0] = str[i];
        tmpchar[1] = '\0';
        strcat(newstr, tmpchar);
    }
    strcpy(str, newstr);
}
int runcommand(char *commandstr)
{
    int i;
    int status = 0;
    char realcommand[MAXCHAR] = "";
    char tmpchar[2];
    for (i = 0; i < strlen(commandstr); i++) {
        if (i < strlen(commandstr) - 6) {
            if (commandstr[i] == 'T' && commandstr[i + 1] == 'M'
                && commandstr[i + 2] == 'P' && commandstr[i + 3] == 'F'
                && commandstr[i + 4] == 'I' && commandstr[i + 5] == 'L'
                && commandstr[i + 6] == 'E') {
                if (spaceline(tmp_file) == 1) {
                    fprintf(stdout, "ERROR: TMPFILE is null string, exit");
                    exit(0);
                }
                strcat(realcommand, tmp_file);
                i += 7;
            }
            if (commandstr[i] == 'M' && commandstr[i + 1] == 'O'
                && commandstr[i + 2] == 'L' && commandstr[i + 3] == 'N'
                && commandstr[i + 4] == 'A' && commandstr[i + 5] == 'M'
                && commandstr[i + 6] == 'E') {
                if (spaceline(molname) == 1) {
                    fprintf(stdout, "ERROR: MOLNAME is null string, exit");
                    exit(0);
                }
                strcat(realcommand, molname);
                i += 7;
            }
        }
        tmpchar[0] = commandstr[i];
        tmpchar[1] = '\0';
        strcat(realcommand, tmpchar);
    }
    status = system(realcommand);
    if (status != 0) {
        fprintf(fperror, "\nNo: %-8ld, Molecualr Name = %s\n", cmpd_num,
                molname);
        fprintf(fperror,
                "Error: cannot run \"%s\" in runcommand() of database.c properly, exit\n",
                realcommand);
        fflush(fperror);
        return 1;
    }
    return 0;
}

void rdef(char *filename)
{
    int i, j, k;
    int num_command = 0;
    int num_omit = 0;
    int num_add = 0;
    int num_quota;
    int index;
    int num;
    char line[MAXCHAR];
    FILE *fpin;

    if ((fpin = fopen(filename, "r")) == NULL) {
        printf("\n Cannot open defintion file: %s , exit", filename);
        exit(1);
    }
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        if (strncmp("SEPARATOR", line, 9) == 0) {
            index = 0;
            num_quota = 0;
            k = 0;
            separator_index = 1;
            for (j = 9; j < strlen(line); j++) {
                if (line[j] == '\"')
                    num_quota++;
                if (index == 0 && num_quota == 1) {
                    index = 1;
                    continue;
                }
                if (num_quota == 2) {
                    separator[k] = '\0';
                    break;
                }
                if (index == 1) {
                    separator[k] = line[j];
                    k++;
                }
            }
        }
        if (strncmp("SEPA_TYPE", line, 9) == 0) {
            sscanf(&line[9], "%d", &separator_type);
            continue;
        }
        if (strncmp("RECORD_BEGIN", line, 12) == 0) {
            index = 0;
            k = 0;
            num_quota = 0;
            rbegin_index = 1;
            rbegin[0] = '\0';
            for (j = 12; j < strlen(line); j++) {
                if (line[j] == '\"')
                    num_quota++;
                if (index == 0 && num_quota == 1) {
                    index = 1;
                    continue;
                }
                if (num_quota == 2) {
                    rbegin[k] = '\0';
                    break;
                }
                if (index == 1) {
                    rbegin[k] = line[j];
                    k++;
                }
            }
        }

        if (strncmp("RECORD_END", line, 10) == 0) {
            index = 0;
            k = 0;
            num_quota = 0;
            rend_index = 1;
            rend[0] = '\0';
            for (j = 10; j < strlen(line); j++) {
                if (line[j] == '\"')
                    num_quota++;
                if (index == 0 && num_quota == 1) {
                    index = 1;
                    continue;
                }
                if (num_quota == 2) {
                    rend[k] = '\0';
                    break;
                }
                if (index == 1) {
                    rend[k] = line[j];
                    k++;
                }
            }
        }
        if (strncmp("MOLNAME_DEF", line, 11) == 0) {
            sscanf(&line[11], "%s", def);
            def_index = 1;
            continue;
        }
        if (strncmp("MOLNAME_ID", line, 10) == 0) {
            sscanf(&line[10], "%ld", &id);
            continue;
        }
        if (strncmp("MOLNAME_FIELD", line, 13) == 0) {
            index = 0;
            k = 0;
            num_quota = 0;
            molname_field_index = 1;
            field[0] = '\0';
            for (j = 13; j < strlen(line); j++) {
                if (line[j] == '\"')
                    num_quota++;
                if (index == 0 && num_quota == 1) {
                    index = 1;
                    continue;
                }
                if (num_quota == 2) {
                    field[k] = '\0';
                    break;
                }
                if (index == 1) {
                    field[k] = line[j];
                    k++;
                }
            }
        }
        if (strncmp("MOLNAME_LD", line, 10) == 0) {
            sscanf(&line[10], "%d", &molname_ld);
            continue;
        }
        if (strncmp("MOLNAME_FL", line, 10) == 0) {
            sscanf(&line[10], "%d", &molname_fl);
            continue;
        }
        if (strncmp("MOLNAME_FR", line, 10) == 0) {
            sscanf(&line[10], "%d", &molname_fr);
            continue;
        }
        if (strncmp("TMPFILE", line, 7) == 0) {
            sscanf(&line[7], "%s", tmp_file);
            continue;
        }
        if (strncmp("PROC_BEGIN", line, 10) == 0) {
            sscanf(&line[10], "%ld", &proc_begin);
            continue;
        }
        if (strncmp("PROC_END", line, 8) == 0) {
            sscanf(&line[8], "%ld", &proc_end);
            continue;
        }
        if (strncmp("PROC_FREQ", line, 9) == 0) {
            sscanf(&line[9], "%d", &proc_freq);
            continue;
        }
        if (strncmp("PROC_PROB", line, 9) == 0) {
            sscanf(&line[9], "%lf", &proc_prob);
            continue;
        }
        if (strncmp("IMMEDIATE_STOP", line, 14) == 0) {
            index = 0;
            k = 0;
            num_quota = 0;
            stop_line_index = 1;
            stop_line[0] = '\0';
            for (j = 14; j < strlen(line); j++) {
                if (line[j] == '\"')
                    num_quota++;
                if (index == 0 && num_quota == 1) {
                    index = 1;
                    continue;
                }
                if (num_quota == 2) {
                    stop_line[k] = '\0';
                    break;
                }
                if (index == 1) {
                    stop_line[k] = line[j];
                    k++;
                }
            }
        }
        if (strncmp("OMIT_LINE", line, 9) == 0) {
            index = 0;
            k = 0;
            num = 0;
            num_quota = 0;
            for (j = 9; j < strlen(line); j++)
                if (line[j] == '\"')
                    num_quota++;
            for (j = 9; j < strlen(line); j++) {
                if (line[j] == '\"') {
                    num++;
                    if (index == 0) {
                        index = 1;
                        continue;
                    }
                    if (num == num_quota && index == 1) {
                        index = 2;
                        omit_line[num_omit][k] = '\0';
                    }
                }
                if (index == 1) {
                    omit_line[num_omit][k] = line[j];
                    k++;
                }
            }
            num_omit++;
            if (num_omit > MAXLINE) {
                printf
                    ("\nError: the number of omitted line exceed MAXLINE, increase MAXLINE, exit");
                exit(1);
            }
            continue;
        }
        if (strncmp("ADD_LINE", line, 8) == 0) {
            index = 0;
            k = 0;
            num = 0;
            num_quota = 0;
            sscanf(&line[8], "%d", &add_line_flag[num_add]);
            for (i = 8; i < strlen(line); i++)
                if (line[i] == '\"')
                    num_quota++;
            for (i = 8; i < strlen(line); i++)
                if (line[i] == '\"') {
                    j = i;
                    break;
                }
            for (j = 8; j < strlen(line); j++) {
                if (line[j] == '\"') {
                    num++;
                    if (index == 0) {
                        index = 1;
                        continue;
                    }
                    if (num == num_quota && index == 1) {
                        index = 2;
                        add_line[num_add][k] = '\0';
                    }
                }
                if (index == 1) {
                    add_line[num_add][k] = line[j];
                    k++;
                }
            }
            num_add++;
            if (num_add > MAXLINE) {
                printf
                    ("\nError: the number of adding line exceed MAXLINE, increase MAXLINE, exit");
                exit(1);
            }
            continue;
        }
        if (strncmp("COMMAND", line, 7) == 0) {
            index = 0;
            k = 0;
            num = 0;
            num_quota = 0;
            for (j = 7; j < strlen(line); j++)
                if (line[j] == '\"')
                    num_quota++;
            for (j = 7; j < strlen(line); j++) {
                if (line[j] == '\"') {
                    num++;
                    if (index == 0) {
                        index = 1;
                        continue;
                    }
                    if (num == num_quota && index == 1) {
                        index = 2;
                        command[num_command][k] = '\0';
                    }
                }
                if (index == 1) {
                    command[num_command][k] = line[j];
                    k++;
                }
            }
            num_command++;
            if (num_command > MAXCOMMAND) {
                printf
                    ("\nError: the number of commands exceed MAXCOMMAND, increase MAXCOMMAND, exit");
                exit(1);
            }
            continue;
        }
    }
    ncom = num_command;
    omit_line_num = num_omit;
    add_line_num = num_add;
    fclose(fpin);

    if (rbegin_index == 0 || rend_index == 0) {
        rbegin_index = 0;
        rend_index = 0;
    }
    if (rbegin_index == 1 && rend_index == 1) {
        separator_index = 0;
        molname_fr = 0;
    }
    if (separator_index == 0 && rbegin_index == 0 && rend_index == 0) {
        printf("\nNo Separation Field, exit\n");
        exit(1);
    }
    if (def_index == 1) {
        molname_field_index = 0;
        molname_fl = 0;
    }
    if (molname_field_index == 1) {
        molname_flag = 0;
        if (rbegin_index == 1)
            if (strcmp(rbegin, field) == 0)
                molname_flag = 1;
        if (rend_index == 1)
            if (strcmp(rend, field) == 0)
                molname_flag = 2;
        if (separator_index == 1)
            if (strcmp(separator, field) == 0)
                molname_flag = 3;
        molname_fl = 0;
        molname_fr = 0;
    }
    if (molname_fl == 0)
        molname_fr = 0;
    if (def_index == 0 && molname_field_index == 0 && molname_fl == 0) {
        printf("\nNo Molecular Name, exit\n");
        exit(1);
    }
    if (molname_fl != 0)
        if (molname_ld < 1) {
            printf
                ("\nMolname_ld must >=1 when the first line of each record is taken as the reference line, exit");
            exit(1);
        }
    if (proc_begin <= 0)
        proc_begin = 1;
    if (proc_end <= 0)
        proc_end = 0;
    if (proc_freq < 1)
        proc_freq = 1;
    if (proc_freq != 1)
        proc_prob = -1;
    if (proc_freq == 1 && (proc_prob <= 0 || proc_prob >= 1))
        proc_prob = -1;
    if (debug) {
        if (rbegin_index == 1)
            printf("\nRecord begin string = %s", rbegin);
        if (rend_index == 1)
            printf("\nRecord end string = %s", rend);
        if (separator_index == 1) {
            printf("\nRecord separator = %s", separator);
            printf("\nRecord separator type = %d", separator_type);
        }
        if (def_index == 1) {
            printf("\nMolname defined = %s", def);
            printf("\nMolname defined starting Id = %ld", id);
        }
        if (molname_field_index == 1)
            printf("\nMolname Field = %s", field);
        if (molname_fl == 1)
            printf
                ("\nFirst line of each record is taken as the reference line for record name reading");
        if (molname_fr == 1)
            printf
                ("\nFor the first record, the first line is taken as the reference line for record name reading");
        if (molname_fr == 2)
            printf
                ("\nFor the first record, the second line is taken as the reference line for record name reading");
        if (molname_fr == 3)
            printf
                ("\nFor the first record, the third line is taken as the reference line for record name reading");
        if (molname_fr > 3)
            printf
                ("\nFor the first record, the %dth line is taken as the reference line for record name reading",
                 molname_fr);
        printf
            ("\nThe distance of the reference line to the line where to read record names: %d",
             molname_ld);
        printf("\nRecord begin: %ld", proc_begin);
        printf("\nRecord end: %ld", proc_end);
        printf("\nRecord read frequency: %d", proc_freq);
        printf("\nRecord read probability: %lf", proc_prob);
        if (stop_line_index == 1)
            printf("\nImmediate stop line %s\n", stop_line);
        printf("\nTmp_file = %s", tmp_file);
        printf("\nNo of omitted line = %d", omit_line_num);
        for (i = 0; i < omit_line_num; i++)
            printf("\nOmitted line (%d) =  %s", i + 1, omit_line[i]);
        printf("\nNo of added line = %d", add_line_num);
        for (i = 0; i < add_line_num; i++)
            printf("\nAdd line (%d) =  %s (%d)", i + 1, add_line[i],
                   add_line_flag[i]);
        printf("\nNo of commands = %d", ncom);
        for (i = 0; i < ncom; i++)
            printf("\nCommand (%d) =  %s", i + 1, command[i]);
        printf("\n");
    }
}
void rdatabase(char *filename)
{
    int i;
    int status;
    int flag;
    int separator_line_flag;
    int rbegin_line_flag;
    int rend_line_flag;
    int tmpint;
    long pos;
    char line[MAXCHAR] = "";
    char linebak5[MAXCHAR];
    char linebak4[MAXCHAR];
    char linebak3[MAXCHAR];
    char linebak2[MAXCHAR];
    char linebak1[MAXCHAR];
    char tmpchar[MAXCHAR];

    fileopen_flag = 0;
    if ((fpin = fopen(filename, "r")) == NULL) {
        printf("\n Cannot open database file: %s ,exit", filename);
        exit(1);
    }
    if (separator_index == 1) {
        if (separator_type == 0)
            cmpd_num = 0;
        if (separator_type == 1) {
            cmpd_num = 1;
            id++;
        }
        process_flag = 1;
        if (proc_prob != -1)
            if (getrandomnum() > proc_prob)
                process_flag = 0;
        if (process_flag == 1)
            if (fileopen_flag == 1) {
                fclose(fptmp);
                fileopen_flag = 0;
            }
        if ((fptmp = fopen(tmp_file, "w")) == NULL) {
            printf("\n Cannot open tmp_file to write: file %s ,exit",
                   tmp_file);
            exit(1);
        }
        fileopen_flag = 1;
/* handle molname for the first record*/
        if (def_index == 1) {
            strcpy(molname, def);
            sprintf(tmpchar, "%ld", id - 1);
            /* newitoa(id - 1, tmpchar); */
            strcat(molname, tmpchar);
        }
/*add lines to the benginning of the file*/
        for (i = 0; i < add_line_num; i++)
            if (add_line_flag[i] == 1)
                fprintf(fptmp, "%s\n", add_line[i]);
    }
/* read in molname for the first record if molname_fr is applied*/
    if (molname_fr != 0 && process_flag == 1) {
        tmpint = 0;
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            if (tmpint == molname_fr - 1) {
                sscanf(&line[molname_col], "%s", molname);
                break;
            }
            tmpint++;
        }
        rewind(fpin);
    }
    for (;;) {
        strcpy(linebak5, linebak4);
        strcpy(linebak4, linebak3);
        strcpy(linebak3, linebak2);
        strcpy(linebak2, linebak1);
        strcpy(linebak1, line);
        if (fgets(line, MAXCHAR, fpin) == NULL) {
//for the last record
            if (process_flag == 1 && separator_index == 1
                && separator_type == 0) {
/*add lines to the end of the file*/
                for (i = 0; i < add_line_num; i++)
                    if (add_line_flag[i] == -1) {
                        clean_addline(add_line[i]);
                        fprintf(fptmp, "%s\n", add_line[i]);
                    }
                if (fileopen_flag == 1) {
                    fclose(fptmp);
                    fileopen_flag = 0;
                }
                if (def_index == 1) {
                    strcpy(molname, def);
                    sprintf(tmpchar, "%ld", id - 1);
                      /* newitoa(id - 1, tmpchar); */
                    strcat(molname, tmpchar);
                }

                printf("\nProcessing: %-8ld id=%7ld, comp. name = %s\n",
                       cmpd_num, id - 1, molname);
                fflush(stdout);
                fprintf(fpprocess,
                        "\nProcessing: %-8ld id=%7ld, comp. name = %s\n",
                        cmpd_num, id - 1, molname);
                fflush(fpprocess);
/*run commands */
                for (i = 0; i < ncom; i++) {
                    status = runcommand(command[i]);
                    if (status == 1)
                        break;
/*if error happens, stop running*/
                }
            }
            break;
        }
/* check stop line */
        if (stop_line_index == 1)
            if (strstr(line, stop_line) != 0)
                break;
/* check omitted lines */
        flag = 0;
        for (i = 0; i < omit_line_num; i++)
            if (strstr(line, omit_line[i]) != 0) {
                flag = 1;
                break;
            }
        if (flag == 1)
            continue;

        rbegin_line_flag = 0;
        rend_line_flag = 0;
        separator_line_flag = 0;

/*get record name if molname_field_index == 1 for most of cases*/
        if (molname_field_index == 1 && strstr(line, field) != 0)
            if (molname_flag <= 1 || (molname_flag == 2 && molname_ld <= 0)
                || (molname_flag == 3 && molname_ld <= 0)) {
                if (molname_ld == -5)
                    sscanf(&linebak5[molname_col], "%s", molname);
                if (molname_ld == -4)
                    sscanf(&linebak4[molname_col], "%s", molname);
                if (molname_ld == -3)
                    sscanf(&linebak3[molname_col], "%s", molname);
                if (molname_ld == -2)
                    sscanf(&linebak2[molname_col], "%s", molname);
                if (molname_ld == -1)
                    sscanf(&linebak1[molname_col], "%s", molname);
                if (molname_ld == 0)
                    sscanf(&line[molname_col], "%s", molname);
                if (molname_ld > 0) {
                    strcpy(linebak5, linebak4);
                    strcpy(linebak4, linebak3);
                    strcpy(linebak3, linebak2);
                    strcpy(linebak2, linebak1);
                    strcpy(linebak1, line);
                }
                if (molname_ld == 1) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 2) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 3) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 4) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 5) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
            }
/*if SEPARATOR LINE is applied*/
        if (separator_index == 1)
            if (strstr(line, separator) != 0) {
                separator_line_flag = 1;
                if (separator_type == 0 && first_record_flag == 1) {
                    first_record_flag = 0;
                    if (process_flag == 1)
                        process_flag = 0;
                }
                if (process_flag == 1) {
                    if (separator_type == 1)
                        fprintf(fptmp, "%s", line);
/*add lines to the end of the file*/
                    for (i = 0; i < add_line_num; i++)
                        if (add_line_flag[i] == -1) {
                            clean_addline(add_line[i]);
                            fprintf(fptmp, "%s\n", add_line[i]);
                        }
/*get record name */
                    if (fileopen_flag == 1) {
                        fclose(fptmp);
                        fileopen_flag = 0;
                    }
                    if (def_index == 1) {
                        strcpy(molname, def);
                        sprintf(tmpchar, "%ld", id - 1);
                        /* newitoa(id - 1, tmpchar); */
                        strcat(molname, tmpchar);
                    }
                    printf("\nProcessing: %-8ld id=%7ld, comp. name = %s",
                           cmpd_num, id - 1, molname);
                    fflush(stdout);
                    fprintf(fpprocess,
                            "\nProcessing: %-8ld id=%7ld, comp. name = %s",
                            cmpd_num, id - 1, molname);
                    fflush(fpprocess);
/*run commands */
/*if error happens, stop running*/
                    for (i = 0; i < ncom; i++) {
                        status = runcommand(command[i]);
                        if (status == 1)
                            break;
                    }
                }
/* prepare next record*/
                cmpd_num++;
                id++;
                if (cmpd_num > proc_end && proc_end > 0)
                    break;
                if (proc_prob == -1) {
                    if (proc_freq == freq) {
                        freq = 1;
                        process_flag = 1;
                    } else {
                        freq++;
                        process_flag = 0;
                    }
                }
                if (proc_prob != -1) {
                    if (getrandomnum() <= proc_prob) {
                        process_flag = 1;
                    } else {
                        process_flag = 0;
                    }
                }
                if (cmpd_num < proc_begin)
                    process_flag = 0;
                if (process_flag == 1) {
                    if (fileopen_flag == 1) {
                        fclose(fptmp);
                        fileopen_flag = 0;
                    }
                    if ((fptmp = fopen(tmp_file, "w")) == NULL) {
                        printf
                            ("\n Cannot open a tmp_file to write: %s ,exit",
                             tmp_file);
                        exit(1);
                    }
                    fileopen_flag = 1;
                    if (molname_fl == 1 && cmpd_num != 1) {
                        strcpy(linebak5, linebak4);
                        strcpy(linebak4, linebak3);
                        strcpy(linebak3, linebak2);
                        strcpy(linebak2, linebak1);
                        strcpy(linebak1, line);
                        if (molname_ld == 1) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 2) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 3) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 4) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 5) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                    }
/*add lines to the benginning of the file*/
                    for (i = 0; i < add_line_num; i++)
                        if (add_line_flag[i] == 1)
                            fprintf(fptmp, "%s\n", add_line[i]);
                }
            }
/*if RECORD_BEGIN and RECORD_END are applied*/
        if (rbegin_index == 1 && rend_index == 1)
            if (strstr(line, rbegin) != 0) {
                rbegin_line_flag = 1;
                cmpd_num++;
                id++;
                if (cmpd_num > proc_end && proc_end > 0)
                    break;
                if (proc_prob == -1){
                    if (proc_freq == freq) {
                        freq = 1;
                        process_flag = 1;
                    } else {
                        freq++;
                        process_flag = 0;
                    }
                }
                if (proc_prob != -1){
                    if (getrandomnum() <= proc_prob){
                        process_flag = 1;
                    } else {
                        process_flag = 0;
                    }
                }
                if (cmpd_num < proc_begin)
                    process_flag = 0;
                if (process_flag == 1) {
                    if (fileopen_flag == 1) {
                        fclose(fptmp);
                        fileopen_flag = 0;
                    }
                    if ((fptmp = fopen(tmp_file, "w")) == NULL) {
                        printf
                            ("\n Cannot open tmp_file to write: %s ,exit",
                             tmp_file);
                        exit(1);
                    }
                    fileopen_flag = 1;
                    if (def_index == 1) {
                        strcpy(molname, def);
                        sprintf(tmpchar, "%ld", id - 1);
                        /* newitoa(id - 1, tmpchar); */
                        strcat(molname, tmpchar);
                    }
                    if (molname_fl == 1) {
                        if (molname_ld > 0) {
                            strcpy(linebak5, linebak4);
                            strcpy(linebak4, linebak3);
                            strcpy(linebak3, linebak2);
                            strcpy(linebak2, linebak1);
                            strcpy(linebak1, line);
                        }
                        if (molname_ld == 1) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 2) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 3) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 4) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                        if (molname_ld == 5) {
                            pos = ftell(fpin);
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            if (fgets(line, MAXCHAR, fpin) == NULL)
                                break;
                            sscanf(&line[molname_col], "%s", molname);
                            fseek(fpin, pos, 0);
                            strcpy(line, linebak1);
                        }
                    }
/*add lines to the benginning of the file*/
                    for (i = 0; i < add_line_num; i++)
                        if (add_line_flag[i] == 1)
                            fprintf(fptmp, "%s\n", add_line[i]);
                }
            }
        if (rbegin_index == 1 && rend_index == 1)
            if (strstr(line, rend) != 0) {
                rend_line_flag = 1;
                if (process_flag == 1) {
                    fprintf(fptmp, "%s", line);
/*add lines to the end of the file*/
                    for (i = 0; i < add_line_num; i++)
                        if (add_line_flag[i] == -1) {
                            clean_addline(add_line[i]);
                            fprintf(fptmp, "%s\n", add_line[i]);
                        }
/*get record name */
                    if (fileopen_flag == 1) {
                        fclose(fptmp);
                        fileopen_flag = 0;
                    }
                    printf("\nProcessing: %-8ld id=%7ld, comp. name = %s",
                           cmpd_num, id - 1, molname);
                    fflush(stdout);
                    fprintf(fpprocess,
                            "\nProcessing: %-8ld id=%7ld, comp. name = %s",
                            cmpd_num, id - 1, molname);
                    fflush(fpprocess);
/*run commands */
                    for (i = 0; i < ncom; i++) {
                        status = runcommand(command[i]);
                        if (status == 1)
                            break;
/*if error happens, stop running*/
                    }
                }
            }

/*get record name if molname_field_index == 1 for some cases*/
        if (molname_field_index == 1 && strstr(line, field) != 0)
            if ((molname_flag == 2 && molname_ld > 0) || (molname_flag == 3
                && molname_ld > 0)) {
                strcpy(linebak5, linebak4);
                strcpy(linebak4, linebak3);
                strcpy(linebak3, linebak2);
                strcpy(linebak2, linebak1);
                strcpy(linebak1, line);
                if (molname_ld == 1) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 2) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 3) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 4) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
                if (molname_ld == 5) {
                    pos = ftell(fpin);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(&line[molname_col], "%s", molname);
                    fseek(fpin, pos, 0);
                    strcpy(line, linebak1);
                }
            }
        if (process_flag == 1){
            if (rend_line_flag == 1 || (separator_line_flag == 1
                && separator_type == 1)){
                continue;
            } else {
                fprintf(fptmp, "%s", line);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int i;


    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        if (argc == 2
            && (strcmp(argv[1], "-h") == 0
                || strcmp(argv[1], "-H") == 0)) {
            printf("[31mUsage: database -i[0m database file name\n"
                   "[31m                -d[0m definition file name\n");
            exit(1);
        }
        if (argc != 5) {
            printf("[31mUsage: database -i[0m database file name\n"
                   "[31m                -d[0m definition file name\n");
            exit(1);
        }
    } else {
        if (argc == 2
            && (strcmp(argv[1], "-h") == 0
                || strcmp(argv[1], "-H") == 0)) {
            printf("Usage: database -i database file name\n"
                   "                -d definition file name\n");
            exit(1);
        }
        if (argc != 5) {
            printf("Usage: database -i database file name\n"
                   "                -d definition file name\n");
            exit(1);
        }

    }

    for (i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0)
            strcpy(ifilename, argv[i + 1]);
        if (strcmp(argv[i], "-d") == 0)
            strcpy(dfilename, argv[i + 1]);
    }
    if ((fperror = fopen("error.log", "w")) == NULL) {
        printf("\n Cannot open file \"error.log\", exit\n");
        exit(1);
    }
    if ((fpprocess = fopen("process.log", "w")) == NULL) {
        printf("\n Cannot open fiel \"process.log\", exit\n");
        exit(1);
    }
    rdef(dfilename);
    if (proc_prob != -1) {
        randomp = (double *) malloc(sizeof(double) * random_array_size);
        if (randomp == NULL) {
            fprintf(stdout, "memory allocation error for *randomp\n");
            exit(1);
        }
        random_id = 0;
        randomp[0] = 0;
        randomnum(random_array_size);
    }
    rdatabase(ifilename);
    fclose(fperror);
    fclose(fpprocess);
    fclose(fpin);
    if (fileopen_flag == 1)
        fclose(fptmp);
    printf("\n");
    return (0);
}
