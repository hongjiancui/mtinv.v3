#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

#define RESULT_FOLDER_FLAG "-"

#define MAX_PATH 256
#define MAX_CMD 512

typedef struct FilenameList
{
    char *filename;
    struct FilenameList *next;
} nameList;

typedef struct PVR_result
{
    double maxvalue;
    char *filename;
    int Index;
} PVR_result;

/***
 * 函数名称：is_contained
 * 函数功能：检测c 是否包含于 s
 ***/
int is_contained(char *s, char *c)
{
    int i = 0, j = 0, flag = 0;
    while (i < strlen(s) && j < strlen(c))
    {
        if (s[i] == c[j])
        {
            i++;
            j++;
        }
        else
        {
            i = i - j + 1;
            j = 0;
        }
        if (j == strlen(c))
        {
            flag = 1;
            break;
        }
    }
    return flag;
}

/***
 * 函数名称：getFiles
 * 函数功能：获取特定目录下包含对应名称的文件名链表
 ***/
nameList *getFiles(char *path, char *perfix)
{
    DIR *dir;
    struct dirent *dir_ent;
    nameList *fnameList = NULL;
    if ((dir = opendir(path)) == NULL)
    {
        printf("路径打开失败\n");
        return NULL;
    }
    else
    {
        printf("Open dir %s\n", path);
        while (dir_ent = readdir(dir)) //read一次指针自动后移
        {
            if (is_contained(dir_ent->d_name, perfix))
            {
                nameList *p = (nameList *)malloc(sizeof(nameList));
                p->filename = dir_ent->d_name;
                p->next = fnameList;
                fnameList = p;
            }
        }
    }
    //closedir(dir);
    return fnameList;
}

/***
 * 函数名称：getNum
 * 函数功能：返回字符串中的数值
 * ***/
double getNum(char *str)
{
    double output;
    char temp[10];
    char *pend;
    int temp_count = 0; //content of temp
    for (int i = 0; i < strlen(str); i++)
    {
        if ((str[i] > 47 && str[i] < 58) || str[i] == 46)
        {
            if (temp_count == 10)
            {
                return 0; //溢出
            }
            temp[temp_count] = str[i];
            temp_count++;
        }
    }
    output = strtod(temp, &pend);
    return output;
}

/***
 * 函数名称：getLine
 * 函数功能：返回某文件中特定一整行字符串
 * ***/
char *getLine(char *str, int maxcol, char *path, int line)
{
    FILE *fp;
    fp = fopen(path, "r");
    if (fp == NULL)
    {
        perror("打开文件时发生错误");
        return ("error");
    }
    printf("Checking file %s\n", path);
    int count = 1;
    while (count <= line)
    {
        fgets(str, maxcol, fp);
        count++;
    }
    fclose(fp);
    return str; //函数返回字符串的四种方法
}

/***
 * 函数名称：getMaxPVR
 * 函数功能：找出path下计算文件中PVR值最高的一个，并返回PVR值
 * ***/
PVR_result getMaxPVR(char *path)
{
    double maxpvr_value = 0;
    double temp = 0;
    char *maxpvr_filename;
    char filepath[MAX_PATH];
    char line[MAX_PATH];
    char perfix[] = "email_";
    nameList *list;
    PVR_result result;
    result.maxvalue = 0;
    list = getFiles(path, perfix);

    if (NULL == list)
    {
        result.Index = 0;
        return result;
    }

    while (list != NULL)
    {
        memset(filepath, 0, sizeof(char) * MAX_PATH);
        sprintf(filepath, "%s/%s", path, list->filename);
        getLine(line, MAX_PATH, filepath, 12);
        temp = getNum(line);
        printf("pvr value is %1.4lf\n", temp);
        if (temp > result.maxvalue)
        {
            //maxpvr_value = temp;
            //maxpvr_filename = list->filename;
            result.maxvalue = temp;
            result.filename = list->filename;
        }
        list = list->next;
    }
    printf("Dir %s. MaxPVR:%1.4lf%%，FileName:%s\n", path, result.maxvalue, result.filename);

    nameList *p1 = list;
    nameList *p2;
    while (p1 != NULL)
    {
        p2 = p1;
        p1 = p1->next;
        free(p2);
    }

    result.Index = 1;
    return result;
}

/***
 * 函数名称：BestResult
 * 函数功能：核心函数：计算路径下所有文件夹中的计算文件中最大PVR值，并以字符串给出最大值和所在文件夹
 * ***/
unsigned int BestResult(char *path, char *pBest)
{
    char perfix_folders[] = RESULT_FOLDER_FLAG;
    char *maxpvr_filename;
    char *maxpvr_foldername;
    nameList *folderlist;
    double maxpvr = 0;
    PVR_result temp;
    folderlist = getFiles(path, perfix_folders);

    if (NULL == folderlist)
    {
        return 0;
    }

    for (nameList *pointer = folderlist; pointer != NULL; pointer = pointer->next)
    {
        char filepath[MAX_PATH];
        memset(filepath, 0, sizeof(char) * MAX_PATH);
        sprintf(filepath, "%s/%s", path, pointer->filename);
        temp = getMaxPVR(filepath);

        if (0 == temp.Index)
        {
            continue;
        }

        if (temp.maxvalue > maxpvr)
        {
            maxpvr = temp.maxvalue;
            maxpvr_filename = temp.filename;
            maxpvr_foldername = pointer->filename;
        }
    }

    printf("%s\n", maxpvr_foldername);

    memset(pBest, 0, sizeof(char) * MAX_PATH);
    memcpy(pBest, maxpvr_foldername, strlen(maxpvr_foldername));

    printf("最大PVR值为:%1.4lf%%，文件名：%s, 所属文件夹：%s\n", maxpvr, maxpvr_filename, maxpvr_foldername);

    nameList *p1 = folderlist;
    nameList *p2;
    while (p1 != NULL)
    {
        p2 = p1;
        p1 = p1->next;
        free(p2);
    }

    return 1;
}

/// set flag to best result folder
void SetFlag(char *path, char *pBest)
{
    char SysCmd[MAX_CMD];

    if (NULL == opendir(path))
    {
        printf("路径打开失败\n");
        return;
    }

    memset(SysCmd, 0, sizeof(char) * MAX_CMD);
    sprintf(SysCmd, "mv ./%s ./%sBEST", pBest, pBest);

    system(SysCmd);
}

void main()
{
    char path[] = "/autoMTInv/autoCompute/20180915011312/MTINV";
    char Best[MAX_PATH];

    unsigned int Index;

    Index = BestResult(path, Best);

    if (1 == Index)
    {
        SetFlag(path, Best);
    }
    else
    {
        printf("Didn't find proper folder.");
    }
}
