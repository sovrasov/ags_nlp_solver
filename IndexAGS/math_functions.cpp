//#include "math_functions.h"
/*
int comp(const void *i, const void *j)
{
    if (((point*)i)->x - ((point*)j)->x<0)
       return -1;
    else
        if (((point*)i)->x - ((point*)j)->x>0)
           return 1;
        else
            return 0;
}

int findmax(double *mem, int size)
{
	int indx=0;
	double max=mem[0];
	for (int i=1; i<size; i++)
	{
		if(mem[i]>max)
		{
			indx=i;
			max=mem[i];
		}
	}
	return indx;
}
int findmin_val(point *mem, int size)
{
	int indx=0;
	double min=mem[0].val;
	for (int i=1; i<size; i++)
		if (mem[i].val<min)
		{
			indx=i;
			min=mem[i].val;
		}
	return indx;
}
void insert(double *mem, int size, int pos, double val)
{
	for (int j = size - 1; j > pos; j--)
		mem[j] = mem[j - 1];
	mem[pos] = val;
}

int sort(point *mem, int size)
{
	point tmp = mem[size - 1];
	int i = 1;
	while (mem[i].x<tmp.x)
		i++;
	for (int j = size - 1; j>i; j--)
		mem[j] = mem[j - 1];
	mem[i] = tmp;
	return i;
}*/
