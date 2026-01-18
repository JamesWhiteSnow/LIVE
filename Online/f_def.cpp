#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "f_def.h"
#include "bnode.h"
#include "bentry.h"

void error(char *msg, bool _exit)
{
	printf(msg);
	if (_exit)
		exit(1);
}

float traverse(BNode *_bn)
{
	if (_bn->level < 0)
		error(const_cast<char *>("found a neg-level node\n"), true);
	//	printf("block=%d, level=%d\n", _bn->block, _bn->level);

	float ret = 0;
	if (_bn->level == 0)
	{
		for (int i = 0; i < _bn->num_entries; i++)
		{
			ret += _bn->entries[i].agg;
			if (_bn->entries[i].son_ptr)
				error(const_cast<char *>("leaf entry has non-null ptr\n"), true);
		}
		return ret;
	}
	else
	{
		for (int i = 0; i < _bn->num_entries; i++)
		{
			BNode *child = _bn->entries[i].get_son();
			float cret = traverse(child);
			if (fabs(cret - _bn->entries[i].agg) > FLOATZERO)
			{
				printf("%f\t%f\n", cret, _bn->entries[i].agg);
				error(const_cast<char *>("found inconsistency\n"), true);
			}
			ret += cret;
		}
		return ret;
	}
}
