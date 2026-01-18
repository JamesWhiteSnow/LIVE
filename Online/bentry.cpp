/* this file defines class BEntry */

#include <stdio.h>
#include <memory.h>
#include "bnode.h"
#include "bentry.h"
#include "btree.h"
#include "f_def.h"
//-----------------------------------------------

BEntry::BEntry()
{
	son = key = -1;
	son_ptr = NULL;
}

//-----------------------------------------------

BEntry::BEntry(BTree *_btree)
{
	son = key = -1;
	my_tree = _btree;
	son_ptr = NULL;
}

//-----------------------------------------------

BEntry::~BEntry()
{
	if (son_ptr)
	{
		if (son_ptr->level < 0)
			printf("testing...\n");
		delete son_ptr;
	}
}

//-----------------------------------------------

int BEntry::get_size()
{
	return sizeof(double) + sizeof(double) + sizeof(int);
}

//-----------------------------------------------

BNode *BEntry::get_son()
{
	if (son_ptr == NULL)
	{
		if (!my_tree->treeloading)
			error(const_cast<char *>("get_son trying to read node from disk...\n"), true);
		son_ptr = new BNode(my_tree, son);
	}
	return son_ptr;
}

void BEntry::del_son()
{
	return;
}

//-----------------------------------------------

void BEntry::read_from_buffer(char *_buf)
{
	int i;
	memcpy(&key, _buf, sizeof(key));
	i = sizeof(key);
	memcpy(&son, &_buf[i], sizeof(son));
	i += sizeof(son);
	memcpy(&agg, &_buf[i], sizeof(agg));
	i += sizeof(agg);
}

//-----------------------------------------------

void BEntry::set_my_tree(BTree *_btree)
{
	my_tree = _btree;
}

//-----------------------------------------------

void BEntry::write_to_buffer(char *_buf)
{
	int i;
	memcpy(_buf, &key, sizeof(key));
	i = sizeof(key);
	memcpy(&_buf[i], &son, sizeof(son));
	i += sizeof(son);
	memcpy(&_buf[i], &agg, sizeof(agg));
	i += sizeof(agg);
}

//-----------------------------------------------

BEntry &BEntry::operator=(BEntry &_e)
{
	key = _e.key;
	my_tree = _e.my_tree;
	son = _e.son;
	son_ptr = _e.son_ptr;
	agg = _e.agg;

	return *this;
}

bool BEntry::operator==(BEntry &_e)
{
	if (key != _e.key)
		return false;
	if (son != _e.son)
		return false;
	if (son_ptr != _e.son_ptr)
		return false;
	if (agg != _e.agg)
		return false;

	return true;
}