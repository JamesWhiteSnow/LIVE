#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <fstream>
#include "btree.h"
#include "bnode.h"
#include "bentry.h"
#include "f_def.h"
#include "blk_file.h"

using namespace std;

BTree::BTree(char *_trfname, int _blen, Cache *_cache)
{
	strcpy(tree_name, _trfname);
	char trfname[100];
	strcpy(trfname, _trfname);
	//	strcat(trfname, ".b");
	file = new BlockFile(trfname, _blen);
	cache = _cache;

	// init the first node-----------------------------------------
	root_ptr = new BNode(this);
	root_ptr->level = 0;
	root = root_ptr->block;
	delroot();
	//------------------------------------------------------------
	BNode::instance_cnt = 0;
	num_of_nodes = 1;
	treeloading = false;
}

BTree::BTree(char *_tree_name, Cache *_cache)
{
	strcpy(tree_name, _tree_name);
	char trfname[100];
	strcpy(trfname, _tree_name);
	//	strcat(trfname, ".b");
	file = new BlockFile(trfname, 0);
	cache = _cache;

	root_ptr = NULL;
	char *header = new char[file->get_blocklength()];
	file->read_header(header);
	read_header(header);
	delete[] header;

	BNode::instance_cnt = 0;
	treeloading = false;
}

BTree::BTree(char *_trfname, char *_dsfname, int _blen, Cache *_cache)
{
	strcpy(tree_name, _trfname);
	char trfname[100];
	strcpy(trfname, _trfname);
	//	strcat(trfname, ".b");
	file = new BlockFile(trfname, _blen);
	cache = _cache;

	// init the first node-----------------------------------------
	root_ptr = new BNode(this);
	root_ptr->level = 0;
	root = root_ptr->block;
	delroot();
	//------------------------------------------------------------
	BNode::instance_cnt = 0;
	num_of_nodes = 1;
	treeloading = false;

	// FILE *fp;
	ifstream fp;
	fp.open(_dsfname);
	// if((fp = fopen(_dsfname,"r")) == NULL)
	if (!fp.is_open())
	{
		delete this;
		error(const_cast<char *>("Cannot open text file"), TRUE);
	}
	else
	{
		int cnt = 0;
		// while (!feof(fp))
		while (!fp.eof())
		{
			cnt++;
			if (cnt % 1000 == 0)
			{
				for (int i = 0; i < 30; i++)
					printf("\b");
				printf("Inserted %d", cnt);
				BEntry *e = find_smallest();
				printf("smallest key=%.3f\n", e->key);
				e->son_ptr = NULL;
				delete e;
			}

			// if (cnt==58427)
			//	printf("testing...\n");

			int id;
			double key, tmp;
			// fscanf(fp, "%d %f\n", &id, &key);
			fp >> id >> key >> tmp;
			BEntry *e = new BEntry(this);
			e->son = id;
			e->key = key;
			e->agg = 1;
			insert(e);
			// if (cnt>10724 && !find_key(495.712041))
			// printf("cnt=%d, caught a bug\n", cnt);
		}
	}

	// fclose(fp);
	fp.close();
	printf("\n");
}

//-----------------------------------------------

BTree::~BTree()
{
	if (root_ptr)
		delete root_ptr;
	char *header = new char[file->get_blocklength()];
	write_header(header);
	file->set_header(header);
	delete[] header;

	if (cache)
		cache->flush();

	delete file;
}

void BTree::delroot()
{
	return;
	// delete root_ptr; root_ptr=NULL;  //commented for the main-memory version
}

//-----------------------------------------------

bool BTree::find_key(double _k)
{
	load_root();
	bool ret = root_ptr->find_key(_k);
	// delete root_ptr; root_ptr = NULL;
	return ret;
}

//-----------------------------------------------

void BTree::imt_qry_agg(int _t1, int _t2, int _d1, int _d2, int _daily_tm)
{
	load_root();
	root_ptr->imt_qry_agg(_t1, _t2, _d1, _d2, _daily_tm, B_MAXINT);
	delete root_ptr;
	root_ptr = NULL;
}

void BTree::insert(BEntry *_new_e)
{
	load_root();
	BNode *new_nd = NULL;
	BINSRT c_ret = root_ptr->insert(_new_e, &new_nd);
	if (c_ret == B_OVRFLW)
	{
		BNode *new_root = new BNode(this);
		new_root->level = root_ptr->level + 1;
		new_root->add_new_child(root_ptr);
		new_root->add_new_child(new_nd);
		root = new_root->block;
		root_ptr = new_root;
		num_of_nodes++;
	}
	delroot();
}

//-----------------------------------------------

void BTree::load_root()
{
	if (root_ptr == NULL)
	{
		if (!treeloading)
			error(const_cast<char *>("load_root trying to read from the disk"), true);
		root_ptr = new BNode(this, root);
	}
}

//-----------------------------------------------

void BTree::qry_agg(int _k1, int _k2)
{
	load_root();
	root_ptr->qry_agg(_k1, _k2, B_MAXINT);
	delete root_ptr;
	root_ptr = NULL;
}

//-----------------------------------------------

void BTree::read_header(char *_buf)
{
	int i = 0;
	memcpy(&root, _buf, sizeof(root));
	i += sizeof(root);
	memcpy(&num_of_nodes, &_buf[i], sizeof(num_of_nodes));
	i += sizeof(num_of_nodes);
}

//-----------------------------------------------

void BTree::write_header(char *_buf)
{
	int i = 0;
	memcpy(_buf, &root, sizeof(root));
	i += sizeof(root);
	memcpy(&_buf[i], &num_of_nodes, sizeof(num_of_nodes));
	i += sizeof(num_of_nodes);
}

BEntry *BTree::rank_find(double _rank)
{
	load_root();
	BEntry *cret = root_ptr->rank_find(_rank);
	delroot();
	return cret;
}

bool BTree::delete_entry(BEntry *_be)
{
	bool ret = true;
	load_root();
	BDEL cret = root_ptr->delete_entry(_be);
	if (cret == B_UNDRFLW)
	{
		root = root_ptr->entries[0].son;
		BNode *new_rootptr = root_ptr->entries[0].son_ptr;
		root_ptr->num_entries = 0;
		delete root_ptr;
		root_ptr = new_rootptr;
		num_of_nodes--;
	}
	if (cret == B_NOTFOUND)
		ret = false;
	return ret;
}

void BTree::materialize_to_disk()
{
	traverse_tree();
	printf("materializing the tree to disk\n");
	load_root();
	//	root_ptr->materialize_to_disk();
	delete root_ptr;
	root_ptr = NULL;
}

void BTree::traverse_tree()
{
	load_root();
	double aggrslt = traverse(root_ptr);
	printf("total agg=%.1f\n", aggrslt);
	treeloading = false;
}

BEntry *BTree::find_smallest()
{
	BNode *nd = root_ptr;
	while (nd->level > 0)
	{
		nd = nd->entries[0].son_ptr;
	}
	if (nd->num_entries == 0)
		return NULL;
	BEntry *be = new BEntry(this);
	*be = nd->entries[0];
	return be;
}