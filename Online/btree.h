/* this file defines class BTree */

#ifndef __BTREE
#define __BTREE

#include "f_def.h"
#include "cache.h"
//-----------------------------------------------

class BlockFile;
class BEntry;
class BNode;
class Cache;

class BTree : public Cacheable
{
public:
	//--==write to disk==--
	int root;

	//--==others==--
	bool treeloading;
	char tree_name[100];
	int num_of_nodes;
	BNode *root_ptr;

	//--==functions==--
	BTree(char *_tree_name, Cache *_cache);
	BTree(char *_trfname, int _blen, Cache *_cache);
	BTree(char *_trfname, char *_dsfname, int _blen, Cache *_cache);
	~BTree();
	void delroot();
	bool find_key(double _k);
	void imt_qry_agg(int _t1, int _t2, int _d1, int _d2, int _daily_tm);
	void insert(BEntry *_new_e);
	void load_root();
	void qry_agg(int _k1, int _k2);
	void read_header(char *_buf);
	void write_header(char *_buf);

	BEntry *rank_find(double _rank);
	bool delete_entry(BEntry *_be);

	//for memory version -------------------
	void materialize_to_disk();
	void traverse_tree();
	BEntry *find_smallest();
};

//-----------------------------------------------

#endif