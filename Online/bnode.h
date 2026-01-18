#ifndef _BNODE
#define _BNODE

#include "f_def.h"
#include <vector>
//-----------------------------------------------

class BTree;
class BEntry;

class BNode
{
public:
	//--==disk residential variables==--
	char level;
	int num_entries;
	int sibling;
	BEntry *entries;

	//--==others==--
	bool dirty;
	int block;
	int capacity;
	BTree *my_tree;

	static int instance_cnt;

	//--==functions==--
	BNode(BTree *_btree);
	BNode(BTree *_btree, int _block);
	~BNode();
	void add_new_child(BNode *_cnd);
	bool chk_undrflw();
	bool chk_ovrflw();
	int choose_subtree(BEntry *_new_e);
	int choose_subtree(double _key);
	void enter(BEntry *_new_e);
	bool find_key(double _k);
	int get_header_size();
	BEntry *get_entries(int _cap);
	void imt_qry_agg(int _t1, int _t2, int _d1, int _d2, int _daily_tm, int _r_bnd);
	BINSRT insert(BEntry *_new_e, BNode **_new_nd);
	int max_lesskey_pos(double _key);
	void qry_agg(int _k1, int _k2, int _r_bnd);
	void read_from_buffer(char *_buf);
	void rmv_entry(int _pos);
	void trt_ovrflw(BNode **_new_nd);
	void write_to_buffer(char *_buf);

	BDEL delete_entry(BEntry *_be);
	BEntry *rank_find(double _rank);
	double sumagg();
	void trt_undrflw(int _follow);

	// added for main-memory b-tree
	// void materialize_to_disk();
};

//-----------------------------------------------

#endif