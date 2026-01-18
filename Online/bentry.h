#ifndef __BENTRY
#define __BENTRY

//-----------------------------------------------

class BTree;
class BNode;

class BEntry
{
public:
	//--==write to disk==--
	double agg;
	double key;
	int son;

	//--==others==--
	BTree *my_tree;
	BNode *son_ptr;

	//--==functions==--
	BEntry();
	BEntry(BTree *_btree);
	~BEntry();
	void del_son();
	static int get_size();
	BNode *get_son();
	void read_from_buffer(char *_buf);
	void set_my_tree(BTree *_btree);
	void write_to_buffer(char *_buf);

	BEntry &operator=(BEntry &_e);
	bool operator==(BEntry &_e);
};

//-----------------------------------------------

#endif