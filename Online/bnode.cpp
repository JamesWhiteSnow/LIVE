#include <memory.h>
#include "bnode.h"
#include "bentry.h"
#include "btree.h"
#include "f_def.h"
#include "blk_file.h"
#include "cache.h"
#include <algorithm>
//-----------------------------------------------

BNode::BNode(BTree *_btree)
{
	my_tree = _btree;
	int b_length = my_tree->file->get_blocklength();
	int entry_size = BEntry::get_size();
	capacity = (b_length - get_header_size()) / entry_size;

	level = -1;
	sibling = -1;
	num_entries = 0;
	entries = get_entries(capacity);

	char *blk = new char[b_length];
	block = my_tree->file->append_block(blk);
	delete[] blk;

	dirty = true;
	BNode::instance_cnt++;
}

//-----------------------------------------------

BNode::BNode(BTree *_btree, int _block)
{
	my_tree = _btree;
	num_entries = 0;
	dirty = false;

	int header_size = get_header_size();
	int b_len = my_tree->file->get_blocklength();
	capacity = (b_len - header_size) / BEntry::get_size();

	entries = get_entries(capacity);

	block = _block;
	char *blk = new char[b_len];
	if (my_tree->cache == NULL) // no cache
		my_tree->file->read_block(blk, block);
	else
		my_tree->cache->read_block(blk, block, my_tree);

	read_from_buffer(blk);
	delete[] blk;

	BNode::instance_cnt++;
}

//-----------------------------------------------

BNode::~BNode()
{
	//	printf("destronying node %d at level %d\n", block, level);
	// if (level<0)
	//	printf("testing...\n");

	char *buf;
	if (dirty)
	{
		buf = new char[my_tree->file->get_blocklength()];
		write_to_buffer(buf);

		if (my_tree->cache == NULL) // no cache
			my_tree->file->write_block(buf, block);
		else
			my_tree->cache->write_block(buf, block, my_tree);

		delete[] buf;
	}
	for (int i = num_entries; i < capacity; i++)
		entries[i].son_ptr = NULL;
	delete[] entries;

	BNode::instance_cnt--;
}

//-----------------------------------------------
void BNode::add_new_child(BNode *_cnd)
{
	BEntry *e = new BEntry(my_tree);
	e->son = _cnd->block;
	e->key = _cnd->entries[0].key;
	e->agg = _cnd->sumagg();
	e->son_ptr = _cnd;
	enter(e);
	e->son_ptr = NULL; // necessary... deleting e without this line would also free the memory for cnd
	delete e;
	// delete _cnd; //commented for the main-memory version
}

//-----------------------------------------------

bool BNode::chk_ovrflw()
{
	if (num_entries == capacity - 1)
		return true;
	else
		return false;
}

//-----------------------------------------------

bool BNode::chk_undrflw()
{
	if (my_tree->root_ptr->level == level)
		if (num_entries == 1 && level > 0)
			return true;
		else
			return false;
	if (num_entries <= (capacity - 2) / 2)
		return true;
	else
		return false;
}

//-----------------------------------------------

int BNode::choose_subtree(BEntry *_new_e)
{
	double key = _new_e->key;
	int follow = max_lesskey_pos(key);
	if (follow == -1)
	{
		entries[0].key = key;
		dirty = true;
		follow = 0;
	}
	return follow;
}

//-----------------------------------------------

int BNode::choose_subtree(double _key)
{
	double key = _key;
	int follow = max_lesskey_pos(key);
	if (follow == -1)
	{
		entries[0].key = key;
		dirty = true;
		follow = 0;
	}
	return follow;
}

//-----------------------------------------------

void BNode::enter(BEntry *_new_e)
{
	double key = _new_e->key;
	int pos = max_lesskey_pos(key);
	pos++; // this is the position of the new key
	for (int i = num_entries; i > pos; i--)
		entries[i] = entries[i - 1];

	entries[pos] = *_new_e;
	num_entries++;
	dirty = true;
}

//-----------------------------------------------

bool BNode::find_key(double _k)
{
	bool ret = false;

	if (level == 0)
	{
		for (int i = 0; i < num_entries; i++)
			if (entries[i].key == _k)
			{
				// if (entries[i].son==10724)
				ret = true;
				i = num_entries;
			}
		return ret;
	}

	int follow = max_lesskey_pos(_k);
	if (follow != -1)
	{
		BNode *succ = entries[follow].get_son();
		ret = succ->find_key(_k);
		entries[follow].del_son();
	}

	return ret;
}

//-----------------------------------------------

BEntry *BNode::get_entries(int _cap)
{
	BEntry *en = new BEntry[_cap];
	for (int i = 0; i < _cap; i++)
	{
		en[i].set_my_tree(my_tree);
		en[i].son_ptr = NULL;
	}
	return en;
}

//-----------------------------------------------

int BNode::get_header_size()
{
	return sizeof(char) + sizeof(int) + sizeof(int);
}

//-----------------------------------------------

void BNode::imt_qry_agg(int _t1, int _t2, int _d1, int _d2,
						int _daily_tm, int _r_bnd)
{
	if (level == 0)
		return;
	for (int i = 0; i < num_entries; i++)
	{
		bool visit = false;
		int st = entries[i].key, ed;
		if (i < num_entries - 1)
			ed = entries[i + 1].key;
		else
			ed = _r_bnd;
		for (int j = _d1; j <= _d2; j++)
		{
			int k1 = j * _daily_tm + _t1,
				k2 = j * _daily_tm + _t2;
			if (!(k1 <= st && ed <= k2) &&
				(std::max(k1, st) <= std::min(k2, ed)))
			{
				visit = true;
				j = _d2;
			}
		}
		if (visit)
		{
			BNode *succ = entries[i].get_son();
			succ->imt_qry_agg(_t1, _t2, _d1, _d2, _daily_tm, ed);
			entries[i].del_son();
		}
	}
}

//-----------------------------------------------

BINSRT BNode::insert(BEntry *_new_e, BNode **_new_nd)
{
	BINSRT ret = B_NRML;
	if (level == 0)
	{
		enter(_new_e);
		delete _new_e;
		if (chk_ovrflw())
		{
			trt_ovrflw(_new_nd);
			ret = B_OVRFLW;
		}
		return ret;
	}

	int follow = choose_subtree(_new_e);
	BNode *succ = entries[follow].get_son();
	BNode *new_nd = NULL;
	BINSRT c_ret = succ->insert(_new_e, &new_nd);
	entries[follow].agg = succ->sumagg();
	entries[follow].del_son();
	dirty = true;
	if (c_ret == B_OVRFLW)
		add_new_child(new_nd);
	if (chk_ovrflw())
	{
		trt_ovrflw(_new_nd);
		ret = B_OVRFLW;
	}
	return ret;
}

//-----------------------------------------------

int BNode::max_lesskey_pos(double _key)
{
	int pos = -1;
	for (int i = num_entries - 1; i >= 0; i--)
		if (entries[i].key <= _key)
		{
			pos = i;
			i = -1;
		}
	return pos;
}

//-----------------------------------------------

void BNode::read_from_buffer(char *_buf)
{
	int i;
	memcpy(&level, _buf, sizeof(level));
	i = sizeof(level);
	memcpy(&num_entries, &_buf[i], sizeof(num_entries));
	i += sizeof(num_entries);
	memcpy(&sibling, &_buf[i], sizeof(sibling));
	i += sizeof(sibling);

	for (int j = 0; j < num_entries; j++)
	{
		entries[j].read_from_buffer(&_buf[i]);
		i += BEntry::get_size();
	}
}

//-----------------------------------------------

void BNode::qry_agg(int _k1, int _k2, int _r_bnd)
{
	if (level == 0)
		return;
	for (int i = 0; i < num_entries; i++)
	{
		int st = entries[i].key, ed;
		if (i < num_entries - 1)
			ed = entries[i + 1].key;
		else
			ed = _r_bnd;
		if (!(_k1 <= st && ed <= _k2) &&
			(std::max(_k1, st) <= std::min(_k2, ed)))
		{
			BNode *succ = entries[i].get_son();
			succ->qry_agg(_k1, _k2, ed);
			entries[i].del_son();
		}
	}
}

//-----------------------------------------------

void BNode::rmv_entry(int _pos)
{
	for (int i = _pos; i < num_entries - 1; i++)
		entries[i] = entries[i + 1];
	num_entries--;
	dirty = true;
}

//-----------------------------------------------

void BNode::trt_ovrflw(BNode **_new_nd)
{
	*_new_nd = new BNode(my_tree);
	(*_new_nd)->level = level;
	(*_new_nd)->sibling = sibling;
	sibling = (*_new_nd)->block;
	my_tree->num_of_nodes++;

	int i = (capacity - 1) / 2;

	//	while (i < num_entries)
	//	{
	//		(*_new_nd) -> enter(&(entries[i]));
	//		rmv_entry(i);
	//	}

	//	while (entries[i].key==entries[i+1].key && i<num_entries-1)
	//		i++;
	//	if (i==num_entries-1)
	//		error("too many duplicates in the dataset\n", true);

	for (int j = i + 1; j < num_entries; j++)
		(*_new_nd)->enter(&(entries[j]));
	num_entries = i + 1;
	for (int j = num_entries; j < capacity; j++)
		entries[j].son_ptr = NULL;
}

//-----------------------------------------------

void BNode::write_to_buffer(char *_buf)
{
	int i;
	memcpy(_buf, &level, sizeof(level));
	i = sizeof(level);
	memcpy(&_buf[i], &num_entries, sizeof(num_entries));
	i += sizeof(num_entries);
	memcpy(&_buf[i], &sibling, sizeof(sibling));
	i += sizeof(sibling);

	for (int j = 0; j < num_entries; j++)
	{
		entries[j].write_to_buffer(&_buf[i]);
		i += BEntry::get_size();
	}
}

//-----------------------------------------------

double BNode::sumagg()
{
	double ret = 0;
	for (int i = 0; i < num_entries; i++)
		ret += entries[i].agg;
	return ret;
}

BEntry *BNode::rank_find(double _rank)
{
	double acc_rank = 0;
	if (level == 0)
	{
		for (int i = 0; i < num_entries; i++)
		{
			acc_rank += entries[i].agg;
			if (acc_rank >= _rank)
			{
				BEntry *e = new BEntry();
				*e = entries[i];
				e->agg = acc_rank;
				return e;
			}
		}
		error(const_cast<char *>("the entry not found\n"), true);
		return NULL;
	}
	for (int i = 0; i < num_entries; i++)
	{
		acc_rank += entries[i].agg;
		if (acc_rank >= _rank)
		{
			acc_rank -= entries[i].agg;
			BNode *succ = entries[i].get_son();
			BEntry *cret = succ->rank_find(_rank - acc_rank);
			entries[i].del_son();
			if (cret)
			{
				cret->agg += acc_rank;
				return cret;
			}
		}
	}
	return NULL;
}

BDEL BNode::delete_entry(BEntry *_be)
{
	BDEL ret = B_NOTFOUND;
	if (level == 0)
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == *_be)
			{
				delete _be;
				rmv_entry(i);
				if (chk_undrflw())
					ret = B_UNDRFLW;
				else
					ret = B_NONE;
				break; // break for
			}
		}
		return ret;
	}

	int follow = choose_subtree(_be->key);
	BDEL cret = B_NOTFOUND;
	bool again = true, smallonce = false;
	while (cret == B_NOTFOUND && again)
	{
		BNode *succ = entries[follow].get_son();
		cret = succ->delete_entry(_be);
		entries[follow].key = succ->entries[0].key;
		entries[follow].agg = succ->sumagg();
		entries[follow].del_son();
		dirty = true;

		if (follow > 0 && entries[follow - 1].key == entries[follow].key)
			follow--;
		else if (follow > 0 && !smallonce)
		{
			follow--;
			smallonce = true;
		}
		else
			again = false;
	}

	if (cret == B_NONE)
		ret = B_NONE;
	if (cret == B_UNDRFLW)
	{
		trt_undrflw(follow);
		if (chk_undrflw())
			ret = B_UNDRFLW;
		else
			ret = B_NONE;
	}
	return ret;
}

void BNode::trt_undrflw(int _follow)
{
	int mergesub = _follow + 1; // the subscript of the non-leaf entry to merge with
	if (_follow == num_entries - 1)
		mergesub = _follow - 1;

	BNode *succ1 = entries[mergesub].get_son();
	BNode *succ2 = entries[_follow].get_son();

	int totalnum = succ1->num_entries + succ2->num_entries;
	if (totalnum >= capacity - 1) // the resulting entries do not fit in one node
	{
		int n = succ1->num_entries;
		int avgnum = totalnum / 2;

		if (mergesub > _follow) // moving from the node with big keys to the one with small keys
		{
			for (int i = 0; i < totalnum / 2 - succ2->num_entries; i++)
			// totalnum/2-succ2->num_entries gives how many entries to be moved from succ1 to succ2
			{
				succ2->enter(&(succ1->entries[0]));
				succ1->rmv_entry(0);
			}
			//			while (succ1->entries[0].key==succ2->entries[succ2->num_entries-1].key)
			//			{
			//				succ2->enter(&(succ1->entries[0]));
			//				succ1->rmv_entry(0);
			//			}
		}
		else // moving from the node with small keys to the one with big keys
		{
			for (int i = totalnum / 2; i < n; i++)
			{
				succ2->enter(&(succ1->entries[totalnum / 2]));
				succ1->rmv_entry(totalnum / 2);
			}
			//			while (succ1->entries[succ1->num_entries-1].key==succ2->entries[0].key)
			//			{
			//				succ2->enter(&(succ1->entries[succ1->num_entries-1]));
			//				succ1->rmv_entry(succ1->num_entries-1);
			//			}
		}
		entries[mergesub].key = succ1->entries[0].key;
		entries[_follow].key = succ2->entries[0].key;
		entries[mergesub].agg = succ1->sumagg();
		entries[_follow].agg = succ2->sumagg();
		entries[mergesub].del_son();
		entries[_follow].del_son();
	}
	else
	{
		// copy all the entries from succ2 to succ1
		for (int i = 0; i < succ2->num_entries; i++)
		{
			succ1->enter(&(succ2->entries[i]));
		}
		entries[mergesub].key = succ1->entries[0].key;
		entries[mergesub].agg = succ1->sumagg();
		entries[mergesub].del_son();
		entries[_follow].del_son();

		// free the memory for the child node of follow
		for (int i = 0; i < entries[_follow].son_ptr->num_entries; i++)
			entries[_follow].son_ptr->entries[i].son_ptr = NULL;
		entries[_follow].son_ptr->dirty = false; // no need to write disk
		delete entries[_follow].son_ptr;
		my_tree->num_of_nodes--;

		rmv_entry(_follow);
	}
}

/*
void BNode::materialize_to_disk()
{
//	printf("level=%d\n", level);
//	if (level==1)
//		 printf("testing...\n");
	for (int i=0; i<num_entries; i++)
	{
		if (entries[i].son_ptr)
		{
			if (level>1)
			{
				BNode *cnd=entries[i].get_son();
				cnd->materialize_to_disk();
			}
			delete entries[i].son_ptr;
			entries[i].son_ptr=NULL;
		}
	}
}
*/

int BNode::instance_cnt;