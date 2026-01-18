#include "blk_file.h"
#include "cache.h"
#include "linlist.h"
#include "rand.h"
#include "cdf.h"

#define NOMINMAX
#undef min
#undef max

#include "custom.h"
#include "CLI11.hpp"
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <string>
#include <map>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <filesystem>

using namespace std;
using namespace chrono;

int main(int argc, char *argv[])
{
	string dataset_path = "../Dataset";
	string query_graph_path = "../Dataset/query_graph.graph";

	CLI::App app{"App description"};
	app.add_option("-d,--data", dataset_path, "dataset path")->required(false);
	app.add_option("-q,--query", query_graph_path, "query graph path")->required(false);
	app.add_option("-a,--alpha", alpha_value, "alpha value")->required(false);
	app.add_option("-b,--beta", beta_value, "beta value")->required(false);
	app.add_option("-e,--dimension", spur_dim, "embedding dimension")->required(false);
	app.add_option("-k,--khop", num_hops, "k hop neighbor")->required(false);

	CLI11_PARSE(app, argc, argv);

	span_dim = spur_dim;
	vde_dim = spur_dim;

	string data_graph_path = dataset_path + "/data_graph.graph";

	Static_Graph *data_graph_s = new Static_Graph(true);
	data_graph_s->loadGraphFromFile(data_graph_path);
	data_graph_s->printGraphMetaData();

	vector<vector<double>> learned_spurs(data_graph_s->getLabelsCount(), vector<double>(spur_dim));
	ifstream fin(dataset_path + "/spur_emb.txt");
	for (ui i = 0; i < data_graph_s->getLabelsCount(); i++)
	{
		for (ui j = 0; j < spur_dim; j++)
		{
			fin >> learned_spurs[i][j];
		}
	}
	fin.close();
	vector<vector<double>> learned_spurs_norm(data_graph_s->getLabelsCount(), vector<double>(spur_dim));
	fin.open(dataset_path + "/spur_emb_norm.txt");
	for (ui i = 0; i < data_graph_s->getLabelsCount(); i++)
	{
		for (ui j = 0; j < spur_dim; j++)
		{
			fin >> learned_spurs_norm[i][j];
		}
	}
	fin.close();

	auto start_time = chrono::high_resolution_clock::now();
	vector<Vertex> data_vertices = gen_vde_data(data_graph_s, learned_spurs, learned_spurs_norm);
	auto end_time = chrono::high_resolution_clock::now();
	double data_embedding_time = NANOSECTOMSEC(double(chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count()));

	int block_len = 4096;
	Cache *cache = NULL;
	string btree_path = dataset_path + "/btree";
	std::filesystem::remove(btree_path);
	BTree *btree = new BTree(const_cast<char *>(btree_path.c_str()), block_len, cache);

	start_time = chrono::high_resolution_clock::now();
	for (ui i = 0; i < data_graph_s->getVerticesCount(); i++)
	{
		BEntry *e = new BEntry(btree);
		e->key = data_vertices[i].key;
		e->son = data_vertices[i].vertex_id;
		e->agg = 1;
		btree->insert(e);
	}
	end_time = chrono::high_resolution_clock::now();
	double index_build_time = NANOSECTOMSEC(double(chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count()));

	string query_path = query_graph_path;

	Static_Graph *query_graph_s = new Static_Graph(true);
	query_graph_s->loadGraphFromFile(query_path);

	start_time = chrono::high_resolution_clock::now();
	vector<Vertex> query_vertices = gen_vde_query(query_graph_s, learned_spurs, learned_spurs_norm);
	end_time = chrono::high_resolution_clock::now();
	double query_embedding_time = NANOSECTOMSEC(double(chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count()));

	vector<set<ui>> candidate_sets(query_graph_s->getVerticesCount());

	vector<double> vertex_search_times;

#pragma omp parallel for
	for (ui qvid = 0; qvid < query_graph_s->getVerticesCount(); qvid++)
	{
		btree->load_root();
		BNode *root_node = btree->root_ptr;

		start_time = chrono::high_resolution_clock::now();
		btree_range_search(query_vertices[qvid].key, (query_vertices[qvid].label + 1) * alpha_value, root_node, query_vertices[qvid], data_vertices, candidate_sets[qvid]);
		end_time = chrono::high_resolution_clock::now();
		vertex_search_times.push_back(NANOSECTOMSEC(double(chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count())));

		btree->delroot();
	}
	double candidate_search_time = NANOSECTOMSEC(*max_element(vertex_search_times.begin(), vertex_search_times.end()));

	ui answer_num = 0;
	double refine_time = NANOSECTOMSEC(refinement(data_graph_s, query_graph_s, candidate_sets, answer_num));

	double query_time = query_embedding_time + candidate_search_time + refine_time;

	delete query_graph_s;

	cout << "index build time: " << index_build_time + data_embedding_time << " ms" << endl;
	cout << "query time: " << query_time << " ms" << endl;

	return 0;
}