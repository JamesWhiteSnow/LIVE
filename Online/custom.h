#include "graph.h"
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
#include <algorithm>
#include <set>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <chrono>
#include <limits>
#include <omp.h>
#include <unordered_map>
#include <iomanip>
#include <unordered_set>
#include <functional>
#include "graphflow.h"
#include "blk_file.h"
#include "cache.h"
#include "linlist.h"
#include "cdf.h"
#include "gendef.h"
#include "heap.h"
#include "btree.h"
#include "bentry.h"
#include "bnode.h"

#define NANOSECTOSEC(elapsed_time) ((elapsed_time) / (double)1000000000)
#define NANOSECTOMSEC(elapsed_time) ((elapsed_time) / (double)1000000)

using namespace std;
using namespace chrono;

ui spur_dim = 2;
ui span_dim = spur_dim;
ui vde_dim = spur_dim;

ui MAX_LIMIT = 100000;

ui num_hops = 2;

double alpha_value = 1000;
double beta_value = 0.01;

class Vertex
{
public:
	ui vertex_id;
	ui label;
	ui degree;

	vector<ui> hop_neighbor_num;

	vector<double> spur;
	vector<double> span;
	vector<double> spur_norm;
	vector<double> vde;
	vector<double> vde_opt;

	vector<vector<double>> neighbor_spur;
	vector<vector<double>> neighbor_spur_sort;

	vector<vector<double>> hop_range;
	vector<vector<double>> degree_group;

	double key;
};

static inline size_t bytes_of_vector_ui(const std::vector<ui> &v)
{
	return v.size() * sizeof(ui);
}

static inline size_t bytes_of_vector_double(const std::vector<double> &v)
{
	return v.size() * sizeof(double);
}

static inline size_t bytes_of_vector_vector_double(const std::vector<std::vector<double>> &vv)
{
	size_t bytes = vv.size() * sizeof(std::vector<double>);

	for (const auto &inner : vv)
	{
		bytes += inner.size() * sizeof(double);
	}
	return bytes;
}

vector<vector<ui>> BFS_k_hops(const Static_Graph *graph, ui start, int max_hop)
{
	vector<vector<ui>> result(max_hop);

	const ui n = graph->getVerticesCount();
	vector<bool> visited(n, false);

	queue<ui> q;
	queue<int> depth;

	visited[start] = true;
	q.push(start);
	depth.push(0);

	while (!q.empty())
	{
		ui u = q.front();
		q.pop();
		int d = depth.front();
		depth.pop();

		if (d == max_hop)
			continue;

		ui cnt = 0;
		const ui *nbrs = graph->getVertexNeighbors(u, cnt);

		for (ui i = 0; i < cnt; i++)
		{
			ui v = nbrs[i];
			if (!visited[v])
			{
				visited[v] = true;

				if (d < max_hop)
					result[d].push_back(v);

				q.push(v);
				depth.push(d + 1);
			}
		}
	}

	return result;
}

vector<Vertex> gen_vde_data(Static_Graph *data_graph, const vector<vector<double>> &spur_learned, const vector<vector<double>> &spur_learned_norm)
{
	vector<Vertex> data_vertices(data_graph->getVerticesCount());

	for (ui vid = 0; vid < data_graph->getVerticesCount(); vid++)
	{
		data_vertices[vid].vertex_id = vid;
		data_vertices[vid].label = data_graph->getVertexLabel(vid);
		data_vertices[vid].degree = data_graph->getVertexDegree(vid);

		data_vertices[vid].hop_neighbor_num.push_back(0);
		data_vertices[vid].hop_neighbor_num.push_back(data_vertices[vid].degree);

		data_vertices[vid].spur = spur_learned[data_graph->getVertexLabel(vid)];

		data_vertices[vid].spur_norm = spur_learned_norm[data_graph->getVertexLabel(vid)];
	}

	for (ui vid = 0; vid < data_graph->getVerticesCount(); vid++)
	{
		ui nbr_cnt;
		const VertexID *nbrs = data_graph->getVertexNeighbors(vid, nbr_cnt);

		data_vertices[vid].span = vector<double>(span_dim, 0);
		for (ui j = 0; j < nbr_cnt; j++)
		{
			for (ui k = 0; k < span_dim; k++)
			{
				data_vertices[vid].span[k] += data_vertices[nbrs[j]].spur[k];
			}
			data_vertices[vid].neighbor_spur.push_back(data_vertices[nbrs[j]].spur);
		}

		data_vertices[vid].vde = vector<double>();
		for (ui j = 0; j < spur_dim; j++)
		{
			data_vertices[vid].vde.push_back(10 * data_vertices[vid].spur[j]);
		}
		for (ui j = 0; j < span_dim; j++)
		{
			data_vertices[vid].vde[j] += 1 * data_vertices[vid].span[j];
		}
		data_vertices[vid].vde_opt = vector<double>();
		for (ui j = 0; j < spur_dim; j++)
		{
			data_vertices[vid].vde_opt.push_back(alpha_value * data_vertices[vid].spur_norm[j]);
		}
		for (ui j = 0; j < span_dim; j++)
		{
			data_vertices[vid].vde_opt[j] += beta_value * data_vertices[vid].span[j];
		}

		double distance = 0;
		for (ui j = 0; j < span_dim; j++)
		{
			distance += pow(data_vertices[vid].span[j], 2);
		}
		data_vertices[vid].key = alpha_value * data_vertices[vid].label + beta_value * sqrt(distance);

		if (data_vertices[vid].neighbor_spur.size() != 0)
		{
			data_vertices[vid].neighbor_spur_sort.resize(data_vertices[vid].neighbor_spur.size());
			for (ui col = 0; col < spur_dim; col++)
			{
				vector<double> column;
				for (ui row = 0; row < data_vertices[vid].neighbor_spur.size(); row++)
				{
					column.push_back(data_vertices[vid].neighbor_spur[row][col]);
				}
				sort(column.begin(), column.end());
				for (ui row = 0; row < data_vertices[vid].neighbor_spur.size(); row++)
				{
					data_vertices[vid].neighbor_spur_sort[row].push_back(column[row]);
				}
			}
		}

		data_vertices[vid].degree_group.resize(data_vertices[vid].degree + 1);
		for (ui deg = 0; deg <= data_vertices[vid].degree; deg++)
		{
			vector<double> mbr_min;
			vector<double> mbr_max;
			for (ui dim = 0; dim < spur_dim; dim++)
			{
				double min_value = 0;
				double max_value = 0;
				for (ui i = 0; i < deg; i++)
				{
					min_value += data_vertices[vid].neighbor_spur_sort[i][dim];
				}
				for (ui i = data_vertices[vid].neighbor_spur_sort.size() - deg; i < data_vertices[vid].neighbor_spur_sort.size(); i++)
				{
					max_value += data_vertices[vid].neighbor_spur_sort[i][dim];
				}
				mbr_min.push_back(min_value);
				mbr_max.push_back(max_value);
			}
			for (ui i = 0; i < spur_dim; i++)
			{
				data_vertices[vid].degree_group[deg].push_back(mbr_min[i]);
				data_vertices[vid].degree_group[deg].push_back(mbr_max[i]);
			}
		}

		if (data_vertices[vid].neighbor_spur.size() != 0)
		{
			vector<vector<ui>> k_hop_nbrs = BFS_k_hops(data_graph, vid, num_hops);

			data_vertices[vid].hop_range.resize(num_hops + 1);

			for (ui i = 1; i < num_hops; i++)
			{
				if (k_hop_nbrs[i].size() == 0)
				{
					data_vertices[vid].hop_neighbor_num.push_back(0);
					continue;
				}

				data_vertices[vid].hop_neighbor_num.push_back(k_hop_nbrs[i].size());

				for (ui j = 0; j < spur_dim; j++)
				{
					data_vertices[vid].hop_range[i].push_back(data_vertices[k_hop_nbrs[i][0]].spur[j]);
					data_vertices[vid].hop_range[i].push_back(data_vertices[k_hop_nbrs[i][0]].spur[j]);
				}

				for (ui j = 1; j < k_hop_nbrs[i].size(); j++)
				{
					for (ui k = 0; k < spur_dim; k++)
					{
						if (data_vertices[k_hop_nbrs[i][j]].spur[k] < data_vertices[vid].hop_range[i][2 * k])
						{
							data_vertices[vid].hop_range[i][2 * k] = data_vertices[k_hop_nbrs[i][j]].spur[k];
						}
						if (data_vertices[k_hop_nbrs[i][j]].spur[k] > data_vertices[vid].hop_range[i][2 * k + 1])
						{
							data_vertices[vid].hop_range[i][2 * k + 1] = data_vertices[k_hop_nbrs[i][j]].spur[k];
						}
					}
				}
			}
		}
	}

	return data_vertices;
}

vector<Vertex> gen_vde_query(Static_Graph *query_graph, const vector<vector<double>> &spur_learned, const vector<vector<double>> &spur_learned_norm)
{
	vector<Vertex> query_vertices(query_graph->getVerticesCount());

	for (ui vid = 0; vid < query_graph->getVerticesCount(); vid++)
	{
		query_vertices[vid].vertex_id = vid;
		query_vertices[vid].label = query_graph->getVertexLabel(vid);
		query_vertices[vid].degree = query_graph->getVertexDegree(vid);

		query_vertices[vid].hop_neighbor_num.push_back(0);
		query_vertices[vid].hop_neighbor_num.push_back(query_vertices[vid].degree);

		query_vertices[vid].spur = spur_learned[query_graph->getVertexLabel(vid)];

		query_vertices[vid].spur_norm = spur_learned_norm[query_graph->getVertexLabel(vid)];
	}

	for (ui vid = 0; vid < query_graph->getVerticesCount(); vid++)
	{
		ui nbr_cnt;
		const VertexID *nbrs = query_graph->getVertexNeighbors(vid, nbr_cnt);

		query_vertices[vid].span = vector<double>(span_dim, 0);
		for (ui j = 0; j < nbr_cnt; j++)
		{
			for (ui k = 0; k < span_dim; k++)
			{
				query_vertices[vid].span[k] += query_vertices[nbrs[j]].spur[k];
			}
		}
		query_vertices[vid].vde = vector<double>();
		for (ui j = 0; j < spur_dim; j++)
		{
			query_vertices[vid].vde.push_back(10 * query_vertices[vid].spur[j]);
		}
		for (ui j = 0; j < span_dim; j++)
		{
			query_vertices[vid].vde[j] += 1 * query_vertices[vid].span[j];
		}
		query_vertices[vid].vde_opt = vector<double>();
		for (ui j = 0; j < spur_dim; j++)
		{
			query_vertices[vid].vde_opt.push_back(alpha_value * query_vertices[vid].spur_norm[j]);
		}
		for (ui j = 0; j < span_dim; j++)
		{
			query_vertices[vid].vde_opt[j] += beta_value * query_vertices[vid].span[j];
		}

		double distance = 0;
		for (ui j = 0; j < span_dim; j++)
		{
			distance += pow(query_vertices[vid].span[j], 2);
		}
		query_vertices[vid].key = alpha_value * query_vertices[vid].label + beta_value * sqrt(distance);

		if (nbr_cnt != 0)
		{
			vector<vector<ui>> k_hop_nbrs = BFS_k_hops(query_graph, vid, num_hops);

			query_vertices[vid].hop_range.resize(num_hops + 1);

			for (ui i = 1; i < num_hops; i++)
			{
				if (k_hop_nbrs[i].size() == 0)
				{
					query_vertices[vid].hop_neighbor_num.push_back(0);
					continue;
				}

				query_vertices[vid].hop_neighbor_num.push_back(k_hop_nbrs[i].size());

				for (ui j = 0; j < spur_dim; j++)
				{
					query_vertices[vid].hop_range[i].push_back(query_vertices[k_hop_nbrs[i][0]].spur[j]);
					query_vertices[vid].hop_range[i].push_back(query_vertices[k_hop_nbrs[i][0]].spur[j]);
				}

				for (ui j = 1; j < k_hop_nbrs[i].size(); j++)
				{
					for (ui k = 0; k < spur_dim; k++)
					{
						if (query_vertices[k_hop_nbrs[i][j]].spur[k] < query_vertices[vid].hop_range[i][2 * k])
						{
							query_vertices[vid].hop_range[i][2 * k] = query_vertices[k_hop_nbrs[i][j]].spur[k];
						}
						if (query_vertices[k_hop_nbrs[i][j]].spur[k] > query_vertices[vid].hop_range[i][2 * k + 1])
						{
							query_vertices[vid].hop_range[i][2 * k + 1] = query_vertices[k_hop_nbrs[i][j]].spur[k];
						}
					}
				}
			}
		}
	}

	return query_vertices;
}

void btree_range_search(double a, double b,
						BNode *node,
						const Vertex &query_vertex,
						const std::vector<Vertex> &data_vertices,
						std::set<ui> &candidates)
{
	if (node == nullptr)
		return;

	if (node->level == 0)
	{
		for (int i = 0; i < node->num_entries; i++)
		{
			BEntry *entry = &node->entries[i];
			double key = entry->key;

			if (key >= a && key < b)
			{
				ui vid = entry->son;

				bool pruning_success = false;

				if (data_vertices[vid].label != query_vertex.label)
				{
					pruning_success = true;
				}

				if (!pruning_success)
				{
					for (ui hop = num_hops; hop > 1; hop--)
					{
						if (data_vertices[vid].hop_neighbor_num[hop] < query_vertex.hop_neighbor_num[hop])
						{
							pruning_success = true;
							break;
						}
						const auto &dh = data_vertices[vid].hop_range[hop];
						const auto &qh = query_vertex.hop_range[hop];
						if (!dh.empty() && !qh.empty())
						{
							for (ui j = 0; j < vde_dim; j++)
							{
								if (qh[2 * j] > dh[2 * j + 1] || qh[2 * j + 1] < dh[2 * j])
								{
									pruning_success = true;
									break;
								}
							}
							if (pruning_success)
								break;
						}
					}
				}

				if (!pruning_success)
				{
					ui deg = query_vertex.degree;
					const auto &dg = data_vertices[vid].degree_group;
					const auto &qspan = query_vertex.span;
					if (query_vertex.degree > data_vertices[vid].degree)
					{
						pruning_success = true;
					}
					else if (!dg.empty())
					{
						for (ui j = 0; j < vde_dim; j++)
						{
							if (qspan[j] > dg[deg][2 * j + 1] || qspan[j] < dg[deg][2 * j])
							{
								pruning_success = true;
								break;
							}
						}
					}
				}

				if (!pruning_success)
				{
					candidates.insert(vid);
				}
			}
			else if (key >= b)
			{
				break;
			}
		}
		return;
	}

	for (int i = 0; i < node->num_entries; i++)
	{
		double st = node->entries[i].key;
		double ed = (i < node->num_entries - 1)
						? node->entries[i + 1].key
						: std::numeric_limits<double>::infinity();

		if (ed < a || st >= b)
			continue;

		BNode *child = node->entries[i].get_son();
		btree_range_search(a, b, child, query_vertex, data_vertices, candidates);
		node->entries[i].del_son();
	}
}

VertexID selectGQLStartVertex(const Static_Graph *query_graph, ui *candidates_count)
{
	ui start_vertex = 0;

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID cur_vertex = i;

		if (candidates_count[cur_vertex] < candidates_count[start_vertex])
		{
			start_vertex = cur_vertex;
		}
		else if (candidates_count[cur_vertex] == candidates_count[start_vertex] && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex))
		{
			start_vertex = cur_vertex;
		}
	}

	return start_vertex;
}

void updateValidVertices(const Static_Graph *query_graph, VertexID query_vertex, std::vector<bool> &visited,
						 std::vector<bool> &adjacent)
{
	visited[query_vertex] = true;
	ui nbr_cnt;
	const ui *nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i)
	{
		ui nbr = nbrs[i];
		adjacent[nbr] = true;
	}
}

double generateGQLQueryPlan(const Static_Graph *data_graph, const Static_Graph *query_graph, ui *candidates_count, ui *&order, ui *&pivot)
{
	std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
	std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
	order = new ui[query_graph->getVerticesCount()];
	pivot = new ui[query_graph->getVerticesCount()];

	auto start = chrono::high_resolution_clock::now();
	VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
	order[0] = start_vertex;
	updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID next_vertex;
		ui min_value = data_graph->getVerticesCount() + 1;
		for (ui j = 0; j < query_graph->getVerticesCount(); ++j)
		{
			VertexID cur_vertex = j;

			if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex])
			{
				if (candidates_count[cur_vertex] < min_value)
				{
					min_value = candidates_count[cur_vertex];
					next_vertex = cur_vertex;
				}
				else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex))
				{
					next_vertex = cur_vertex;
				}
			}
		}
		updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
		order[i] = next_vertex;
	}

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID u = order[i];
		for (ui j = 0; j < i; ++j)
		{
			VertexID cur_vertex = order[j];
			if (query_graph->checkEdgeExistence(u, cur_vertex))
			{
				pivot[i] = cur_vertex;
				break;
			}
		}
	}
	auto end = chrono::high_resolution_clock::now();
	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

void generateBN(const Static_Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count)
{
	ui query_vertices_num = query_graph->getVerticesCount();
	bn_count = new ui[query_vertices_num];
	std::fill(bn_count, bn_count + query_vertices_num, 0);
	bn = new ui *[query_vertices_num];
	for (ui i = 0; i < query_vertices_num; ++i)
	{
		bn[i] = new ui[query_vertices_num];
	}

	std::vector<bool> visited_vertices(query_vertices_num, false);
	visited_vertices[order[0]] = true;
	for (ui i = 1; i < query_vertices_num; ++i)
	{
		VertexID vertex = order[i];

		ui nbrs_cnt;
		const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
		for (ui j = 0; j < nbrs_cnt; ++j)
		{
			VertexID nbr = nbrs[j];

			if (visited_vertices[nbr] && nbr != pivot[i])
			{
				bn[i][bn_count[i]++] = nbr;
			}
		}

		visited_vertices[vertex] = true;
	}
}

void generateValidCandidates(const Static_Graph *query_graph, const Static_Graph *data_graph, ui depth, vector<ui> &embedding,
							 ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot)
{
	VertexID u = order[depth];
	LabelID u_label = query_graph->getVertexLabel(u);
	ui u_degree = query_graph->getVertexDegree(u);

	idx_count[depth] = 0;

	VertexID p = embedding[pivot[depth]];
	ui nbr_cnt;
	const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i)
	{
		VertexID v = nbrs[i];

		if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
			u_degree <= data_graph->getVertexDegree(v))
		{
			bool valid = true;

			for (ui j = 0; j < bn_cnt[depth]; ++j)
			{
				VertexID u_nbr = bn[depth][j];
				VertexID u_nbr_v = embedding[u_nbr];

				if (!data_graph->checkEdgeExistence(v, u_nbr_v))
				{
					valid = false;
					break;
				}
			}

			if (valid)
			{
				valid_candidate[depth][idx_count[depth]++] = v;
			}
		}
	}
}

double exploreQuickSIStyle(const Static_Graph *data_graph, const Static_Graph *query_graph, ui **candidates,
						   ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &embedding_cnt)
{
	int cur_depth = 0;
	int max_depth = query_graph->getVerticesCount();
	VertexID start_vertex = order[0];

	ui **bn;
	ui *bn_count;

	ui *idx;
	ui *idx_count;
	VertexID **valid_candidate;
	bool *visited_vertices;

	idx = new ui[max_depth];
	idx_count = new ui[max_depth];
	vector<ui> embedding(max_depth);
	visited_vertices = new bool[data_graph->getVerticesCount()];
	std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
	valid_candidate = new ui *[max_depth];

	ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
	for (ui i = 0; i < max_depth; ++i)
	{
		valid_candidate[i] = new VertexID[max_candidate_count];
	}

	idx[cur_depth] = 0;
	idx_count[cur_depth] = candidates_count[start_vertex];
	std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
			  valid_candidate[cur_depth]);

	auto start = chrono::high_resolution_clock::now();

	generateBN(query_graph, order, pivot, bn, bn_count);

	while (true)
	{
		while (idx[cur_depth] < idx_count[cur_depth])
		{
			VertexID u = order[cur_depth];
			VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
			embedding[u] = v;
			visited_vertices[v] = true;
			idx[cur_depth] += 1;

			if (cur_depth == max_depth - 1)
			{
				embedding_cnt += 1;

				visited_vertices[v] = false;
				if (embedding_cnt >= output_limit_num)
				{
					goto EXIT;
				}
			}
			else
			{
				cur_depth += 1;
				idx[cur_depth] = 0;
				generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
										visited_vertices, bn, bn_count, order, pivot);
			}
		}

		cur_depth -= 1;
		if (cur_depth < 0)
			break;
		else
			visited_vertices[embedding[order[cur_depth]]] = false;
	}

EXIT:
	auto end = chrono::high_resolution_clock::now();
	delete[] bn_count;
	delete[] idx;
	delete[] idx_count;
	delete[] visited_vertices;
	for (ui i = 0; i < max_depth; ++i)
	{
		delete[] bn[i];
		delete[] valid_candidate[i];
	}

	delete[] bn;
	delete[] valid_candidate;

	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

double refinement(const Static_Graph *data_graph, const Static_Graph *query_graph, const vector<set<ui>> &candidate_set, ui &answer_num)
{
	size_t call_count = 0;
	size_t output_limit = MAX_LIMIT;
	bool chaoshi = false;

	ui **candidates = new ui *[query_graph->getVerticesCount()];
	ui *candidates_count = new ui[query_graph->getVerticesCount()];
	ui *matching_order = NULL;
	ui *pivots = NULL;

	for (ui i = 0; i < query_graph->getVerticesCount(); i++)
	{
		candidates[i] = new ui[candidate_set[i].size()];
		candidates_count[i] = candidate_set[i].size();
		ui count = 0;
		for (set<ui>::iterator iter = candidate_set[i].begin(); iter != candidate_set[i].end(); iter++)
		{
			candidates[i][count] = *iter;
			count++;
		}
	}

	vector<vector<ui>> answers;

	double query_plan_time = generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);

	size_t embedding_cnt = 0;

	double refinement_time = exploreQuickSIStyle(data_graph, query_graph, candidates, candidates_count, matching_order, pivots, output_limit, embedding_cnt);

	answer_num = embedding_cnt;

	for (ui i = 0; i < query_graph->getVerticesCount(); i++)
	{
		delete[] candidates[i];
	}
	delete[] candidates;
	delete[] candidates_count;
	delete[] matching_order;
	delete[] pivots;

	return query_plan_time + refinement_time;
}
