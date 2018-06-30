#include <bits/stdc++.h>

#include <pthread.h>

using namespace std;

//#define DEBUG
#define debug
#define BETTER
//#define PRINT_UPD
//#define PRINT_THREAD

#define optH
#define TIMELIMIT 30
#define EARLY_TERM 0.1
#define PARALLEL
#define NUM_THREADS 4
#define PARALLEL_TASK_LIMIT 100

#define DELETE -2
#define NOT_MATCH -1
#define PURE_COST 0
#define PREDICT_COST 1

#define PREDICT_OPT

#define deln(x) cerr << #x << " = " << x << endl

const int MAX_NODE = 80, MAX_EDGE = 1110, INF = int(1e9), BIAS = int(1e5);
int sb = 0;
int min_cost = 10000;

pthread_mutex_t ans_mutex = PTHREAD_MUTEX_INITIALIZER;

int cost_node_sub, cost_node_di, cost_edge_sub, cost_edge_di;

int random(int x) {
    return rand() % x;
}

struct node {
    int index, attr;
};

struct edge {
    int x, y, attr;
};

map<string, int> attr_str_to_id;
map<int, string> id_to_attr_str;
int attr_num = 0;
vector<int> sorted_list;

struct graph {
    int n; // #nodes
    int m; // #edges
    int k; // #attr
    vector<node> nodes;
    vector<edge> edges;
    //vector<vector<int> > to_nodes; // for each node the vector contains all the edge's index.
    vector<int> adj[MAX_NODE]; // for each node the vector contains all the edges.
    map<string, int> name_to_id;
    int adj_mat[MAX_NODE][MAX_NODE]; // adjacent matrix
    vector<int> adj_H[MAX_NODE];
    bool is_H[MAX_NODE];

    void read_from_gxl(const string file) {

        //freopen(file.c_str(), "r", stdin);
        memset(adj_mat, -1, sizeof(adj_mat));

        ifstream input(file);

        n = m = 0; 
        string tmp;

#ifdef DEBUG
        cout << "reading " << file << endl;
        getline(cin, tmp);
        cout << tmp << endl;
#endif
        while (getline(input, tmp)) {
#ifdef DEBUG
            cout << tmp << endl;
#endif
            if (tmp.length() < 5) continue;
            if (tmp.substr(1, 4) == "node") {
                // get name
                int name_pos1 = tmp.find('\"');
                int name_pos2 = tmp.find('\"', name_pos1 + 1);
                string name = tmp.substr(name_pos1 + 1, name_pos2 - name_pos1 - 1);


                int x = 0; // attr-id

                // get attr
                if (tmp.find("string") != string::npos) {
                    int attr_pos1 = tmp.find("<string>") + 8;
                    int attr_pos2 = tmp.find("</string>");
                    string attr_str = tmp.substr(attr_pos1, attr_pos2 - attr_pos1);
                    if (attr_str_to_id[attr_str] == 0) attr_str_to_id[attr_str] = ++attr_num, id_to_attr_str[attr_num] = attr_str;
                    x = attr_str_to_id[attr_str];
#ifdef DEBUG
                    cout << "string=" << attr_str << " id = " << x << endl;
#endif
                } else {
                    getline(input, tmp);
                    int i;
                    for (i = 0; i < tmp.length(); ++i) if (tmp[i] >= '0' && tmp[i] <= '9') break;
                    for (; i < tmp.length() && tmp[i] >= '0' && tmp[i] <= '9'; i++) x = x * 10 + tmp[i] - '0';
                }


                int id = n++;
                nodes.push_back((node) {
                    id, x
                });
                name_to_id[name] = id;

#ifdef DEBUG
                cout << "node " << id << ' ' << name << endl;
                cout << "attr " << x << endl;
#endif
            }

            if (tmp.substr(1, 4) == "edge") {
                // get name1
                int sb = tmp.find("from");
                int name_pos1 = tmp.find('\"', sb);
                int name_pos2 = tmp.find('\"', name_pos1 + 1);
                string name1 = tmp.substr(name_pos1 + 1, name_pos2 - name_pos1 - 1);

                // get name2
                name_pos1 = tmp.find('\"', name_pos2 + 1);
                name_pos2 = tmp.find('\"', name_pos1 + 1);
                string name2 = tmp.substr(name_pos1 + 1, name_pos2 - name_pos1 - 1);

                int id1 = name_to_id[name1], id2 = name_to_id[name2];

                // get attr
                if (tmp.find("attr") == string::npos) getline(input, tmp);
                int x = 0;
                int i = tmp.find("attr");
                for (; i < tmp.length(); ++i) if (tmp[i] >= '0' && tmp[i] <= '9') break;
                for (; i < tmp.length() && tmp[i] >= '0' && tmp[i] <= '9'; i++) x = x * 10 + tmp[i] - '0';

                int edge_id = m++;
                edges.push_back((edge) {
                    id1, id2, x
                });

                adj[id1].push_back(edge_id);
                adj[id2].push_back(edge_id);

                adj_mat[id1][id2] = adj_mat[id2][id1] = edge_id; // update adjacent matrix

#ifdef DEBUG
                cout << "edge " << name1 << ' ' << name2 << endl;
                cout << "attr " << x << endl;
#endif
            }

        }

        // find H !
        for (int i = 0; i < nodes.size(); i++)
            if (id_to_attr_str[nodes[i].attr] == "H") {
                is_H[i] = true;
                auto ed = edges[adj[i][0]];
                int j = (ed.x == i) ? ed.y : ed.x;
                adj_H[j].push_back(i);
            } else 
                is_H[i] = false;

    }
};

graph origin, target;

const int Flow_V = MAX_NODE * 2 + 10, Flow_E = MAX_NODE * MAX_NODE * 2 + 100;
struct flow {
    int edge,S,T,N,fir[Flow_V],e[Flow_E],b[Flow_E],c[Flow_E],w[Flow_E],dis[Flow_V],de[Flow_V];
    int totflow,totcost;
    int cur[Flow_V], q[Flow_E*10];
    bool v[Flow_V],o[Flow_V];

    void init() {
        edge = 1;
        for (int i = 0; i < N; ++i) fir[i] = 0;
    }
    void add2(int x,int y,int z,int q) {
        e[++edge] = y;
        c[edge] = z;
        w[edge] = q;
        b[edge] = fir[x];
        fir[x] = edge;
    }
    void add(int x,int y,int z,int q) {
        add2(x,y,z,q);
        add2(y,x,0,-q);
    }
    void spfa() {
        int i,j,k,u;
        for(int i = 0; i < N; ++i) dis[i]=INF, v[i]=0;
        q[1]=S;
        dis[S]=0;
        v[S]=1;
        for(i=j=1; u=q[i],i<=j; v[u]=0,i++)
            for(k=fir[u]; k; k=b[k])if(c[k])
                    if(dis[u]+w[k]<dis[e[k]]) {
                        dis[e[k]]=dis[u]+w[k];
                        if(!v[e[k]]) {
                            v[e[k]]=1;
                            q[++j]=e[k];
                        }
                    }
    }
    int zkw(int i,int flow) {
        int d, r=flow, l;
        if(i==T) {
            totcost+=dis[i]*flow;
            return flow;
        }
        v[i]=o[i]=1;
        for(int&k=cur[i]; k; k=b[k])
            if(c[k]) {
                l=dis[i]+w[k]-dis[e[k]];
                de[e[k]]=min(de[e[k]],l);
                if(l==0&&!o[e[k]]) {
                    d=zkw(e[k],min(c[k],r));
                    c[k]-=d;
                    c[k^1]+=d;
                    r-=d;
                    if(r==0)break;
                }
            }
        o[i]=0;
        return flow-r;
    }
    int solve() {
        spfa();
        for (int i = 0; i < N; ++i) v[i] = 0, o[i] = 0;
        totcost=totflow=0;
        while(dis[T] < 0) {
            for (int i = 0; i < N; ++i) de[i]=INF, v[i]=0, cur[i]=fir[i];
            totflow+=zkw(S, INF);
            int tmp = INF;
            for (int i = 0; i < N; ++i) if(!v[i]) tmp=min(tmp,de[i]);
            if(tmp == INF)break;
            for (int i = 0; i < N; ++i) if(!v[i]) dis[i]+=tmp;
        }
        return totcost;
    }
} flow_g[NUM_THREADS];

struct answer {
    int cur_cost, eval_cost, appro_final_cost, last_match;
    vector<int> match; // match: NOT_MATCH: -1, DELETE: -2, other: 0 ~ n - 1
    vector<int> target_map; // for each node in target graph mark the matched node's index in origin graph.

    bool finish() {
        for (int p = 0; p < match.size(); ++p) {
            if (match[p] == NOT_MATCH) {
                return false;
            }
        }
        return true;
    }

    int full_match_cost(int is_final = 0) {
        // a solution's full cost !
        // full match
        int node_sub = 0;
        int node_del = 0;
        int node_ins = 0;
        int node_match = 0;
        int edge_sub = 0;
        int edge_del = 0;
        int edge_ins = 0;
        int edge_match = 0;

        map<string, int> match_kind;
        map<pair<int, int>, int> edge_match_list;

        for (int i = 0; i < match.size(); i++) {

            if (match[i] == NOT_MATCH) return NOT_MATCH;

            if (match[i] == DELETE) {
                node_del++;

                for (int j = 0; j < origin.adj[i].size(); j++) {
                    auto ed = origin.edges[origin.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;
                    if (match[k] != DELETE || (match[k] == DELETE && i < k)) // i < k : avoid duplicated calculation
                        edge_del++;
                }
            } else {
#ifdef debug
                if (origin.nodes[i].attr == target.nodes[match[i]].attr)
                    match_kind[id_to_attr_str[origin.nodes[i].attr]]++;
#endif

                node_sub += (origin.nodes[i].attr != target.nodes[match[i]].attr);
                node_match += (origin.nodes[i].attr == target.nodes[match[i]].attr);

                for (int j = 0; j < origin.adj[i].size(); j++)
                {
                    auto ed = origin.edges[origin.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;
                    if (match[k] != DELETE && i < k) // i < k : avoid duplicated calculation
                    {
                        if (target.adj_mat[match[i]][match[k]] != NOT_MATCH)
                        {
                            if (ed.attr == target.edges[target.adj_mat[match[i]][match[k]]].attr)
                                edge_match_list[make_pair(origin.nodes[i].attr, origin.nodes[k].attr)]++;
                            // edge_sub
                            edge_sub += (ed.attr != target.edges[target.adj_mat[match[i]][match[k]]].attr);
                            edge_match += (ed.attr == target.edges[target.adj_mat[match[i]][match[k]]].attr);
                        } else
                            edge_del++;
                    }
                }
            }
        }

        for (int i = 0; i < target_map.size(); i++) {
            if (target_map[i] == NOT_MATCH) {

                node_ins++;

                for (int j = 0; j < target.adj[i].size(); j++) {
                    auto ed = target.edges[target.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;

                    if (target_map[k] != NOT_MATCH) // match - ins : insert an edge!
                        edge_ins++;
                    if (target_map[k] == NOT_MATCH && i < k) // insert - insert : insert an edge !
                        edge_ins++;
                }
            } else {
                //matched node
                for (int j = 0; j < target.adj[i].size(); j++) {
                    auto ed = target.edges[target.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;

                    if (target_map[k] != NOT_MATCH) {
                        if (origin.adj_mat[target_map[i]][target_map[k]] == NOT_MATCH)  // match-match, not matched edge : edge_ins !
                            edge_ins++;
                    }
                }
            }
        }
#ifdef DEBUG

        if (is_final == 1) {
            cout << "node_sub = " << node_sub << endl;
            cout << "node_del = " << node_del << endl;
            cout << "node_ins = " << node_ins << endl;
            cout << "node_match = " << node_match << endl;
            cout << "edge_sub = " << edge_sub << endl;
            cout << "edge_del = " << edge_del << endl;
            cout << "edge_ins = " << edge_ins << endl;
            cout << "edge_match = " << edge_match << endl;

            for (auto it = match_kind.begin(); it != match_kind.end(); ++it) {
                cout << it -> first << "_match = " << it -> second << endl;
            }
            for (auto it = edge_match_list.begin(); it != edge_match_list.end(); ++it) {
                cout << id_to_attr_str[it -> first.first] << '-' << id_to_attr_str[it -> first.second]<< " match = " << it -> second << endl;
            }
        }
#endif
        int all = node_sub * cost_node_sub + (node_ins + node_del) * cost_node_di + edge_sub * cost_edge_sub + (edge_ins + edge_del) * cost_edge_di;
#ifdef DEBUG
        if (all < min_cost) {
            min_cost = all;
            int mc = min_cost / 2; 
            deln(mc);
        }
#endif
        return all;
    }

    void upd_eval_cost(const int thread_id, answer& final_ans) {
        if (finish()) {
            cur_cost = full_match_cost();
            appro_final_cost = cur_cost;
            eval_cost = 0;
            if (*this < final_ans) {
                pthread_mutex_lock(&ans_mutex);
                if (*this < final_ans) final_ans = *this;
                #ifdef PRINT_UPD
                    printf("upd cost = %d\n", final_ans.cur_cost);
                #endif
                pthread_mutex_unlock(&ans_mutex);
            }
            return;
        }

        // Build the flow graph
        vector<int> target_not_match_list;
        vector<int> origin_not_match_list;
        for (int i = 0; i < match.size(); ++i) {
            if (match[i] == NOT_MATCH) {
#ifdef optH
                if (origin.is_H[i]) continue; // 
#endif 
                
                origin_not_match_list.push_back(i);
            }
        }
        for (int i = 0; i < target_map.size(); ++i) {
            if (target_map[i] == NOT_MATCH) {
#ifdef optH
                if (target.is_H[i]) continue; // 
#endif 

                target_not_match_list.push_back(i);
            }
        }

        int lenx, leny;

        origin_not_match_list.push_back(DELETE);
        lenx = origin_not_match_list.size();

        target_not_match_list.push_back(DELETE);
        leny = target_not_match_list.size();

        int x_del = lenx - 1;
        int y_del = lenx + leny - 1;

        flow* g = &flow_g[thread_id];

        g->S = lenx + leny + 1;
        g->T = lenx + leny + 2;
        g->N = g->T + 1;
        g->init();

        for (int i = 0; i < lenx; ++i) {
            g->add(g->S, i, (i == x_del) ? INF : 1, (i == x_del) ? 0 : -BIAS);
        }
        for (int i = 0; i < leny; ++i) {
            g->add(lenx + i, g->T, (lenx + i == y_del) ? INF: 1, (lenx + i == y_del) ? 0 : -BIAS);
        }

        for (int i = 0; i < leny - 1; ++i) {
            g->add(x_del, lenx + i, INF, calc_edit_cost(-1, target_not_match_list[i], PREDICT_COST));
        }

        for (int i = 0; i < lenx - 1; ++i) {
            g->add(i, y_del, INF, calc_edit_cost(origin_not_match_list[i], -1, PREDICT_COST));
        }

        for (int i = 0; i < lenx - 1; ++i) {
            for (int j = 0; j < leny - 1; ++j) {
                g->add(i, lenx + j, INF, calc_edit_cost(origin_not_match_list[i], target_not_match_list[j], PREDICT_COST));
            }
        }

        eval_cost = g->solve() + (lenx - 1 + leny - 1) * BIAS;

        answer appro_sol = *this;
        for (int i = 0; i < lenx - 1; ++i) {
            appro_sol.match[origin_not_match_list[i]] = DELETE;
        }
        for (int i = 0; i < lenx - 1; ++i) {
            for (int k = g->fir[i]; k; k = g->b[k]) if (g->c[k^1]) {
                    int j = g->e[k] - lenx;
                    if (j >= 0 && j < leny) {
                        int x = origin_not_match_list[i], y = target_not_match_list[j];
                        appro_sol.match[x] = y;
                        if (y != DELETE) {
                            appro_sol.target_map[y] = x;
                        }
                    }
                }
        }


#ifdef optH

        int cost_of_H = 0;

        vector<int> not_match_H;
        bool label[MAX_NODE];
        memset(label, 0, sizeof(label));
        //add H !
        for (int i = 0; i < origin.nodes.size(); i++)
            if (!origin.is_H[i]) {
                int v = appro_sol.match[i];
                for (int j = 0; j < origin.adj_H[i].size(); j++) {
                    if (v != DELETE && j < target.adj_H[v].size()) {
                        appro_sol.match[origin.adj_H[i][j]] = target.adj_H[v][j];
                        appro_sol.target_map[target.adj_H[v][j]] = origin.adj_H[i][j];
                    }
                    else 
                        not_match_H.push_back(origin.adj_H[i][j]);
                }
            }

        int j = 0;
        for (int i = 0; i < not_match_H.size(); i++) {
            for (; j < target.nodes.size(); j++)
                if (target.is_H[j] && target_map[j] == NOT_MATCH)
                    break;

            if (j < target.nodes.size()) {
                appro_sol.match[not_match_H[i]] = j, 
                target_map[j] = not_match_H[i];
                cost_of_H += 2 * cost_edge_di;
            } else 
            {
                appro_sol.match[not_match_H[i]] = DELETE;
                cost_of_H += cost_edge_di + cost_node_di;
            }
        } 

        for (int i = 0; i < target.nodes.size(); i++) {
            if (target.is_H[i] && appro_sol.target_map[i] == NOT_MATCH)
                cost_of_H += cost_edge_di + cost_node_di;
        }

        eval_cost += cost_of_H;
#endif 

        appro_sol.cur_cost = appro_sol.full_match_cost();
        appro_sol.appro_final_cost = appro_sol.cur_cost;
        appro_sol.eval_cost = 0;

        if (match.size() >= 40) {
        	eval_cost = max((appro_sol.cur_cost - cur_cost) * 5 / 10, eval_cost);
        }

        appro_final_cost = appro_sol.cur_cost;
        

#ifdef DEBUG
        printf("appro=");
        appro_sol.print();
#endif

        if (appro_sol < final_ans) {
            pthread_mutex_lock(&ans_mutex);
            if (appro_sol < final_ans) final_ans = appro_sol;
            #ifdef PRINT_UPD
                printf("upd cost = %d\n", final_ans.cur_cost);
            #endif
            pthread_mutex_unlock(&ans_mutex);

        }

        
      	//remove a some matching
        answer better = appro_sol;
        //better.fix(thread_id, final_ans);
        bool changed = false;
        for (int i = 0; i < match.size(); i++) {
        	if (match[i] == NOT_MATCH) {
                if (random(3) < 2) { // remove node : 2/3 probability
                    changed = true;
                    if (better.match[i] == DELETE)
                        better.match[i] = NOT_MATCH;
                    else if (better.match[i] >= 0) {
                        better.target_map[better.match[i]] = NOT_MATCH;
                        better.match[i] = NOT_MATCH;
                    }
                }
            }
        }
        //upd cost
        if (changed) {
        	better.upd_eval_cost(thread_id, final_ans);
        }

    }

    int calc_edit_cost(const int p, const int q, const int cost_kind) const {
        // p = -1 (insert q)
        // q = -1 (delete p)
        // p != -1 and q != -1 (p -> q mapping)
        // cost_kind : a) PURE_COST b) PREDICT_COST
        // PURE_COST : only count the known edge : two nodes are both known
        // PREDICT_COST : consider potential cost
        //
        // VERSION 2 : predict_cost lower bound

        int pure_cost = 0;
        int predict_cost = 0;

#ifdef DEBUG
		cerr << "calc_edit_cost (" << p << ", " << q << ")" << endl;
#endif

        if (p == -1) {


            pure_cost += cost_node_di; // insert

            for (int k = 0; k < target.adj[q].size(); ++k) {

                auto ed = target.edges[target.adj[q][k]];
                int v = (ed.x == q) ? ed.y : ed.x; // v : q's adjacent node

#ifdef optH
                if (target.is_H[v]) continue;
#endif

                if (target_map[v] == NOT_MATCH) {
                    // if adjacent node has not been matched, assign a half 'edge_insert' cost to current node
                    predict_cost += cost_edge_di / 2;
                } else {
                    // if adjacent node has been matched or delete, assign all 'edge_insert' cost to current node
                    pure_cost += cost_edge_di;
                }
            }
        } else if (q == -1) {

            pure_cost += cost_node_di; // delete

            for (int k = 0; k < origin.adj[p].size(); ++k) {

                auto ed = origin.edges[origin.adj[p][k]];
                int v = (ed.x == p) ? ed.y : ed.x; // v : p's adjacent node

#ifdef optH
                if (origin.is_H[v]) continue;
#endif

                if (match[v] == NOT_MATCH) {
                    // if adjacent node has not been matched, assign a half 'edge_delete' cost to current node
                    predict_cost += cost_edge_di / 2;
                } else {
                    // if adjacent node has been matched, assign all 'edge_delete' cost to current node
                    pure_cost += cost_edge_di;
                }
            }
        } else {

            pure_cost += (origin.nodes[p].attr == target.nodes[q].attr) ? 0 : cost_node_sub;

            // count matched edge!
            //int matched_edge = 0;
            //int matched_edge_cost = 0;
            int deg_p = 0;
            int deg_q = 0;

            map<int, int> p_edge_set;
            map<int, int> q_edge_set;

            for (int k = 0; k < origin.adj[p].size(); ++k) { // scan p's adjacent edge.
                auto ed = origin.edges[origin.adj[p][k]];
                int u = (ed.x == p) ? ed.y : ed.x; // u : p's adjacent node
#ifdef optH
                if (origin.is_H[u]) continue; // 
#endif

                if (match[u] != NOT_MATCH) {

                    int v = match[u]; // v : u's match node

                    if (v != DELETE && target.adj_mat[v][q] != NOT_MATCH) { // v is adjacent to q : find a matched edge!

                        //matched_edge++;

                        if (target.edges[target.adj_mat[v][q]].attr == ed.attr)
                            pure_cost += 0;
                        else
                            pure_cost += min(cost_edge_sub, 2 * cost_edge_di); // a matched edge , but attr is not the same : sub or del

                    } else {
                        // adjacent to a matched node, but have to delete edge : full cost
                        pure_cost += cost_edge_di;
                    }
                } else { 
                    // not matched edge -> delete : a half cost
                    //predict_cost += cost_edge_di / 2;
                    deg_p ++;
                    p_edge_set[ed.attr] += 1;
                }
            }

            for (int k = 0; k < target.adj[q].size(); ++k) {
                auto ed = target.edges[target.adj[q][k]];
                int v = (ed.x == q) ? ed.y : ed.x; // v : p's adjacent node
#ifdef optH
                if (target.is_H[v]) continue;
#endif

                if (target_map[v] != NOT_MATCH) {

                    int u = target_map[v];
                    if (origin.adj_mat[u][p] == NOT_MATCH)  // a matched node, but not adjacent to p, full cost
                        pure_cost += cost_edge_di;
                } else {
                    //predict_cost += cost_edge_di / 2;
                    deg_q++;
                    q_edge_set[ed.attr] += 1;
                }
            }

            int zero_cost_edge = 0;
            for (auto it = p_edge_set.begin(); it != p_edge_set.end(); ++it) { 
                int attr = it -> first;
                zero_cost_edge += min(it -> second, q_edge_set[attr]);
            }

            if (deg_p > deg_q) { 
                predict_cost = (deg_q - zero_cost_edge) * cost_edge_sub + (deg_p - deg_q) * cost_edge_di;
            } else { 
                predict_cost = (deg_p - zero_cost_edge) * cost_edge_sub + (deg_q - deg_p) * cost_edge_di;
            }

            predict_cost /= 2;
        }

        if (cost_kind == PURE_COST) return pure_cost;
        if (cost_kind == PREDICT_COST) return predict_cost + pure_cost;
    }

    void print() {
        printf("cost = %d\n", cur_cost);
        printf("Match List:");
        for(int i = 0; i < match.size(); ++i) {
            printf("%d, ", match[i] == DELETE ? -1 : match[i]);
        }
        printf("\n"); 
    }

    bool operator<(const answer& x) const  {
        return cur_cost + eval_cost < x.cur_cost + x.eval_cost;
    }
    bool operator<=(const answer& x) const {
        return cur_cost + eval_cost <= x.cur_cost + x.eval_cost;
    }
};

answer final_ans;

vector<answer> get_next_list(const int thread_id, const answer now) {
    vector<answer> v;
    for (int i = now.last_match + 1; i < sorted_list.size(); ++i) { 

        int p = sorted_list[i];

        if (origin.is_H[p] == 1) continue;

        if (now.match[p] == NOT_MATCH) {
            auto ret = now;

            // p -> empty
            ret.match[p] = DELETE;
            ret.cur_cost += now.calc_edit_cost(p, -1, PURE_COST); // delete cost
            ret.last_match = i;
            ret.upd_eval_cost(thread_id, final_ans);
            v.push_back(ret);

            for (int q = 0; q < now.target_map.size(); ++q) {
                if (now.target_map[q] == NOT_MATCH) {
#ifdef optH
                    if (target.is_H[q]) continue;
#endif                        
                    // p -> q
                    ret = now;
                    ret.match[p] = q;
                    ret.target_map[q] = p;
                    ret.cur_cost += now.calc_edit_cost(p, q, PURE_COST); // subtitute cost
                    ret.last_match = i;
                    ret.upd_eval_cost(thread_id, final_ans);
                    v.push_back(ret);
                }
            }
            break;
        }
    }
    return v;
}

vector<int> generate_sorted_list(answer &now) {

    vector<pair<int, int> > list;
    for (int i = 0; i < now.match.size(); i++)  {
        if (now.match[i] == DELETE)
            list.push_back(make_pair(now.calc_edit_cost(i, -1, PURE_COST), i));
        else 
            list.push_back(make_pair(now.calc_edit_cost(i, now.match[i], PURE_COST), i));
    }

    sort(list.begin(), list.end());

    vector<int> ret;

    for (int i = 0; i < list.size(); i++) {
        ret.push_back(list[i].second);
    }
    return ret;
}

struct cmp {
    bool operator() (const answer& a, const answer& b) const {
        return b < a; // turn to less heap.
    }
};

struct task_args {
    int thread_id;
    vector<answer> tasks;
};

void* run(void* args) {
    task_args* st = (task_args *)args;
    priority_queue<answer, vector<answer>, cmp> que;
    for (auto& x : (*st).tasks) {
        que.push(x);
    }

    int thread_id = (*st).thread_id;

#ifdef PRINT_THREAD
    printf("thread %d start!\n", thread_id);
#endif

    while (!que.empty()) {
        auto now = que.top();
        que.pop();

        if (final_ans <= now) {
            continue;
        }

        for (auto& next : get_next_list(thread_id, now)) {
            if (final_ans <= next) {
                continue;
            }
            que.push(next);
        }
    }

#ifdef PRINT_THREAD
    // cerr << "early stop!" << endl;
    printf("thread %d is finished!\n", thread_id);
#endif
}


int main(int argc, char* argv[]) {
    auto start_point = chrono::system_clock::now();

    //srand(12345);
    srand(time(0));

    cost_node_sub = atoi(argv[1]) * 2;  // convenient for divide 2
    cost_node_di = atoi(argv[2]) * 2;
    cost_edge_sub = atoi(argv[3]) * 2;
    cost_edge_di = atoi(argv[4]) * 2;

    cost_edge_sub = min(cost_edge_sub, 2 * cost_edge_di); // del-ins < sub

    string input1(argv[5]);
    string input2(argv[6]);

    origin.read_from_gxl(input1);
    target.read_from_gxl(input2);

    final_ans.cur_cost = INF;

    priority_queue<answer, vector<answer>, cmp> que;

    answer empty_answer = (answer) {
        0, 0, 0, -1, vector<int>((int)origin.nodes.size(), NOT_MATCH), vector<int>((int)target.nodes.size(), NOT_MATCH)
    };
    empty_answer.upd_eval_cost(0, final_ans);

    // generate sorted list
    sorted_list = generate_sorted_list(final_ans);

    que.push(empty_answer);

    int main_iter_times = 0;
    
    while (!que.empty()) {
#ifdef PARALLEL
        if (que.size() >= PARALLEL_TASK_LIMIT) {
            break;
        }
#else
        if ((++main_iter_times) % 100 == 0) {
            double spend_time = (chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start_point)).count() / 1e3;
            if (abs(TIMELIMIT - spend_time) < EARLY_TERM) {
                break;
            }
            printf("%d : %d\n", main_iter_times, final_ans.cur_cost);
        }
#endif

        auto now = que.top();
        que.pop();

        if (final_ans <= now) {
            continue;
        }

        for (auto& next : get_next_list(0, now)) {
            if (final_ans <= next) {
                continue;
            }
            que.push(next);
        }
    }

#ifdef PARALLEL
    double time_before_parallel = (chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start_point)).count() / 1e3;
    // printf("time_before_parallel time = %.3lfs\n",time_before_parallel);
    // printf("start parallel search....\n");

    pthread_t tids[NUM_THREADS];
    task_args task[NUM_THREADS];

    int p = 0;
    while(!que.empty()) {
        task[p].tasks.push_back(que.top());
        que.pop();
        p = (p + 1) % NUM_THREADS;
    }

    for(int i = 0 ; i < NUM_THREADS; ++i) {
        task[i].thread_id = i;
        int ret = pthread_create(&tids[i], NULL, run, &task[i]);
        if (ret != 0) {
            printf("pthread_create error: error_code = %d", ret);
        }
    }

    this_thread::sleep_for(std::chrono::milliseconds(int((TIMELIMIT - time_before_parallel - EARLY_TERM) * 1e3)));

#endif

    pthread_mutex_lock(&ans_mutex);
    final_ans.cur_cost /= 2;
    final_ans.print();
    final_ans.full_match_cost(1);
    pthread_mutex_unlock(&ans_mutex);

    cout << "run time = " << (chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start_point)).count() / 1e3 << "s" << endl;

    return 0;
}
