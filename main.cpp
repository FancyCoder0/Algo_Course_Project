#include <bits/stdc++.h>

using namespace std;

// #define DEBUG
#define HEAP_OPT

#define DELETE -2
#define NOT_MATCH -1
#define PURE_COST 0
#define PREDICT_COST 1

const int MAX_NODE = 110, MAX_EDGE = 1110, INF = int(1e9), BIAS = int(1e5);

int cost_node_sub, cost_node_di, cost_edge_sub, cost_edge_di;

struct node {
    int index, attr;
};

struct edge {
    int x, y, attr;
};

struct graph {
    int n; // #nodes
    int m; // #edges
    int k; // #attr
    vector<node> nodes;
    vector<edge> edges;
    //vector<vector<int> > to_nodes; // for each node the vector contains all the edge's index.
    vector<int> adj[MAX_NODE]; // for each node the vector contains all the edges.
    map<string, int> name_to_id;
    map<string, int> attr_str_to_id;

    int adj_mat[MAX_NODE][MAX_NODE]; // adjacent matrix


    void read_from_gxl(const string file) {

        //freopen(file.c_str(), "r", stdin);
        memset(adj_mat, -1, sizeof(adj_mat));

        ifstream input(file);
        cout << "reading " << file << endl;

        n = m = k = 0;
        string tmp;
        //getline(cin, tmp);
        //cout << tmp << endl;
        while (getline(input, tmp)) {
            //cout << tmp << endl;
            if (tmp.length() < 5) continue;
            if (tmp.substr(1, 4) == "node") {
                // get name
                int name_pos1 = tmp.find('\"');
                int name_pos2 = tmp.find('\"', name_pos1 + 1);
                string name = tmp.substr(name_pos1 + 1, name_pos2 - name_pos1 - 1);


                int x = 0; // attr-id

                // get attr
                if (tmp.find("string") != string::npos)
                {
                    int attr_pos1 = tmp.find("<string>") + 8;
                    int attr_pos2 = tmp.find("</string>");
                    string attr_str = tmp.substr(attr_pos1, attr_pos2 - attr_pos1);
                    if (attr_str_to_id[attr_str] == 0) attr_str_to_id[attr_str] = ++k;
                    x = attr_str_to_id[attr_str];
#ifdef DEBUG
                    cout << "string=" << attr_str << endl;
#endif
                } else
                {
                    getline(input, tmp);
                    //cout << tmp << endl;
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

        //fclose(stdin);
    }
};

graph origin, target;

const int Flow_V = MAX_NODE * 2 + 10, Flow_E = MAX_NODE * MAX_NODE * 2 + 100;
namespace flow {
int edge,S,T,N,fir[Flow_V],e[Flow_E],b[Flow_E],c[Flow_E],dis[Flow_V],de[Flow_V],w[Flow_E];
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
    while(1) {
        for (int i = 0; i < N; ++i) de[i]=INF, v[i]=0, cur[i]=fir[i];
        totflow+=zkw(S,int(1e9));
        int tmp = INF;
        for (int i = 0; i < N; ++i) if(!v[i]) tmp=min(tmp,de[i]);
        if(tmp ==	INF)break;
        for (int i = 0; i < N; ++i) if(!v[i])dis[i]+=tmp;
    }
    return totcost;
}
};

/*
namespace KM {
		int lenx,leny;
		int w[MAX_NODE][MAX_NODE];
		int slack[MAX_NODE],lx[MAX_NODE],ly[MAX_NODE],maty[MAX_NODE];
		bool vx[MAX_NODE],vy[MAX_NODE]; //S集合、Y集合

		bool search(int u) {
		    int i,t;
		    vx[u]=1;
		    for(i=0;i<leny;++i)
		        if(!vy[i]) {
		            t=lx[u]+ly[i]-w[u][i];
		            if (t==0) {
		                vy[i]=1;
		                if(maty[i]==-1||search(maty[i])){
		                    maty[i]=u;
		                    return 1;
		                }
		            }
		            else if(slack[i]>t)
		                slack[i]=t;
		        }
		    return 0;
		}
		int get() {
		    int i,j,ans=0;
		    for(i=0;i<lenx;++i)
		        for(lx[i]=-INF,j=0;j<leny;++j)
		            lx[i]=max(lx[i],w[i][j]);
		    memset(maty,-1,sizeof(maty));
		    memset(ly,0,sizeof(ly));
		    for(i=0;i<lenx;++i) { //找增广路
		        for(j=0;j<leny;++j)
		            slack[j]=INF;
		        while(1) {
		            memset(vx,0,sizeof(vx));
		            memset(vy,0,sizeof(vy));
		            if(search(i))//找到i对应的增广路，不再找
		                break;
		            //没找到增广路，修正
		            int d=INF;
		            for(j=0;j<leny;++j)
		                if(!vy[j]&&d>slack[j])
		                    d=slack[j];
		            for(j=0;j<lenx;++j)
		                if(vx[j])
		                    lx[j]-=d;
		            for(j=0;j<leny;++j)
		                if(vy[j])
		                    ly[j]+=d;
		        }
		    }
		    for(i=0;i<leny;++i)
		        if(maty[i]!=-1)
		            ans+=w[maty[i]][i];
		    return ans;
		}

		void init() {
			lenx = 0;
			leny = 0;
		}
};
*/

struct answer {
    int cur_cost, eval_cost;
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

    int full_match_cost() {
        // a solution's full cost !
        // full match
        //
        int node_sub = 0;
        int node_del = 0;
        int node_ins = 0;
        int edge_sub = 0;
        int edge_del = 0;
        int edge_ins = 0;

        for (int i = 0; i < match.size(); i++)
        {
            if (match[i] == NOT_MATCH) return NOT_MATCH;
            if (match[i] == DELETE)
            {
                node_del++;

                for (int j = 0; j < origin.adj[i].size(); j++)
                {
                    auto ed = origin.edges[origin.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;
                    if (match[k] != DELETE || (match[k] == DELETE && i < k)) // i < k : avoid duplicated calculation
                        edge_del++;
                }
            } else
            {
                node_sub += (origin.nodes[i].attr != target.nodes[match[i]].attr);

                for (int j = 0; j < origin.adj[i].size(); j++)
                {
                    auto ed = origin.edges[origin.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;
                    if (match[k] != DELETE && i < k) // i < k : avoid duplicated calculation
                    {
                        if (target.adj_mat[match[i]][match[k]] != NOT_MATCH)
                        {
                            // edge_sub
                            edge_sub += (ed.attr != target.edges[target.adj_mat[match[i]][match[k]]].attr);
                        } else
                            edge_del++;
                    }
                }
            }
        }

        for (int i = 0; i < target_map.size(); i++)
        {
            if (target_map[i] == NOT_MATCH) // insert
            {
                node_ins++;
                for (int j = 0; j < target.adj[i].size(); j++)
                {
                    auto ed = target.edges[target.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;

                    if (target_map[k] != NOT_MATCH) // match - ins : insert an edge!
                        edge_ins++;
                    if (target_map[k] == NOT_MATCH && i < k) // insert - insert : insert an edge !
                        edge_ins++;
                }
            } else
            {
                //matched node
                for (int j = 0; j < target.adj[i].size(); j++)
                {
                    auto ed = target.edges[target.adj[i][j]];
                    int k = (ed.x == i) ? ed.y : ed.x;

                    if (target_map[k] != NOT_MATCH)
                    {
                        if (origin.adj_mat[target_map[i]][target_map[k]] == NOT_MATCH)  // match-match, not matched edge : edge_ins !
                            edge_ins++;
                    }
                }
            }
        }


        return node_sub * cost_node_sub + (node_ins + node_del) * cost_node_di + edge_sub * cost_edge_sub + (edge_ins + edge_del) * cost_edge_di;
    }


    void upd_eval_cost(answer& final_ans) {
        // Build the flow graph
        vector<int> target_not_match_list;
        vector<int> origin_not_match_list;
        for (int i = 0; i < match.size(); ++i) {
            if (match[i] == NOT_MATCH) {
                origin_not_match_list.push_back(i);
            }
        }
        for (int i = 0; i < target_map.size(); ++i) {
            if (target_map[i] == NOT_MATCH) {
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

        flow::S = lenx + leny + 1;
        flow::T = lenx + leny + 2;
        flow::N = flow::T + 1;
        flow::init();

        for (int i = 0; i < lenx; ++i) {
            flow::add(flow::S, i, (i == x_del) ? INF : 1, (i == x_del) ? 0 : -BIAS);
        }
        for (int i = 0; i < leny; ++i) {
            flow::add(lenx + i, flow::T, (lenx + i == y_del) ? INF: 1, (lenx + i == y_del) ? 0 : -BIAS);
        }

        for (int i = 0; i < leny - 1; ++i) {
            flow::add(x_del, lenx + i, INF, calc_edit_cost(-1, target_not_match_list[i], PREDICT_COST));
        }

        for (int i = 0; i < lenx - 1; ++i) {
            flow::add(i, y_del, INF, calc_edit_cost(origin_not_match_list[i], -1, PREDICT_COST));
        }

        for (int i = 0; i < lenx - 1; ++i) {
            for (int j = 0; j < leny - 1; ++j) {
                flow::add(i, lenx + j, INF, calc_edit_cost(origin_not_match_list[i], target_not_match_list[j], PREDICT_COST));
            }
        }


        eval_cost = flow::solve() + flow::totflow * BIAS;

        assert(flow::totflow == (lenx - 1 + leny - 1));

        answer appro_sol = *this;
        for (int i = 0; i < lenx - 1; ++i) {
            appro_sol.match[origin_not_match_list[i]] = DELETE;
        }
        for (int i = 0; i < lenx - 1; ++i) {
            for (int k = flow::fir[i]; k; k = flow::b[k]) if (flow::c[k^1]) {
                    int j = flow::e[k] - lenx;
                    if (j >= 0 && j < leny) {
                        int x = origin_not_match_list[i], y = target_not_match_list[j];
                        appro_sol.match[x] = y;
                        if (y != DELETE) {
                            appro_sol.target_map[y] = x;
                        }
                    }
                }
        }

        appro_sol.cur_cost = appro_sol.full_match_cost();
        appro_sol.eval_cost = 0;

#ifdef DEBUG
        printf("appro=");
        appro_sol.print();
#endif

        if (appro_sol < final_ans) {
            final_ans = appro_sol;
        }

    }
    /*
    void upd_eval_cost() {
        // Build the KM graph
    		// todo: not consider node -> DELETE
            // todo: consider node DELETE and INSERT -- add more empty node !

        KM::init();
        vector<int> target_not_match_list;
        vector<int> origin_not_match_list;

    			for (int i = 0; i < target_map.size(); ++i) {
    				if (target_map[i] == NOT_MATCH) {
    					KM::leny ++;
    					target_not_match_list.push_back(i);
    				}
    			}
    			for (int i = 0; i < match.size(); ++i) {
    				 if (match[i] == NOT_MATCH) {
    				 		KM::lenx ++;
                            origin_not_match_list.push_back(i);

    				 		for (int j = 0; j < target_not_match_list.size(); ++j) {
                                // min - max
    				 			KM::w[KM::lenx - 1][j] = -calc_edit_cost(i, target_not_match_list[j], PREDICT_COST);


                                #ifdef DEBUG
                                    // if (KM::lenx - 1 == j)
                                    //     cout << "w[" << KM::lenx - 1 << "][" << j << "]=" << KM::w[KM::lenx - 1][j] << endl;
                                #endif
    				 		}
    				 }
    			}
                //origin node -> empty
                if (KM::lenx > KM::leny) {
                    for (int i = 0; i < origin_not_match_list.size(); ++i) {
                        for (int j = KM::leny + 1; j <= KM::lenx; j++)
                            // min - max
                            KM::w[i][j - 1] = -calc_edit_cost(origin_not_match_list[i], -1, PREDICT_COST);
                    }
                    KM::lenx = KM::leny;
                }

                //empty -> target node
                if (KM::lenx < KM::leny) {
                    for (int i = 0; i < target_not_match_list.size(); ++i) {
                        for (int j = KM::lenx + 1; j <= KM::leny; j++)
                            KM::w[j - 1][i] = -calc_edit_cost(-1, target_not_match_list[i], PREDICT_COST);
                    }
                    KM::lenx = KM::leny;
                }

        // Calculate the KM result
        eval_cost = -KM::get();
    }
    */


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

        if (p == -1) {

            pure_cost += cost_node_di; // insert

            for (int k = 0; k < target.adj[q].size(); ++k) {

                auto ed = target.edges[target.adj[q][k]];
                int v = (ed.x == q) ? ed.y : ed.x; // v : q's adjacent node

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
            //
            //
            int deg_p = 0;
            int deg_q = 0;

            map<int, int> p_edge_set;
            map<int, int> q_edge_set;

            for (int k = 0; k < origin.adj[p].size(); ++k) { // scan p's adjacent edge.
                auto ed = origin.edges[origin.adj[p][k]];
                int u = (ed.x == p) ? ed.y : ed.x; // u : p's adjacent node

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
                } else
                    // not matched edge -> delete : a half cost
                    //predict_cost += cost_edge_di / 2;
                {
                    deg_p ++;
                    p_edge_set[ed.attr] += 1;
                }
            }

            for (int k = 0; k < target.adj[q].size(); ++k) {
                auto ed = target.edges[target.adj[q][k]];
                int v = (ed.x == q) ? ed.y : ed.x; // v : p's adjacent node

                if (target_map[v] != NOT_MATCH) {

                    int u = target_map[v];
                    if (origin.adj_mat[u][p] == NOT_MATCH)  // a matched node, but not adjacent to p, full cost
                        pure_cost += cost_edge_di;
                } else
                    //predict_cost += cost_edge_di / 2;
                {
                    deg_q++;
                    q_edge_set[ed.attr] += 1;
                }
            }

            int zero_cost_edge = 0;
            for (auto it = p_edge_set.begin(); it != p_edge_set.end(); ++it)
            {
                int attr = it -> first;
                zero_cost_edge += min(it -> second, q_edge_set[attr]);
            }

            if (deg_p > deg_q)
            {
                predict_cost = (deg_q - zero_cost_edge) * cost_edge_sub + (deg_p - deg_q) * cost_edge_di;
            } else
            {
                predict_cost = (deg_p - zero_cost_edge) * cost_edge_sub + (deg_q - deg_p) * cost_edge_di;
            }

            predict_cost /= 2;
        }

        if (cost_kind == PURE_COST) return pure_cost;
        if (cost_kind == PREDICT_COST) return predict_cost + pure_cost;
    }

    void print() {
        printf("cur_cost = %d, eval_cost = %d, total_cost = %d\n", cur_cost, eval_cost, cur_cost + eval_cost);
        printf("Match List:");
        for(int i = 0; i < match.size(); ++i) {
            printf("%d, ", match[i]);
        }
        printf("\n");
        /*
        printf("Target Map List:");
        for(int i = 0; i < target_map.size(); ++i) {
            printf("%d, ", target_map[i]);
        }
        printf("\n");
        */
    }

    bool operator<(const answer& x) const  {
        return cur_cost + eval_cost < x.cur_cost + x.eval_cost;
    }
    bool operator<=(const answer& x) const {
        return cur_cost + eval_cost <= x.cur_cost + x.eval_cost;
        // return !(x < *this);
    }
};

answer final_ans;

vector<answer> get_next_list(const answer now) {
    vector<answer> v;
    for (int p = 0; p < now.match.size(); ++p) {
        if (now.match[p] == NOT_MATCH) {
            auto ret = now;

            // p -> empty
            ret.match[p] = DELETE;
            ret.cur_cost += now.calc_edit_cost(p, -1, PURE_COST); // delete cost
            ret.upd_eval_cost(final_ans);
            v.push_back(ret);

            for (int q = 0; q < now.target_map.size(); ++q) {
                if (now.target_map[q] == NOT_MATCH) {
                    // p -> q
                    ret = now;
                    ret.match[p] = q;
                    ret.target_map[q] = p;
                    ret.cur_cost += now.calc_edit_cost(p, q, PURE_COST); // subtitute cost
                    ret.upd_eval_cost(final_ans);
                    v.push_back(ret);
                }
            }
            break;
        }
    }
    return v;
}


struct cmp {
    bool operator() (const answer& a, const answer& b) const {
        return b < a; // turn to less heap.
    }
};

int main(int argc, char* argv[]) {
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
        0, 0, vector<int>((int)origin.nodes.size(), NOT_MATCH), vector<int>((int)target.nodes.size(), NOT_MATCH)
    };
    empty_answer.upd_eval_cost(final_ans);
    que.push(empty_answer);

    // empty_answer.print();

    while (!que.empty()) {
        auto now = que.top();
        que.pop();

#ifdef DEBUG
        printf("now:");
        now.print();
        bool ok=1;
        for(int i = 0; i < now.match.size(); ++i) {
            if (now.match[i] != NOT_MATCH)
                if (now.match[i] != i) {
                    ok = 0;
                }
        }
        if (ok) {
            printf("now:");
            now.print();
            printf("full=%d\n",now.full_match_cost());
        }
#endif

#ifdef HEAP_OPT
        if (final_ans.cur_cost <= now.cur_cost) {
            continue;
        }
#endif


        if (now.finish()) {
            // now.cur_cost =  now.full_match_cost();
            // now.eval_cost = 0;
            if (now < final_ans) {
                final_ans = now;
            }
        }


        auto next_list = get_next_list(now);
        for (auto next : next_list) {
            que.push(next);
        }
    }

    final_ans.cur_cost /= 2;
    final_ans.print();

    return 0;
}
