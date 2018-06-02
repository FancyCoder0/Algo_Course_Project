#include <bits/stdc++.h>
using namespace std;

const int MAX_NODE = 110, MAX_EDGE = 1110, INF = int(1e9);

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
    vector<int> adj[101]; // for each node the vector contains all the edges.
    map<string, int> name_to_id;
    map<string, int> attr_str_to_id;
    void read_from_gxl(const string file) {

        //freopen(file.c_str(), "r", stdin);
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
                    cout << "string=" << attr_str << endl;
                } else 
                {
                    getline(input, tmp);
                    //cout << tmp << endl;
                    int i;
                    for (i = 0; i < tmp.length(); ++i) if (tmp[i] >= '0' && tmp[i] <= '9') break;
                    for (; i < tmp.length() && tmp[i] >= '0' && tmp[i] <= '9'; i++) x = x * 10 + tmp[i] - '0';
                }


                int id = ++n;
                nodes.push_back((node){id, x});
                name_to_id[name] = id;

                //debug
                cout << "node " << id << ' ' << name << endl;
                cout << "attr " << x << endl;
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

                int edge_id = ++m;
                edges.push_back((edge){id1, id2, x});

                adj[id1].push_back(edge_id);
                adj[id2].push_back(edge_id);

                //debug
                cout << "edge " << name1 << ' ' << name2 << endl;
                cout << "attr " << x << endl;
            }

        }

        fclose(stdin);
    }
};

graph origin, target;

namespace KM {
    int get() {
        // todo
        return 0;
    }
};

struct answer {
    int cur_cost, eval_cost;
    vector<int> match; // match: default: 0, empty: -1, other: 1-n
    vector<int> target_map; // for each node in target graph mark the matched node's index in origin graph.

    bool finish() {
        for (int p = 0; p < match.size(); ++p) {
            if (match[p] == 0) {
                return false;
            }
        }
        return true;
    }

    void upd_eval_cost() {
        // todo

        // Build the KM graph
        
        // Calculate the KM result
        eval_cost = KM::get();
    }

    int calc_edit_cost(const int p, const int q) const {
        for (int k = 0; k < origin.adj[p].size(); ++k) {
            auto ed = origin.edges[origin.adj[p][k]];
            // todo
        }
        for (int k = 0; k < target.adj[q].size(); ++k) {
            auto ed = target.edges[target.adj[q][k]];
            // todo
        }
        return 0;
    }

    void print() {
        printf("cur_cost = %d, eval_cost = %d\n", cur_cost, eval_cost);
        printf("Match List:");
        for(int i = 0; i < match.size(); ++i) {
            printf("%d, ", match[i]);
        }
        printf("\n");
        printf("Target Map List:");
        for(int i = 0; i < target_map.size(); ++i) {
            printf("%d, ", target_map[i]);
        }
        printf("\n");
    }

    bool operator<(const answer x) const  {
        return cur_cost + eval_cost > x.cur_cost + x.eval_cost;
    }
};

vector<answer> get_next_list(const answer now) {
    vector<answer> v = {};
    for (int p = 0; p < now.match.size(); ++p) {
        if (now.match[p] == 0) {
            auto ret = now;

            // p -> empty
            ret.match[p] = -1;
            ret.cur_cost += now.calc_edit_cost(p, -1);
            ret.upd_eval_cost();
            v.push_back(ret);

            for (int q = 0; q < now.target_map.size(); ++q) {
                if (now.target_map[q] == 0) {
                    // p -> q
                    ret = now;
                    ret.match[p] = q;
                    ret.target_map[q] = p;
                    ret.cur_cost += now.calc_edit_cost(p, q);
                    ret.upd_eval_cost();
                    v.push_back(ret);
                }
            }
            break;
        }
    }
    return v;
}

int main(int argc, char* argv[]) {
    cost_node_sub = atoi(argv[1]);
    cost_node_di = atoi(argv[2]);
    cost_edge_sub = atoi(argv[3]);
    cost_edge_di = atoi(argv[4]);

    string input1(argv[5]);
    string input2(argv[6]);

    origin.read_from_gxl(input1);
    target.read_from_gxl(input2);

    return 0;

    answer final_ans;
    final_ans.cur_cost = INF;

    priority_queue<answer> que;

    answer empty_answer = (answer) {0, 0, vector<int>((int)origin.nodes.size() + 1), vector<int>{(int)target.nodes.size() + 1}};
    empty_answer.upd_eval_cost();
    que.push(empty_answer);

    while (!que.empty()) {
        auto now = que.top();
        que.pop();

        if (final_ans < now) {
            continue;
        }

        if (now.finish()) {
            final_ans = now;
            continue;
        }

        auto next_list = get_next_list(now);
        for (auto next : next_list) {
            que.push(next);
        }
    }

    // final_ans.print();

    return 0;
}
