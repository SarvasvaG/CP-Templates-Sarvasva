#include "macros_debug.hpp"

// TC=O(n+m)
ll n;
vector<vector<ll>> adj;
vector<bool> visited;
vector<ll> d, p;
void bfs(ll s)
{
    queue<ll> q;
    q.push(s);
    visited[s] = true;
    p[s] = -1;
    while (!q.empty())
    {
        ll v = q.front();
        q.pop();
        for (auto u : adj[v])
        {
            if (!visited[u])
            {
                visited[u] = true;
                q.push(u);
                d[u] = d[v] + 1;
                p[u] = v;
            }
        }
    }
    // Restoring path from u to source node s using parent vector
}

// TC=O(n+m)
vector<vector<ll>> adj;
vector<bool> visited;
void dfs(ll v)
{
    visited[v] = true;
    for (auto u : adj[v])
    {
        if (!visited[u])
        {
            dfs(u);
        }
    }
}

// Generic DFS implementation
vector<vector<ll>> adj;
vector<ll> color;
vector<ll> time_in, time_out;
ll dfs_timer = 0;
// Color = 0 if unvisited, 1 if visited but not exited, 2 if exited.
void dfs(ll v)
{
    time_in[v] = dfs_timer++;
    color[v] = 1;
    for (auto u : adj[v])
    {
        if (color[u] == 0)
        {
            dfs(u);
        }
    }
    color[v] = 2;
    time_out[v] = dfs_timer++;
}

// Get the number of components of an UNDIRECTED graph
// TC=O(n+m)
vector<ll> comp;
void findComps()
{
    fill(visited.begin(), visited.end(), 0);
    for (ll v = 0; v < n; v++)
    {
        if (!visited[v])
        {
            comp.clear();
            dfs(v);
            cout << "Component:";
            for (auto u : comp)
            {
                cout << ' ' << u;
            }
            cout << endl;
        }
    }
}

// Finding all bridges in an UNDIRECTED graph
// TC=O(n+m)
vector<vector<ll>> adj;
vector<ll> tin;
vector<ll> low;
vector<bool> visited;
ll dfs_timer;
set<pair<ll, ll>> bridges;
void dfsBridge(int v, int par = -1)
{
    visited[v] = 1;
    tin[v] = low[v] = dfs_timer++;
    for (auto &child : adj[v])
    {
        if (child == par)
        {
            continue;
        }
        if (!visited[child])
        {
            dfsBridge(child, v);
            low[v] = min(low[v], low[child]);
            if (tin[v] < low[child])
                bridges.insert({v, child});
        }
        else
        {
            low[v] = min(low[v], tin[child]);
        }
    }
}

void findBridges()
{
    dfs_timer = 0;
    visited.assign(n, false);
    tin.assign(n, -1);
    low.assign(n, -1);
    for (ll i = 0; i < n; ++i)
    {
        if (!visited[i])
        {
            dfsBridge(i);
        }
    }
}

// Using vector may result in repeated vertices.
set<ll> articulation;
void dfsArticulation(int v, int par = -1)
{
    visited[v] = 1;
    tin[v] = low[v] = dfs_timer++;
    int children = 0;
    for (auto &child : adj[v])
    {
        if (child == par)
        {
            continue;
        }
        if (!visited[child])
        {
            dfsArticulation(child, v);
            low[v] = min(low[v], low[child]);
            if (low[child] >= tin[v] && par != -1)
                articulation.insert(v);
            children++;
        }
        else
        {
            low[v] = min(low[v], tin[child]);
        }
    }
    if (par == -1 && children >= 2)
    {
        articulation.insert(v);
    }
}

void findCutpoints()
{
    dfs_timer = 0;
    visited.assign(n, false);
    tin.assign(n, -1);
    low.assign(n, -1);
    for (ll i = 0; i < n; ++i)
    {
        if (!visited[i])
        {
            dfsArticulation(i);
        }
    }
}

/*DIJKSTRA'S ALGORITHM*/
/*Time Complexity: O(mlogn)
Input: Adjacency List, Source Vertex
Output: Distance Vector, Parent Vector to track the shortest path
NOTE: Doesn't work with NEGATIVE-WEIGHT EDGES because we NEVER LOOK BACK*/
pair<vector<ll>, vector<ll>> dijkstra(vector<vector<pair<ll, ll>>> &adj, ll n, ll src)
{
    priority_queue<pair<ll, ll>, vector<pair<ll, ll>>, greater<pair<ll, ll>>> q;
    vector<ll> dist(n, 1e15);
    vector<ll> p(n, -1);

    dist[src] = 0;
    q.push({0, src});

    while (!q.empty())
    {
        ll v = q.top().second;
        ll d = q.top().first;
        q.pop();
        if (dist[v] < d)
        {
            continue;
        }
        for (auto &to : adj[v])
        {
            ll child = to.first;
            ll wt = to.second;
            if (dist[child] > dist[v] + wt)
            {
                dist[child] = dist[v] + wt;
                q.push({dist[child], child});
                p[child] = v;
            }
        }
    }
    return {dist, p};
}

/*BELLMAN-FORD ALGORITHM*/
/*Time Complexity: O(nm)
Input: Edges Vector, Number of vertices, Number of edges, Source Vertex
Output:
vll - dist[i], -INF if i is involved in/ affected by negative cycle
bool - true if negtaive cycle exists
NOTE: For undirected graphs denote (u,v,w) and (v,u,w) in edges vector*/
pair<vector<ll>, bool> bellman_ford(vector<vector<ll>> &edges, ll n, ll m, ll s)
{

    vector<ll> dist(n, INF);
    bool negative_cycle = 0;
    dist[s] = 0;
    for (ll i = 1; i < n; i++)
    {
        for (ll j = 0; j < m; j++)
        {
            ll u = edges[j][0], v = edges[j][1], w = edges[j][2];
            if (dist[u] != INF)
            {
                dist[v] = min(dist[v], dist[u] + w);
            }
        }
    }

    for (ll j = 0; j < m; j++)
    {
        int u = edges[j][0], v = edges[j][1], w = edges[j][2];
        if (dist[u] != INF && dist[u] + w < dist[v])
        {
            negative_cycle = 1;
            dist[v] = -INF;
        }
    }

    return {dist, negative_cycle};
}

// 0-1 BFS
// TC=O(n+m)
vector<vector<pair<ll, ll>>> adj;
vi d;
void bfs01(int s)
{
    deque<ll> q;
    q.push_front(s);

    while (!q.empty())
    {
        ll v = q.front();
        q.pop_front();

        for (auto edge : adj[v])
        {
            ll child = edge.first;
            ll wt = edge.second;

            if (d[v] + wt < d[child])
            {
                d[child] = d[v] + wt;

                if (wt == 1)
                {
                    q.push_back(child);
                }
                else
                {
                    q.push_front(child);
                }
            }
        }
    }
}

/*FLOYD-WARSHALL ALGORITHM*/
/*Time Complexity: O(n^3)
Input: adjacency MATRIX, number of vertices
Output:
vvll - dist[i][j], 1e15 if no path.
bool - true if negative cycle exists */
pair<vector<vector<ll>>, bool> floydWarshall(vector<vector<ll>> &adj, ll n, unordered_map<ll, ll> &mp)
{
    bool negativeCycle = false;
    vector<vector<ll>> dist(n, vector<ll>(n, 0));
    vector<ll> ans;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                dist[i][j] = 0;
            else if (adj[i][j] != -1)
                dist[i][j] = adj[i][j];
            else
                dist[i][j] = INF;
        }
    }
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (dist[i][k] != INF && dist[k][j] != INF)
                {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (dist[i][i] < 0)
        {
            negativeCycle = true;
        }
    }
    return {dist, negativeCycle};
}
/*CYCLE DETECTION IN DIRECTED GRAPH*/
/*Returns True if cycle is found.
NOTE:
-Color=0 for not visited, Color=1 for visited and in current DFS, Color=2 for visited and not in current dfs
-cycle_start and cycle_end are start and end points of any one cycle in the graph.*/
vi color, parent;
vvi adj;
int cycle_start, cycle_end, n, m;
bool dfs_cycle_find(int v)
{
    color[v] = 1;
    bool found = 0;
    for (auto &child : adj[v])
        if (color[child] == 0)
        {
            parent[child] = v;
            if (dfs_cycle_find(child))
                found = 1;
        }
        else if (color[child] == 1)
        {
            cycle_start = child, cycle_end = v;
            found = 1;
        }
    color[v] = 2;
    return found;
}

/*HELPER FUNCTION FOR CYCLE PATH DIRECTED GRAPH*/
/*Returns: {-1} if no cycle is present, [a..b..c..a] if cycle is found*/
vi find_cycle()
{
    parent.assign(n, -1);
    color.assign(n, 0);
    cycle_start = -1;
    fe(i, 0, n, 1) if (color[i] == 0 && dfs_cycle_find(i)) break;
    if (cycle_start == -1)
        return {-1};
    vi path(1, cycle_start);
    while (cycle_end != cycle_start)
    {
        path.pb(cycle_end);
        cycle_end = parent[cycle_end];
    }
    path.pb(cycle_start);
    reverse(all(path));
    return path;
}

// Detect Cycle in UNDIRECTED Graph
// TC=O(n+m)
int n;
vvi adj;
vi color;
vi parent;
int cycle_start, cycle_end;

bool dfs(int v, int par = -1)
{

    visited[v] = 1;
    for (auto child : adj[v])
    {
        if (!visited[child])
        {
            parent[child] = v;
            if (dfs(u, parent[u]))
                return true;
        }
        else if (par != child)
        {
            cycle_end = v;
            cycle_start = child;
            return true;
        }
    }
    return false;
}

void find_cycle()
{
    visited.assign(n, false);
    parent.assign(n, -1);
    cycle_start = -1;

    for (int v = 0; v < n; v++)
    {
        if (!visited[v] && dfs(v, parent[v]))
            break;
    }

    if (cycle_start == -1)
    {
        cout << "Acyclic" << endl;
    }
    else
    {
        vector<int> cycle;
        cycle.push_back(cycle_start);
        for (int v = cycle_end; v != cycle_start; v = parent[v])
            cycle.push_back(v);
        cycle.push_back(cycle_start);

        cout << "Cycle found: ";
        for (int v : cycle)
            cout << v << " ";
        cout << endl;
    }
}

vector<bool> vis;
vi order;
void dfs_topo(int v)
{
    vis[v] = 1;
    for (auto &child : adj[v])
        if (!vis[child])
            dfs_topo(child);
    order.pb(v);
}

/*TOPOLOGICAL-SORT*/
/*Returns: if(DAG) returns Topological Sort, else returns {-1}*/
vector<ll> topologicalSort(ll n)
{
    queue<ll> q;
    vector<ll> indegree(n, 0);
    vector<ll> order(n, 0);

    for (int i = 0; i < n; i++)
    {
        for (auto x : adj[i])
        {
            indegree[x]++;
        }
    }

    for (int i = 0; i < n; i++)
    {
        if (indegree[i] == 0)
        {
            q.push(i);
        }
    }

    while (q.size())
    {
        ll node = q.front();
        q.pop();
        order.push_back(node);
        for (auto child : adj[node])
        {
            if (--indegree[child] == 0)
            {
                q.push(child);
            }
        }
    }
    if (order.size() == n)
        return order;
    return {-1};
}

/*TOPOLOGICAL-SORT*/
/*Returns: Topological sort, if DAG else returns -1*/
vi topological_sort()
{
    color.assign(n, 0);
    parent.assign(n, -1);
    cycle_start = -1;
    fe(i, 0, n, 1) if (color[i] == 0 && dfs_cycle_find(i)) break;
    if (cycle_start != -1)
        return {-1};

    vis.assign(n, 0);
    order.clear();
    fe(i, 0, n, 1) if (!vis[i])
        dfs_topo(i);
    reverse(all(order));
    return order;
}

// Check if the gievn graph is bipartite
// TC=O(m+n)
int n;
vvi adj;
vi side(n, -1);
bool is_bipartite = true;
void bfs(int v)
{

    side[v] = 0;
    queue<int> q;
    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (auto child : adj[u])
        {
            if (side[child] == -1)
            {
                side[child] = side[u] ^ 1;
                q.push(child);
            }

            else
            {
                is_bipartite &= (side[u] != side[child]);
            }
        }
    }
}
void check_bipartite(int n)
{
    fe(i, 0, n, 1)
    {
        if (side[i] == -1)
            bfs(i);
    }
}

/*UNION BY SIZE AND PATH COMPRESSION*/
/*Time Complexity of find: O(alpha(n)) where n = number of vertices and alpha is the inverse-ackermann function*/
class UnionFind
{
    int numberOfNodes, numberOfComponents;
    vector<int> parentOfNode;
    vector<int> sizeOfComponent;

public:
    UnionFind(int numberOfNodes)
        : numberOfNodes(numberOfNodes), numberOfComponents(numberOfNodes)
    {
        parentOfNode.resize(numberOfNodes + 1);
        sizeOfComponent.resize(numberOfNodes + 1);
        for (int node = 0; node < numberOfNodes; node++)
        {
            parentOfNode[node] = node;
            sizeOfComponent[node] = 1;
        }
    }

    void makeSet(int v)
    {
        parentOfNode[v] = v;
        sizeOfComponent[v] = 1;
    }

    int findSet(int v)
    {
        if (v == parentOfNode[v])
            return v;
        return parentOfNode[v] = findSet(parentOfNode[v]);
    }

    /*Returns true, in case union has occured. False, in case they were already in same component */
    bool unionSets(int a, int b)
    {
        a = findSet(a);
        b = findSet(b);
        if (a != b)
        {
            if (sizeOfComponent[a] < sizeOfComponent[b])
                swap(a, b);
            parentOfNode[b] = a;
            sizeOfComponent[a] += sizeOfComponent[b];
            numberOfComponents--;
            return true;
        }
        return false;
    }

    bool isUnited() const { return (numberOfComponents == 1); }
    vector<int> getSizeOfComponent() const { return sizeOfComponent; }
    int componentSize(int node) { return sizeOfComponent[findSet(node)]; }
    int getNumberOfComponents() const { return numberOfComponents; }
};

/*KRUSKAL'S ALGORITHM FOR MST*/
/*Time Complexity: O(mlogm)
Input: Edges Vector - edges[i][0]=u, edges[i][1]=v, edges[i][2]=wt, Number of vertices, Number of edges
Output: MST Edges, Weight of MST
NOTE: If number of MST Edges is NOT n-1 then no MST is possible.*/
pair<vvi, ll> minimum_spanning_tree(vvi edges, int n, int m)
{
    sort(all(edges), [&](vi &a, vi &b)
         { return a[2] < b[2]; });
    UnionFind mst(n);
    ll mst_weight = 0;
    vvi mst_edges;
    fe(i, 0, m, 1)
    {
        int u = edges[i][0], v = edges[i][1], w = edges[i][2];
        if (mst.unionSets(u, v))
        {
            mst_weight += w;
            mst_edges.pb(edges[i]);
        }
    }
    return {mst_edges, mst_weight};
}

int n, l;
vvi adj;
int timer;
vi tin, tout;
vvi up;
void dfs_binary_lifting(int v, int p = -1)
{
    tin[v] = ++timer;
    up[v][0] = p;
    fe(i, 1, l + 1, 1) if (up[v][i - 1] == -1) break;
    else up[v][i] = up[up[v][i - 1]][i - 1];
    for (int child : adj[v])
        if (child != p)
            dfs(child, v);
    tout[v] = ++timer;
}

bool is_ancestor(int u, int v)
{
    if (u == -1 || v == -1)
        return 1;
    return tin[u] <= tin[v] && tout[u] >= tout[v];
}

int lca(int u, int v)
{
    if (is_ancestor(u, v))
        return u;
    if (is_ancestor(v, u))
        return v;
    for (int i = l; i >= 0; --i)
        if (!is_ancestor(up[u][i], v))
            u = up[u][i];
    return up[u][0];
}
void preprocess(int root)
{
    tin.assign(n, 0);
    tout.assign(n, 0);
    timer = 0;
    l = ceil(log2(n));
    up.assign(n, vi(l + 1, -1));
    dfs_binary_lifting(root);
}

/*SUCCESSOR-TABLE FOR FUNCTIONAL GRAPHS*/
/*Time Complexity of Construction: O(n*table_size)
Space Complexity: O(n*table_size)
NOTE: 0-BASED INDEXING is used*/
vvi generate_succesor_table(vi &succ, int n)
{
    int table_size = 30;
    vvi succ_table(n, vi(table_size, 0));
    fe(i, 0, n, 1)
        succ_table[i][0] = succ[i];

    fe(k, 1, table_size, 1)
        fe(x, 0, n, 1)
            succ_table[x][k] = succ_table[succ_table[x][k - 1]][k - 1];
    return succ_table;
}

/*QUERIES FOR FUNCTIONAL GRAPHS*/
/*NOTE: 0-BASED INDEXING is used*/
int query(vvi &succ_table, int node, int jump)
{
    if (jump == 0)
        return node;
    int table_size = sz(succ_table[0]);
    fer(i, table_size - 1, 0, 1) if ((jump >> i) & 1)
        node = succ_table[node][i];
    return node;
}

/*STORNGLY CONNTECTED COMPONENTS UTILITY CLASS*/
/*n,m-> Number of vertices and edges in the ORIGINAL graph.
comp-> Number of strongly connected components in the ORIGINAL graph.
adj-> Adjacency List of the original graph.
adj_rev-> Edges reversed Adjacency List of the original graph.
adj_dag-> Adjacency List for the CONDENSED graph (condensed vertices are labelled starting from 0)
component-> component[i] denotes the component to which vertex i in the original graph belongs ([0,comp-1])
NOTE: TOPOLOGICAL ORDER OF CONDENSED DAG IS 0,1,...,comp-1
*/
class SCC
{
    ll n, m, comp;
    vector<vector<ll>> adj, adj_rev;
    vector<set<ll>> adj_dag;
    vector<bool> vis;
    vector<ll> order, component;
    vector<vector<ll>> component_to_vertex;
    void dfs_topo(ll v)
    {
        vis[v] = 1;
        for (auto &child : adj[v])
            if (!vis[child])
                dfs_topo(child);
        order.push_back(v);
    }
    void dfs_scc(ll v, ll comp)
    {
        vis[v] = 1;
        component[v] = comp;
        for (auto child : adj_rev[v])
            if (!vis[child])
                dfs_scc(child, comp);
    }

public:
    SCC(vector<vector<ll>> &adj, vector<vector<ll>> &adj_rev, ll n, ll m) : adj(adj), adj_rev(adj_rev), n(n), m(m) {};

    /*KOSARAJU'S ALGORITHM*/
    /*Time Complexity: O(n+m)
    Returns: component vector and comp*/
    pair<vector<ll>, ll> strongly_connected_components()
    {
        vis.assign(n, 0);
        order.clear();
        component.assign(n, 0);
        for (int i = 0; i < n; i++)
            if (!vis[i])
                dfs_topo(i);
        reverse(all(order));
        vis.assign(n, 0);
        comp = 0;
        for (auto &node : order)
            if (!vis[node])
                dfs_scc(node, comp++);
        return {component, comp};
    }

    /*CREATES ADJ_DAG LIST*/
    /*Time Complexity: O(m)*/
    vector<set<ll>> condensed_dag()
    {
        adj_dag.assign(comp, set<ll>());
        component_to_vertex.assign(comp, vector<ll>(0));
        for(int node=0;node<n;node++)
        {
            component_to_vertex[component[node]].push_back(node);
            for (auto &child : adj[node])
                if (component[child] != component[node])
                    adj_dag[component[node]].insert(component[child]);
        }
        return adj_dag;
    }
    vector<vector<ll>> get_component_to_vertex()
    {
        return component_to_vertex;
    }
};

int n, m;
vi parent;
vvi adj;
vvll residual_graph, graph;
vector<bool> vis;

/*HELPER FUNCTION FOR AUGMENTING PATH*/
/*Time Complexity: O(n+m)
Input: Source and Destination
Output: Bottleneck in the SHORTEST NUMBER OF EDGES PATH from source to destination. +ve if present, 0 if no more augmenting paths.*/
ll bfs_max_flow(int source, int destination)
{
    parent.assign(n, -1);
    parent[source] = -2;
    queue<pair<int, ll>> q;
    q.push({source, INF});
    while (!q.empty())
    {
        auto [current_vertex, current_flow] = q.front();
        q.pop();
        for (auto &next_vertex : adj[current_vertex])
        {
            if (parent[next_vertex] == -1 && residual_graph[current_vertex][next_vertex])
            {
                parent[next_vertex] = current_vertex;
                ll new_flow = min(current_flow, residual_graph[current_vertex][next_vertex]);
                if (next_vertex == destination)
                    return new_flow;
                q.push({next_vertex, new_flow});
            }
        }
    }
    return 0;
}

/*EDMOND KARP'S MAX-FLOW ALGORITHM*/
/*Time Complexity: O(m^2n)
Input: Source and Destination
Output: Max-Flow in 'graph' from source to destination.
NOTE:
-Upper Bound on the number of augmenting paths = O(nm)*/
ll edmond_max_flow(int source, int destination)
{
    residual_graph = graph;
    ll flow = 0, bottle_neck = 0;
    while (bottle_neck = bfs_max_flow(source, destination))
    {
        flow += bottle_neck;
        int current_vertex = destination;
        while (current_vertex != source)
        {
            int previous_vertex = parent[current_vertex];
            residual_graph[previous_vertex][current_vertex] -= bottle_neck;
            residual_graph[current_vertex][previous_vertex] += bottle_neck;
            current_vertex = previous_vertex;
        }
    }
    return flow;
}

/*HELPER FUNCTION FOR MIN_CUT_EDGES*/
/*Visits all the vertices reachable from the source in the residual graph*/
void dfs_residual(int current_vertex)
{
    vis[current_vertex] = 1;
    for (auto &next_vertex : adj[current_vertex])
        if (!vis[next_vertex] && residual_graph[current_vertex][next_vertex])
            dfs_residual(next_vertex);
}

/*MINIMUM CUT EDGES*/
/*Time Complexity: O(m+n)
Input: Source and Destination
Output: Minimum cut edges in form of edges[i][0]->edges[i][1]
NOTE: Sum of min cut edges' weight = Maximum Flow in the graph*/
vvi min_cut_edges(int source, int destination)
{
    edmond_max_flow(source, destination);
    vis.assign(n, 0);
    vvi edges;
    dfs_residual(source);
    fe(current_vertex, 0, n, 1) for (auto &next_vertex : adj[current_vertex]) if (vis[current_vertex] && !vis[next_vertex] && graph[current_vertex][next_vertex])
        edges.pb({current_vertex, next_vertex});
    return edges;
}

/*HELPER FUNCTION FOR EDGE DISJOINT PATH*/
/*Time Complexity : O(n+m)
Input: source,destination,current_vertex(intially set to source), vis_edges (collection of used up edges), final set of edje disjoint paths (source_to_destination_edp)
Output:
-The function returns edge disjoint path starting from current_vertex and terminating at destination vertex (which is unique for all vertices other than source and destination as edge capacity=1).
-For source, the function combines all the edge disjoint paths starting at source.
-For destination the function simply returns it.*/
vi dfs_edp(int current_vertex, int &source, int &destination, vvi &source_to_destination_edp, set<pii> &vis_edge)
{
    if (current_vertex == destination)
        return {destination};
    for (auto &next_vertex : adj[current_vertex])
    {
        if (graph[current_vertex][next_vertex] && !residual_graph[current_vertex][next_vertex] && vis_edge.find({current_vertex, next_vertex}) == vis_edge.end())
        {
            vis_edge.insert({current_vertex, next_vertex});
            vi next_edp = dfs_edp(next_vertex, source, destination, source_to_destination_edp, vis_edge);
            next_edp.push_back(current_vertex);
            if (current_vertex == source)
            {
                reverse(all(next_edp));
                source_to_destination_edp.push_back(next_edp);
            }
            else
                return next_edp;
        }
    }
    return {};
}

/*EDGE DISJOINT PATHS*/
/*Time Complexity: O(m^2n)+O(m+n)
Input: source and destination
Output: Returns all the edge_disjoint paths*/
vvi edge_disjoint_paths(int source, int destination)
{
    ll max_flow = edmond_max_flow(source, destination);
    set<pii> vis_edge;
    vvi source_to_destination_edp;
    dfs_edp(source, source, destination, source_to_destination_edp, vis_edge);
    return source_to_destination_edp;
}