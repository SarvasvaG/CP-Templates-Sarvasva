#include <bits/stdc++.h>
#include "macros_debug.hpp"

using namespace std;

struct Node
{
    ll m;
    Node() : m(INF) {};
    Node(ll m) : m(m) {};

    Node operator+(const Node &other)
    {
        Node res;
        res.m = min(m, other.m);
        return res;
    }
};

template <typename T>
class SegmentTree
{
public:
    vector<T> segment;
    vector<ll> array;
    ll arraySize;

    SegmentTree(vector<ll> array, ll arraySize) : array(array), arraySize(arraySize)
    {
        segment.assign(4 * arraySize + 1, T());
        build(0, 0, arraySize - 1);
    }

    void build(int index, int low, int high)
    {
        if (low == high)
        {
            segment[index] = T(array[low]);
            return;
        }
        ll mid = low + (high - low) / 2;
        build(2 * index + 1, low, mid);
        build(2 * index + 2, mid + 1, high);
        segment[index] = segment[2 * index + 1] + segment[2 * index + 2];
    }

    T query(int index, int low, int high, int l, int r)
    {
        if (high < l || low > r)
        {
            return T();
        }
        if (l <= low && high <= r)
        {
            return segment[index];
        }
        int mid = low + (high - low) / 2;
        T left = query(2 * index + 1, low, mid, l, r);
        T right = query(2 * index + 2, mid + 1, high, l, r);
        return left + right;
    }

    void update(int index, int low, int high, int pos, int x)
    {
        if (low == high)
        {
            segment[index] = T(x);
            return;
        }
        ll mid = low + (high - low) / 2;
        if (pos <= mid)
        {
            update(2 * index + 1, low, mid, pos, x);
        }
        else
        {
            update(2 * index + 2, mid + 1, high, pos, x);
        }
        segment[index] = segment[2 * index + 1] + segment[2 * index + 2];
    }
};

class SparseTable
{
public:
    ll n, p;
    vector<ll> log_table;
    vector<vector<ll>> sparse_table;
    SparseTable(vector<ll> &a, ll n)
    {
        this->n = n;
        p = floor(log2(n));

        log_table.assign(n + 1, 0);
        for (int i = 2; i <= n; i++)
            log_table[i] = log_table[i / 2] + 1;

        sparse_table.assign(p + 1, vector<ll>(n, 0));
        sparse_table[0] = a;
        for (int i = 1; i <= p; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if ((j + (1ll << i)) > n)
                    break;
                ll first_half = sparse_table[i - 1][j];
                ll second_half = sparse_table[i - 1][j + (1ll << (i - 1))];
                sparse_table[i][j] = min(first_half, second_half);
            }
        }
    }

    ll query_overalapping(ll l, ll r)
    {
        ll width = r - l + 1;
        ll k = log_table[width];
        ll left_ans = sparse_table[k][l];
        ll right_ans = sparse_table[k][r - (1ll << k) + 1];
        return min(left_ans, right_ans);
    }

    ll query_non_overlapping(ll l, ll r)
    {
        ll min_value = INF;
        for (int i = log_table[r - l + 1]; l <= r; i = log_table[r - l + 1])
        {
            min_value = min(min_value, sparse_table[i][l]);
            l += (1ll << i);
        }
        return min_value;
    }
};

struct LNode
{
    ll m, s;
    LNode() : s(0), m(-Constants::INF) {};
    LNode(ll m) : s(m), m(m) {};

    LNode operator+(const LNode &other)
    {
        LNode res;
        res.m = max(m, other.m);
        res.s = s + other.s;
        return res;
    }
};

template <typename T>
class SGTreeLazy
{
public:
    vector<T> seg;
    vector<ll> a, lazy;
    ll n;
    vector<bool> cLazy;

    void build(ll ind, ll low, ll high)
    {
        if (low == high)
        {
            seg[ind] = T(a[low]);
            return;
        }
        ll mid = low + (high - low) / 2;
        build(2 * ind + 1, low, mid);
        build(2 * ind + 2, mid + 1, high);
        seg[ind] = seg[2 * ind + 1] + seg[2 * ind + 2];
    }

    void lazyHelper(ll ind, ll low, ll high)
    {
        seg[ind].s += (high - low + 1) * lazy[ind];
        seg[ind].m += lazy[ind];
        if (low != high)
        {
            ll mid = low + (high - low) / 2;
            cLazy[2 * ind + 1] = cLazy[2 * ind + 2] = 1;
            lazy[2 * ind + 1] += lazy[ind];
            lazy[2 * ind + 2] += lazy[ind];
        }
        lazy[ind] = 0;
        cLazy[ind] = 0;
    }

    T queryUtil(ll ind, ll low, ll high, ll l, ll r)
    {
        if (cLazy[ind])
            lazyHelper(ind, low, high);
        if (high < l || r < low)
            return T();
        if (low >= l && high <= r)
            return seg[ind];

        ll mid = low + (high - low) / 2;
        auto left = queryUtil(2 * ind + 1, low, mid, l, r);
        auto right = queryUtil(2 * ind + 2, mid + 1, high, l, r);
        return left + right;
    }

    void updateUtil(ll ind, ll low, ll high, ll l, ll r, ll val)
    {
        if (cLazy[ind])
            lazyHelper(ind, low, high);
        if (high < l || low > r)
            return;
        if (l <= low && high <= r)
        {
            cLazy[ind] = 1;
            lazy[ind] = val;
            lazyHelper(ind, low, high);
            return;
        }
        ll mid = low + (high - low) / 2;
        updateUtil(2 * ind + 1, low, mid, l, r, val);
        updateUtil(2 * ind + 2, mid + 1, high, l, r, val);
        seg[ind] = seg[2 * ind + 1] + seg[2 * ind + 2];
    }

public:
    SGTreeLazy(vector<ll> &a, ll n)
    {
        this->a = a;
        this->n = n;
        seg.assign(4 * n, T());
        lazy.assign(4 * n, 0);
        cLazy.assign(4 * n, 0);
        build(0, 0, n - 1);
    }

    T query(ll l, ll r)
    {
        return queryUtil(0, 0, n - 1, l, r);
    }

    T pQuery(ll pos)
    {
        return queryUtil(0, 0, n - 1, pos, pos);
    }

    void update(ll l, ll r, ll val)
    {
        updateUtil(0, 0, n - 1, l, r, val);
    }

    void pUpdate(ll pos, ll val)
    {
        updateUtil(0, 0, n - 1, pos, pos, val);
    }
};

/*STRUCTURE FOR OFFLINE QUERIES*/
/*NOTE:
-idx is used to obtain the answers in the original order of queries*/
ll rootn;
struct Query
{
    ll l, r, idx;
    Query() : l(0), r(0), idx(0) {};
    Query(ll l, ll r, ll idx) : l(l), r(r), idx(idx) {};
    bool operator<(Query other) const
    {
        return make_pair(l / rootn, r) <
               make_pair(other.l / rootn, other.r);
    }
};

/*SQUARE-ROOT DECOMPOSITION USING MO'S ALGORITHM*/
/*Time Complexity for processing Q queries: O((N+Q)sqrt(N))*/
/*0-BASED INDEXING IS USED FOR l,r,idx*/
/*The interval [curr_l, curr_r] is the current interval (both the end points included)*/
void moAlgorithm(vector<ll> &a, vector<Query> &queries)
{
    ll n = a.size();
    ll q = queries.size();
    ll MAXI = 2e5 + 5;
    rootn = ceil(sqrt(n + 0.1));

    /*Sort the Queries by Blocks*/
    sort(all(queries));

    ll currCount = 0;
    vector<ll> freq(2 * (*max_element(all(a))) + 2, 0ll);

    vector<ll> answer(q);
    ll curr_l = 0, curr_r = -1, curr_ans = 0, l, r; // Do not change this

    auto insert = [&](ll x)
    {
        freq[x]++;
    };

    auto remove = [&](ll x)
    {
        freq[x]--;
    };

    for (int i = 0; i < q; i++)
    {
        l = queries[i].l, r = queries[i].r;
        while (curr_r < r)
        {
            curr_r++;
            insert(a[curr_r]);
        }

        while (curr_l > l)
        {
            curr_l--;
            insert(a[curr_l]);
        }

        while (curr_l < l)
        {
            remove(a[curr_l]);
            curr_l++;
        }

        while (curr_r > r)
        {
            remove(a[curr_r]);
            curr_r--;
        }
        answer[queries[i].idx] = currCount;
    }
    for (auto &ans : answer)
        cout << ans << endl;
}