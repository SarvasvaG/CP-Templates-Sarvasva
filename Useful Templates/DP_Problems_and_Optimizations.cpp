#include "macros_debug.hpp"

/*PLEASE FOLLOW 1-BASED INDEXING FOR WEIGHTS AND VALUES
KS[i][w]=Max value obtainable using [1..i] items with total weight of bag exactly equal to w*/
ll knapsack1(ll n, ll W, vector<ll> &weights, vector<ll> &values)
{
    vector<vector<ll>> KS(n + 1, vector<ll>(W + 1, 0));
    for (int i = 1; i <= n; i++)
        for (int w = 0; w <= W; w++)
            if (w - weights[i] >= 0)
                KS[i][w] = max(KS[i - 1][w], values[i] + KS[i - 1][w - weights[i]]);
            else
                KS[i][w] = KS[i - 1][w];
    return *max_element(all(KS[n]));
}

/*PLEASE FOLLOW 1-BASED INDEXING FOR WEIGHTS AND VALUES
KS[i][v]=Minimum weight required for obtaining value exactly equal to v using items [1..i] */
ll knapsack2(ll n, ll W, vector<ll> &weights, vector<ll> &values)
{
    ll V = accumulate(all(values), 0LL);
    vector<vector<ll>> KS(n + 1, vector<ll>(V + 1, INF));
    KS[0][0] = 0;
    for (int i = 1; i <= n; i++)
        for (int v = 0; v <= V; v++)
            if (v - values[i] >= 0)
                KS[i][v] = min(KS[i - 1][v], KS[i - 1][v - values[i]] + weights[i]);
            else
                KS[i][v] = KS[i - 1][v];
    for (int v = V; v >= 1; v--)
        if (KS[n][v] <= W)
            return v;
    return -1;
}

void print_lcs(int i, int j, string &s1, string &s2, vector<vector<ll>> &LCS)
{
    if (i >= s1.size() || j >= s2.size())
        return;

    if (s1[i] == s2[j])
    {
        cout << s1[i];
        print_lcs(i + 1, j + 1, s1, s2, LCS);
    }
    else
    {
        ll a = LCS[i][j], b = LCS[i + 1][j], c = LCS[i][j + 1];
        if (a == b)
            print_lcs(i + 1, j, s1, s2, LCS);
        else
            print_lcs(i, j + 1, s1, s2, LCS);
    }
}

/*0-BASED INDEXING IS USED
LCS[i][j]=Length of LCS of s1[i..n1-1] and s2[j..n2-1]
Returns: LCS[0][0]=LCS of s1 and s2
If printlcs=1, lcs of strings is also printed*/
ll lcs(string &s1, string &s2, bool printlcs)
{
    ll n1 = sz(s1), n2 = sz(s2);
    vector<vector<ll>> LCS;
    LCS.assign(n1 + 1, vector<ll>(n2 + 1, 0));
    for (int i = n1 - 1; i >= 0; i--)
        for (int j = n2 - 1; j >= 0; j--)
            if (s1[i] == s2[j])
                LCS[i][j] = 1 + LCS[i + 1][j + 1];
            else
                LCS[i][j] = max(LCS[i + 1][j], LCS[i][j + 1]);
    if (printlcs)
        print_lcs(0, 0, s1, s2, LCS), cout << nl;
    return LCS[0][0];
}

/*0-BASED INDEXING FOR ARRAY A*/
/*Returns: lis of length n+1
lis[i]=smallest integer which ends an increasing sequence of length i
Note: LIS - Longest strictly increasing subsequence*/
vector<ll> getLIS(const vector<ll> &a)
{
    ll n = a.size();
    ll inf = Constants::INF;
    vector<ll> lis(n + 1, inf);
    lis[0] = 0;
    for (int i = 0; i < n; i++)
    {
        ll it = lower_bound(all(lis), a[i]) - lis.begin();
        lis[it] = min(lis[it], a[i]);
    }
    return lis;
}

/*
Convex Hull Trick: 
1. Adding Line in O(logn) (y=mx+c)
2. Querying min/max value among all lines added given the input x in O(logn)
*/
struct ConvexHullDynamic
{
    static const ll INF = 1e18;

    struct Line
    {
        ll a, b;      // y = ax + b
        double xLeft; // Stores the intersection wiith previous line in the convex hull. First line has -INF

        enum Type
        {
            line,
            maxQuery,
            minQuery
        } type;
        ll val;

        explicit Line(ll aa = 0, ll bb = 0) : a(aa), b(bb), xLeft(-INF), type(Type::line), val(0) {}

        ll valueAt(ll x) const
        {
            return a * x + b;
        }
        friend bool isParallel(const Line &l1, const Line &l2)
        {
            return l1.a == l2.a;
        }
        friend double intersectX(const Line &l1, const Line &l2)
        {
            return isParallel(l1, l2) ? INF : 1.0 * (l2.b - l1.b) / (l1.a - l2.a);
        }
        bool operator<(const Line &l2) const
        {
            if (l2.type == line)
                return this->a < l2.a;
            if (l2.type == maxQuery)
                return this->xLeft < l2.val;
            if (l2.type == minQuery)
                return this->xLeft > l2.val;
        }
    };

    bool isMax;
    set<Line> hull;

    bool hasPrev(set<Line>::iterator it)
    {
        return it != hull.begin();
    }
    bool hasNext(set<Line>::iterator it)
    {
        return it != hull.end() && next(it) != hull.end();
    }
    bool irrelevant(const Line &l1, const Line &l2, const Line &l3)
    {
        return intersectX(l1, l3) <= intersectX(l1, l2);
    }
    bool irrelevant(set<Line>::iterator it)
    {
        return hasPrev(it) && hasNext(it) && ((isMax && irrelevant(*prev(it), *it, *next(it))) || (!isMax && irrelevant(*next(it), *it, *prev(it))));
    }
    // Updates xValue of line pointed by it
    set<Line>::iterator updateLeftBorder(set<Line>::iterator it)
    {
        if (isMax && !hasPrev(it) || !isMax && !hasNext(it))
            return it;
        double val = intersectX(*it, isMax ? (*prev(it)) : (*next(it)));
        Line temp(*it);
        it = hull.erase(it);
        temp.xLeft = val;
        it = hull.insert(it, temp);
        return it;
    }

    explicit ConvexHullDynamic(bool isMax) : isMax(isMax) {}

    void addLine(ll a, ll b) // Add ax + b in logN time
    {
        Line l3 = Line(a, b);
        auto it = hull.lower_bound(l3);

        // If parallel liune is already in set, one of the lines becomes irrelevant
        if (it != hull.end() && isParallel(*it, l3))
        {
            if (isMax && it->b < b || !isMax && it->b > b)
                it = hull.erase(it);
            else
                return;
        }

        it = hull.insert(it, l3);
        if (irrelevant(it))
        {
            hull.erase(it);
            return;
        }

        // Remove lines which became irrelevant after inserting
        while (hasPrev(it) && irrelevant(prev(it)))
            hull.erase(prev(it));
        while (hasNext(it) && irrelevant(next(it)))
            hull.erase(next(it));

        // Update xLine
        it = updateLeftBorder(it);
        if (hasPrev(it))
            updateLeftBorder(prev(it));
        if (hasNext(it))
            updateLeftBorder(next(it));
    }

    ll getBest(ll x)
    {
        Line q;
        q.val = x;
        q.type = isMax ? Line::Type::maxQuery : Line::Type::minQuery;

        auto bestLine = hull.lower_bound(q);
        if (isMax)
            --bestLine;
        return bestLine->valueAt(x);
    }
};

/*
SOS DP
F[mask]=Sum of Array Elements with index i such that i is SUBSET of mask (Assuming initially, F[mask]=A[mask])
Time Complexity: O(MAXBITS.2^MAXBITS)
*/
void forwardProp1(vector<ll> &F, ll MAXBITS = Constants::MAXBITS, ll MAXN = Constants::MAXN){
    for(int i=0;i<=MAXBITS;i++){
        for(int mask=0;mask<MAXN;mask++){
            if(mask&(1ll<<i)){
                F[mask]+=F[mask^(1ll<<i)];
            }
        }
    }
}

/*
SOS DP
Revert the previous process
Time Complexity: O(MAXBITS.2^MAXBITS)
*/
void backwardProp1(vector<ll> &F, ll MAXBITS = Constants::MAXBITS, ll MAXN = Constants::MAXN){
    for(int i=0;i<=MAXBITS;i++){
        for(int mask=MAXN-1;mask>=0;mask--){
            if(mask&(1ll<<i)){
                F[mask]-=F[mask^(1ll<<i)];
            }
        }
    }
}

/*
SOS DP
F[mask]=Sum of Array Elements with index i such that i is SUPERSET of mask (Assuming initially, F[mask]=A[mask])
Time Complexity: O(MAXBITS.2^MAXBITS)
*/
void forwardProp2(vector<ll> &F, ll MAXBITS = Constants::MAXBITS, ll MAXN = Constants::MAXN){
    for(int i=0;i<=MAXBITS;i++){
        for(int mask=MAXN-1;mask>=0;mask--){
            if(mask&(1ll<<i)){
                F[mask^(1ll<<i)]+=F[mask];
            }
        }
    }
}

/*
SOS DP
Revert the Previous Process
Time Complexity: O(MAXBITS.2^MAXBITS)
*/
void backwardProp2(vector<ll> &F, ll MAXBITS = Constants::MAXBITS, ll MAXN = Constants::MAXN){
    for(int i=0;i<=MAXBITS;i++){
        for(int mask=0;mask<MAXN;mask++){
            if(mask&(1ll<<i)){
                F[mask^(1ll<<i)]-=F[mask];
            }
        }
    }
}
