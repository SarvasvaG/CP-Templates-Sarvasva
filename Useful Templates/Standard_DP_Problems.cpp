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