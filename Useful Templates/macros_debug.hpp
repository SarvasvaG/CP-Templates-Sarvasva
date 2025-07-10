#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
template <class T>
using oset = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T>
using omset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

#define fastIO()                      \
    ios_base::sync_with_stdio(false); \
    cin.tie(NULL);                    \
    cout.tie(NULL)

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef vector<vector<int>> vvi;
typedef vector<vector<ll>> vvll;
typedef vector<pair<int, int>> vpii;
typedef vector<pair<ll, ll>> vpll;
typedef vector<string> vs;

#define nl "\n"
#define setp(x, n) fixed << setprecision(n) << x
#define bit_count __builtin_popcountll
#define alice cout << "Alice" << nl
#define bob cout << "Bob" << nl
#define yes cout << "YES" << nl
#define no cout << "NO" << nl

#define fe(i, a, b, inc) for (int i = a; i < b; i += inc)
#define fer(i, a, b, dec) for (int i = a; i >= b; i -= dec)

#define all(v) (v).begin(), (v).end()
#define sz(v) ((int)(v).size())
#define lbs(x) lower_bound(x) // Use this for sets and maps
#define ubs(x) upper_bound(x) // Use this for sets and maps
#define lb(v, x) lower_bound(all(v), (x))
#define ub(v, x) upper_bound(all(v), (x))

#define pb push_back
#define ppb pop_back
#define mp make_pair
#define ff first
#define ss second

const int N = 2e5 + 10;
const ll MOD1 = 1000000007;
const ll MOD2 = 998244353;
const ll INF = 1e18;
const double eps = 1e-6;
const double PI = 3.14159265358979;

#ifndef ONLINE_JUDGE
#define debug(x)                           \
    cerr << __LINE__ << ": " << #x << " "; \
    _print(x);                             \
    cerr << nl;
#define custom_solve(t, temp)                 \
    cerr << "Test " << temp - t << ":" << nl; \
    solve();                                  \
    cerr << nl;
#else
#define custom_solve(t, temp) solve();
#define debug(x)
#endif

void _print(int x) { cerr << x; }
void _print(ll x) { cerr << x; }
void _print(ull x) { cerr << x; }
void _print(string x) { cerr << x; }
void _print(char x) { cerr << x; }
void _print(double x) { cerr << x; }

template <class T, class V>
void _print(pair<T, V> p);
template <class T>
void _print(vector<T> v);
template <class T>
void _print(set<T> v);
template <class T, class V>
void _print(map<T, V> v);
template <class T>
void _print(multiset<T> v);

template <class T, class V>
void _print(pair<T, V> p)
{
    cerr << "{";
    _print(p.ff);
    cerr << ",";
    _print(p.ss);
    cerr << "}";
}
template <class T>
void _print(vector<T> v)
{
    cerr << "[ ";
    for (T i : v)
    {
        _print(i);
        cerr << " ";
    }
    cerr << "]";
}
template <class T>
void _print(set<T> v)
{
    cerr << "[ ";
    for (T i : v)
    {
        _print(i);
        cerr << " ";
    }
    cerr << "]";
}
template <class T>
void _print(multiset<T> v)
{
    cerr << "[ ";
    for (T i : v)
    {
        _print(i);
        cerr << " ";
    }
    cerr << "]";
}
template <class T, class V>
void _print(map<T, V> v)
{
    cerr << "[ ";
    for (auto i : v)
    {
        _print(i);
        cerr << " ";
    }
    cerr << "]";
}

ll add(ll a, ll b)
{
    return (a + b) % MOD1;
}

ll sub(ll a, ll b)
{
    return (a - b + MOD1) % MOD1;
}

ll mul(ll a, ll b)
{
    return (a * 1LL * b) % MOD1;
}

ll binPow(ll a, ll b)
{
    ll ans = 1;
    while (b){
        if (b & 1) ans = mul(ans, a);
        a = mul(a, a);
        b >>= 1;
    }
    return ans;
}

ll inv(ll a) { return binPow(a, MOD1 - 2) % MOD1; }
ll divi(ll a, ll b) { return (a * 1LL * inv(b)) % MOD1; }