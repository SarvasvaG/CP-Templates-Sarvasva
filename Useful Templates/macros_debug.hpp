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

typedef unsigned long long ull;
typedef long double ld;
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

namespace Constants
{
    constexpr ll N = 2e6 + 10;
    constexpr ll MOD1 = 1000000007;
    constexpr ll MOD2 = 998244353;
    constexpr ll MOD3 = ll(1e18) + 3;
    constexpr ll MOD4 = ll(1e15) + 37;
    constexpr ll INF = 1e18;
    constexpr double eps = 1e-6;
};

namespace Math
{
    vector<ll> factorials;
    void computeFactorials();

    ll ncr(ll n, ll r);

    void printD(ld a, int prec, bool bankerRound = false); // Prints rounded a to prec decimal places.
};
using namespace Math;

namespace Mod32
{
    constexpr ll mod = Constants::MOD1;
    ll binPow(ll a, ll b, ll mod = mod);
    inline ll add(ll a, ll b, ll mod = mod) { return (a + b) % mod; }
    inline ll sub(ll a, ll b, ll mod = mod) { return (a - b + mod) % mod; }
    inline ll mul(ll a, ll b, ll mod = mod) { return (a * b) % mod; }
    inline ll inv(ll a, ll mod = mod) { return binPow(a, mod - 2, mod) % mod; }
    inline ll divi(ll a, ll b, ll mod = mod) { return (a * 1LL * inv(b, mod)) % mod; }
};
using namespace Mod32;

namespace Math
{
    void computeFactorials()
    {
        factorials.assign(Constants::N + 1, 0);
        factorials[0] = factorials[1] = 1;
        for (int i = 2; i <= Constants::N; i++)
        {
            factorials[i] = mul(i, factorials[i - 1]);
        }
    }

    ll ncr(ll n, ll r)
    {
        ll nr = factorials[n];
        ll dr = mul(factorials[n - r], factorials[r]);
        return divi(nr, dr);
    }

    void printD(ld a, int prec, bool bankerRound)
    {
        cout << fixed << setprecision(prec);
        ld num = pow(10, prec);
        if (bankerRound)
        {
            ld temp = floor(a * num);
            if ((ll)temp % 2 == 0 && a * num - temp == 0.5)
                cout << temp / num << endl;
            else
                cout << round(a * num) / num << endl;
        }
        else
        {
            cout << round(a * num) / num << endl;
        }
    }
};

namespace Mod32
{
    ll binPow(ll a, ll b, ll mod)
    {
        ll result = 1;
        while (b)
        {
            if (b & 1)
            {
                result = mul(result, a, mod);
            }
            a = mul(a, a, mod);
            b >>= 1;
        }
        return result;
    }
};