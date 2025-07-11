#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;
using namespace std;
using ll = long long;
using ull = unsigned long long;
using ld = long double;
template <class T>
using oset = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T>
using omset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

#define all(v) (v).begin(), (v).end()
#define SpeedIsNeeded                 \
    ios_base::sync_with_stdio(false); \
    cin.tie(NULL);                    \
    cout.tie(NULL)

#ifndef ONLINE_JUDGE
#define debug(x)                           \
    cerr << __LINE__ << ": " << #x << " "; \
    _print(x);                             \
    cerr << endl;
#define custom_solve(t, temp)                   \
    cerr << "Test " << temp - t << ":" << endl; \
    solve();                                    \
    cerr << endl;
#define tcn "Case #" << testCaseNumber++ << ": "
#include "debug.h"
#else
#define custom_solve(t, temp) solve();
#define debug(x)
#endif

constexpr int N = 2e5 + 10;
constexpr ll MOD1 = 1000000007;
constexpr ll MOD2 = 998244353;
constexpr ll INF = 1e18;
constexpr double eps = 1e-6;

template <typename T>
istream &operator>>(istream &in, vector<T> &v)
{
    for (auto &x : v)
        in >> x;
    return in;
}

template <typename T>
ostream &operator<<(ostream &out, vector<T> &v)
{
    for (auto &x : v)
        out << x << " ";
    return out;
}

vector<ll> fact;
ld binPow(ld a, ll b);
ll binPow(ll a, ll b);
void factorials();
void printD(ld a, int prec, bool bankerRound = false); // Prints rounded a to prec decimal places.
ll ncr(ll n, ll r);
ll add(ll a, ll b) { return (a + b) % MOD1; }
ll sub(ll a, ll b) { return (a - b + MOD1) % MOD1; }
ll mul(ll a, ll b) { return (a * 1LL * b) % MOD1; }
ll inv(ll a) { return binPow(a, MOD1 - 2) % MOD1; }
ll divi(ll a, ll b) { return (a * 1LL * inv(b)) % MOD1; }
ll testCaseNumber = 1;

void solve()
{
    ll n;
    cin >> n;
    vector<ll> a(n);
    cin >> a;
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("Input.txt", "r", stdin);
    freopen("Output.txt", "w", stdout);
    freopen("Error.txt", "w", stderr);
#endif
    SpeedIsNeeded;
    int t = 1;
    cin >> t;
    int temp = t;
    while (t--)
    {
        custom_solve(t, temp);
    }
    return 0;
}

ld binPow(ld a, ll b)
{
    ld ans = 1;
    while (b)
    {
        if (b & 1)
            ans = ans * a;
        a = a * a;
        b >>= 1;
    }
    return ans;
}

ll binPow(ll a, ll b)
{
    ld ans = 1;
    while (b)
    {
        if (b & 1)
            ans = mul(ans, a);
        a = mul(a, a);
        b >>= 1;
    }
    return ans;
}

void factorials()
{
    fact.assign(N + 1, 0);
    fact[0] = fact[1] = 1;
    for (int i = 2; i <= N; i++)
        fact[i] = mul(i, fact[i - 1]);
}

ll ncr(ll n, ll r)
{
    ll nr = fact[n];
    ll dr = mul(fact[n - r], fact[r]);
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