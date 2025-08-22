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
#include "debug.h"
#else
#define custom_solve(t, temp) solve();
#define debug(x)
#endif

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

namespace Constants
{
    constexpr ll N = 2e6 + 10;
    constexpr ll N2 = 5005;
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
    void computeFactorials(ll n = Constants::N);

    vector<vector<ll>> g;
    void precomputeGCD(ll n = Constants::N2);

    ll ncr(ll n, ll r);

    void printD(ld a, int prec, bool bankerRound = false); // Prints rounded a to prec decimal places.

    namespace Matrix
    {
        vector<vector<ll>> matAdd(const vector<vector<ll>> &a, const vector<vector<ll>> &b);
        vector<vector<ll>> matSub(const vector<vector<ll>> &a, const vector<vector<ll>> &b);
        vector<vector<ll>> matMul(const vector<vector<ll>> &a, const vector<vector<ll>> &b);
        vector<vector<ll>> matPow(vector<vector<ll>> a, ll k);
    };
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

void solve()
{
    ll n;
    cin >> n;
    vector<ll> a(n);
    cin >> a;
    ll gd = 0;
    for (auto x : a)
        gd = g[x][gd];
    ll cnt = 0;
    for (auto &x : a)
    {
        x /= gd;
        if (x == 1)
            cnt++;
    }
    if (cnt >= 1)
    {
        cout << n - cnt << endl;
        return;
    }

    ll s = *max_element(all(a));
    ll inf = Constants::INF;
    vector<ll> dp(s + 1, inf);
    dp[1] = 0;
    for (int i = 2; i <= s; i++)
    {
        for (int j = 0; j < n; j++)
        {
            dp[i] = min(dp[i], 1 + dp[g[i][a[j]]]);
        }
    }
    debug(a);
    debug(dp);
    ll mini = inf;
    for (auto x : a)
        mini = min(mini, dp[x]);
    cout << mini + n - 1 << endl;
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
    precomputeGCD();
    int temp = t;
    while (t--)
    {
        custom_solve(t, temp);
    }
    return 0;
}

namespace Math
{
    void computeFactorials(ll n)
    {
        factorials.assign(n + 1, 0);
        factorials[0] = factorials[1] = 1;
        for (int i = 2; i <= n; i++)
        {
            factorials[i] = mul(i, factorials[i - 1]);
        }
    }

    void precomputeGCD(ll n)
    {
        g.assign(n + 1, vector<ll>(n + 1, 0));
        for (int i = 1; i <= n; i++)
        {
            g[0][i] = g[i][0] = i;
        }
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                g[i][j] = g[j][i % j];
                g[j][i] = g[i][j];
            }
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

    namespace Matrix
    {
        vector<vector<ll>> matAdd(const vector<vector<ll>> &a, const vector<vector<ll>> &b)
        {
            ll n1 = a.size(), n2 = b.size(), m1 = a[0].size(), m2 = b[0].size();
            if (n1 != m1 || n2 != m2)
                return {{}};
            vector<vector<ll>> c(n1, vector<ll>(m1, 0));
            for (ll i = 0; i < n1; i++)
            {
                for (ll j = 0; j < m1; j++)
                {
                    c[i][j] = add(a[i][j], b[i][j]);
                }
            }
            return c;
        }

        vector<vector<ll>> matSub(const vector<vector<ll>> &a, const vector<vector<ll>> &b)
        {
            ll n1 = a.size(), n2 = b.size(), m1 = a[0].size(), m2 = b[0].size();
            if (n1 != m1 || n2 != m2)
            {
                return {{}};
            }
            vector<vector<ll>> c(n1, vector<ll>(m1, 0));
            for (ll i = 0; i < n1; i++)
            {
                for (ll j = 0; j < m1; j++)
                {
                    c[i][j] = sub(a[i][j], b[i][j]);
                }
            }
            return c;
        }

        vector<vector<ll>> matMul(const vector<vector<ll>> &a, const vector<vector<ll>> &b)
        {
            ll n1 = a.size(), n2 = b.size(), m1 = a[0].size(), m2 = b[0].size();
            if (n2 != m1)
            {
                return {{}};
            }
            vector<vector<ll>> c(n1, vector<ll>(m2, 0));
            for (ll i = 0; i < n1; i++)
            {
                for (ll j = 0; j < m2; j++)
                {
                    for (ll k = 0; k < n2; k++)
                    {
                        c[i][j] = add(c[i][j], mul(a[i][k], b[k][j]));
                    }
                }
            }
            return c;
        }

        vector<vector<ll>> matPow(vector<vector<ll>> a, ll k)
        {
            ll n = a.size();
            vector<vector<ll>> ans(n, vector<ll>(n, 0));
            for (ll i = 0; i < n; i++)
                ans[i][i] = 1;
            while (k)
            {
                if (k & 1)
                {
                    ans = matMul(ans, a);
                }
                a = matMul(a, a);
                k >>= 1;
            }
            return ans;
        }
    };
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