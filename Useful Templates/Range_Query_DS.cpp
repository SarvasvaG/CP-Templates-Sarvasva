#include<bits/stdc++.h>
#include "macros_debug.hpp"

using namespace std;

struct Node{
    ll sum;
    Node():sum(0){};
    Node(ll num):sum(num){};
};
 
Node combine(const Node& A, const Node& B){
    Node C(A.sum+B.sum);
    return C;
}
template<typename T>
class SegmentTree{
    public:
    vector<T>  segment;
    vector<ll> array;
    ll arraySize;
 
    SegmentTree(vector<ll> array, ll arraySize):array(array),arraySize(arraySize){
        segment.assign(4*arraySize+1,0);
    }
 
    void build(int index,int low,int high){
        if(low==high){
            segment[index]=T(array[low]);
            return;
        }
        ll mid = low + (high-low)/2;
        build(2*index+1,low,mid);
        build(2*index+2,mid+1,high);
        segment[index]=combine(segment[2*index+1],segment[2*index+2]);
    }
 
    T query(int index, int low, int high, int l, int r){
        if(high<l||low>r){
            return T();
        }
        if(l<=low&&high<=r){
            return segment[index];
        }
        int mid = low + (high-low)/2;
        T left = query(2*index+1,low,mid,l,r);
        T right = query(2*index+2,mid+1,high,l,r);
        return combine(left,right);
    }
    void update(int index,int low,int high,int i,int x){
        if(low==high){
            segment[index]=T(x);
            return;
        }
        int mid=low+(high-low)/2;
        if(i<=mid) update(2*index+1,low,mid,i,x);
        else update(2*index+2,mid+1,high,i,x);
        segment[index]=combine(segment[2*index+1],segment[2*index+2]);
    }
};

class SparseTable{
    public:
    ll n,p;
    vector<ll> log_table;
    vector<vector<ll>> sparse_table;
    SparseTable(vector<ll> &a, ll n){
        this->n=n;
        p=floor(log2(n));
 
        log_table.assign(n+1,0);
        fe(i,2,n+1,1)
            log_table[i]=log_table[i/2]+1;
        
        sparse_table.assign(p+1,vector<ll>(n,0));
        sparse_table[0]=a;
        fe(i,1,p+1,1){
            fe(j,0,n,1){
                if((j+(1ll<<i))>n)
                    break;
                ll first_half=sparse_table[i-1][j];
                ll second_half=sparse_table[i-1][j+(1ll<<(i-1))];
                sparse_table[i][j]=min(first_half,second_half);
            }
        }
    }
 
    ll query_overalapping(ll l, ll r){
        ll width=r-l+1;
        ll k=log_table[width];
        ll left_ans=sparse_table[k][l];
        ll right_ans=sparse_table[k][r-(1ll<<k)+1];
        return min(left_ans,right_ans);
    }
 
    ll query_non_overlapping(ll l, ll r){
        ll min_value=INF;
        for(int i=log_table[r-l+1];l<=r;i=log_table[r-l+1]){
            min_value=min(min_value,sparse_table[i][l]);
            l+=(1ll<<i);
        }
        return min_value;
    }
};


class SGTreeLazy{
    public:
    vector<ll> seg,a,lazy;ll n;
    vector<bool> reset;
    public:
    SGTreeLazy(vector<ll> &a, ll n){
        this->a=a;
        this->n=n;
        seg.assign(4*n,0);
        lazy.assign(4*n,0);
        reset.assign(4*n,0);
    }
    void build(ll ind, ll low, ll high){
        if(low==high){
            seg[ind]=a[low];
            return;
        }
        ll mid = low+(high-low)/2;
        build(2*ind+1,low,mid);
        build(2*ind+2,mid+1,high);
        seg[ind]=seg[2*ind+1]+seg[2*ind+2];
    }
    void lazy_helper(ll ind,ll low,ll high){
        seg[ind]+=(high-low+1)*lazy[ind];
        if(low!=high){
            ll mid=low+(high-low)/2;
            lazy[2*ind+1]+=lazy[ind];
            lazy[2*ind+2]+=lazy[ind];
        }
        lazy[ind]=0;
    }
    ll query(ll ind, ll low, ll high, ll l, ll r){
        if(lazy[ind]!=0)
            lazy_helper(ind,low,high);      
        if(high<l||r<low)
            return 0;
        if(low>=l&&high<=r)
            return seg[ind];
        
        ll mid = low + (high-low)/2;
        ll left = query(2*ind+1, low, mid, l, r);
        ll right = query(2*ind+2, mid+1, high, l, r);
        return left+right;
    }
    void update(ll ind,ll low, ll high, ll l, ll r, ll val){
        if(lazy[ind]!=0)
            lazy_helper(ind,low,high);
        if(high<l||low>r)
            return;
        if(l<=low&&high<=r){
            seg[ind]+=(high-low+1)*val;
            if(low!=high){
                ll mid=low+(high-low)/2;
                lazy[2*ind+1]+=val;
                lazy[2*ind+2]+=val;
            }
            return;
        }
        ll mid = low + (high-low)/2;
        update(2*ind+1,low,mid,l,r,val);
        update(2*ind+2,mid+1,high,l,r,val);
        seg[ind]=seg[2*ind+1]+seg[2*ind+2];
    }
};
 
class SparseTable{
    public:
    ll n,p;
    vector<ll> log_table;
    vector<vector<ll>> sparse_table;
    SparseTable(vector<ll> &a, ll n){
        this->n=n;
        p=floor(log2(n));
 
        log_table.assign(n+1,0);
        fe(i,2,n+1,1)
            log_table[i]=log_table[i/2]+1;
        
        sparse_table.assign(p+1,vector<ll>(n,0));
        sparse_table[0]=a;
        fe(i,1,p+1,1){
            fe(j,0,n,1){
                if((j+(1ll<<i))>n)
                    break;
                ll first_half=sparse_table[i-1][j];
                ll second_half=sparse_table[i-1][j+(1ll<<(i-1))];
                sparse_table[i][j]=min(first_half,second_half);
            }
        }
    }
 
    ll query_overalapping(ll l, ll r){
        ll width=r-l+1;
        ll k=log_table[width];
        ll left_ans=sparse_table[k][l];
        ll right_ans=sparse_table[k][r-(1ll<<k)+1];
        return min(left_ans,right_ans);
    }
 
    ll query_non_overlapping(ll l, ll r){
        ll min_value=INF;
        for(int i=log_table[r-l+1];l<=r;i=log_table[r-l+1]){
            min_value=min(min_value,sparse_table[i][l]);
            l+=(1ll<<i);
        }
        return min_value;
    }
};


/*STRUCTURE FOR OFFLINE QUERIES*/
/*NOTE:
-idx is used to obtain the answers in the original order of queries*/
ll rootn;
struct Query{
    ll l,r,idx;
    Query():l(0),r(0),idx(0){};
    Query(ll l,ll r,ll idx):l(l),r(r),idx(idx){};
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
 
void moAlgorithm(vector<int> &a, vector<Query> &queries){
    int n=a.size();
    int q=queries.size();
    int MAXI=2e5+5;
    rootn=ceil(sqrt(n+0.1));
    
    /*Sort the Queries by Blocks*/
    sort(all(queries));
 
    unordered_map<int,int> freq;
    vector<int> answer(q);
    set<pair<int,int>> st;
    ll curr_l=0,curr_r=-1,curr_ans=0,l,r;// Do not change this
 
    auto insert = [&](ll x){
        freq[x]++;
        if(freq[x]!=1) 
            st.erase({freq[x]-1,x});
        st.insert({freq[x],x});
    };
 
    auto remove = [&](ll x){
        freq[x]--;
        st.erase({freq[x]+1,x});
        if(freq[x]!=0)
            st.insert({freq[x],x});
        if(!freq[x]) freq.erase(x);
    };
 
    for(int i=0;i<q;i++){
        l=queries[i].l, r=queries[i].r;
        while(curr_r<r){
            curr_r++;
            insert(a[curr_r]);
        }
 
        while(curr_l>l){
            curr_l--;
            insert(a[curr_l]);
        }
 
        while(curr_l<l){
            remove(a[curr_l]);
            curr_l++;
        }
 
        while(curr_r>r){
            remove(a[curr_r]);
            curr_r--;
        }
 
        debug(st);
        debug(queries[i].idx);
        answer[queries[i].idx]=(r-l+1)-(st.rbegin()->first);
    }
    for(auto &ans : answer)
        cout<<ans<<endl;
}
 