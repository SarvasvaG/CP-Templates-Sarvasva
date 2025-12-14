#include "macros_debug.hpp"

struct TrieNode
{
    TrieNode *child[2];
    int cnt;
    TrieNode()
    {
        child[0] = child[1] = nullptr;
        cnt = 0;
    }
};

struct BinaryTrie
{
    const int BITS = 31;
    TrieNode *root;
    long long total = 0;

    BinaryTrie()
    {
        root = new TrieNode();
    }

    void insert(int x)
    {
        TrieNode *cur = root;
        for (int b = BITS; b >= 0; b--)
        {
            int bit = (x >> b) & 1;
            if (!cur->child[bit])
                cur->child[bit] = new TrieNode();
            cur = cur->child[bit];
            cur->cnt++;
        }
        total++;
    }

    void erase(int x)
    {
        TrieNode *cur = root;
        for (int b = BITS; b >= 0; b--)
        {
            int bit = (x >> b) & 1;
            cur = cur->child[bit];
            cur->cnt--;
        }
        total--;
    }

    /*Returns Count of elements in the trie whose xor with x is < K*/
    long long countLess(int x, int K)
    {
        TrieNode *curr = root;
        long long res = 0;
        for (int i = BITS; i >= 0; i--)
        {
            if (!curr)
                break;
            int x_bit = (x >> i) & 1;
            int k_bit = (K >> i) & 1;

            if (x_bit == k_bit)
            {
                if (k_bit && curr->child[1])
                    res += curr->child[1]->cnt;
                curr = curr->child[0];
            }
            else
            {
                if (k_bit && curr->child[0])
                    res += curr->child[0]->cnt;
                curr = curr->child[1];
            }
        }
        return res;
    }

    /*Returns Count of elements in the trie whose xor with x is <= K*/
    long long countLessOrEqual(int x, int K)
    {
        return countLess(x, K + 1);
    }

    /*Returns Count of elements in the trie whose xor with x is > K*/
    long long countMore(int x, int K)
    {
        total - countLessOrEqual(x, K);
    }

    /*Returns max xor between n and any of the numbers in Trie*/
    int findMaxXOR(int n)
    {
        TrieNode *curr = root;
        int max_xor = 0;
        for (int i = BITS; i >= 0; i--)
        {
            int bit = (n >> i) & 1;
            if (curr->child[bit ^ 1] && curr->child[bit ^ 1]->cnt > 0)
            {
                max_xor |= (1 << i);
                curr = curr->child[bit ^ 1];
            }
            else
            {
                curr = curr->child[bit];
            }
        }
        return max_xor;
    }
};
