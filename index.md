<head>
    <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
    <script type="text/x-mathjax-config">
        MathJax.Hub.Config({
            tex2jax: {
            skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
            inlineMath: [['$','$']]
            }
        });
    </script>
</head>

下面是本次研究性学习的所有代码。

### 使用方法

+ 在本地下载一个 IDE（如 dev-c++, Vim, VSCode）等。

+ 编译程序后即可运行。

+ 如果懒得下载 IDE，可以使用在线 IDE。如洛谷的在线 IDE (https://www.luogu.com.cn/ide)，将碳原子数输入到输入框中，输出框会返回答案。

### 代码

所有的答案均对 $998244353$ 取模。如果想更改模数可以找到 `const LL MOD = xxx` 的部分，把 `xxx` 改为您想要的模数，需要注意的是模数必须是一个质数。注意第三份代码不可以更改模数，否则会导致答案错误。

#### 烷基计数代码

时间复杂度 $O(n^3)$，可处理 $n \leq 500$ 的数据。

```cpp
#include <iostream>
#include <cstring>
#include <cstdio>

#define Maxn 510
#define LL long long

using namespace std;

const LL MOD = 998244353;

inline int read() {
    int x = 0, f = 1;
    char c = getchar();
    while(c < '0' || c > '9') {
        if(c == '-') f = -1;
        c = getchar();
    }
    while('0' <= c && c <= '9') {
        x = x * 10 + c - '0';
        c = getchar();
    }
    return x * f;
}

LL f[Maxn];

LL Pow(LL a, LL b) {
    LL ans = 1, base = a;
    while(b) {
        if(b & 1) ans = ans * base % MOD;
        base = base * base % MOD;
        b >>= 1;
    }
    return ans;
}

LL inv2 = Pow(2, MOD - 2), inv6 = Pow(6, MOD - 2);

int main() {
    int n = read();
    f[0] = 1, f[1] = 1;
    for(int i = 2; i <= n; ++i) {
        for(int j = 0; j <= n; ++j)
            for(int k = j; k + j <= n; ++k) {
                int l = i - 1 - j - k; if(l < k) break;
                if(j < k && k < l) f[i] += f[j] * f[k] % MOD * f[l] % MOD;
                else if(j != k && k == l) f[i] += f[j] * f[k] % MOD * (f[k] + 1) % MOD * inv2 % MOD;
                else if(j == k && k != l) f[i] += f[l] * f[k] % MOD * (f[k] + 1) % MOD * inv2 % MOD;
                else if(j == k && k == l) f[i] += (f[j] + 2) * (f[j] + 1) % MOD * f[j] % MOD * inv6 % MOD;
                f[i] %= MOD;
            }
    }
    cout << (f[n] % MOD + MOD) % MOD << endl;
    return 0;
}
```

#### 烷烃计数代码

时间复杂度 $O(n^3)$，可处理 $n \leq 500$ 的数据。

```cpp
#include <iostream>
#include <cstring>
#include <cstdio>

#define Maxn 510
#define LL long long

using namespace std;

const LL MOD = 998244353;

inline int read() {
    int x = 0, f = 1;
    char c = getchar();
    while(c < '0' || c > '9') {
        if(c == '-') f = -1;
        c = getchar();
    }
    while('0' <= c && c <= '9') {
        x = x * 10 + c - '0';
        c = getchar();
    }
    return x * f;
}

LL n, m, f[Maxn];

LL Pow(LL a, LL b) {
    LL ans = 1, base = a;
    while(b) {
        if(b & 1) ans = ans * base % MOD;
        base = base * base % MOD;
        b >>= 1;
    }
    return ans;
}

LL inv2 = Pow(2, MOD - 2), inv6 = Pow(6, MOD - 2), inv24 = Pow(24, MOD - 2);

LL calc() {
	LL ans = 0;
	for(int i = 0; i <= m; ++i)
		for(int j = i; j <= m; ++j)
			for(int k = j; k <= m; ++k) {
				int l = n - 1 - i - j - k; if(l > m || l < k) continue;
				if(i < j  && j < k  && k < l ) ans += f[i] * f[j] % MOD * f[k] % MOD * f[l] % MOD;
				if(i == j && j < k  && k < l ) ans += f[i] * (f[i] + 1) % MOD * inv2 % MOD * f[k] % MOD * f[l] % MOD;
				if(i < j  && j == k && k < l ) ans += f[j] * (f[j] + 1) % MOD * inv2 % MOD * f[i] % MOD * f[l] % MOD;
				if(i < j  && j < k  && k == l) ans += f[k] * (f[k] + 1) % MOD * inv2 % MOD * f[i] % MOD * f[j] % MOD;
				if(i == j && j == k && k < l ) ans += f[i] * (f[i] + 1) % MOD * (f[i] + 2) % MOD * inv6 % MOD * f[l] % MOD;
				if(i < j  && j == k && k == l) ans += f[l] * (f[l] + 1) % MOD * (f[l] + 2) % MOD * inv6 % MOD * f[i] % MOD; 
				if(i == j && j <  k && k == l) ans += f[i] * (f[i] + 1) % MOD * inv2 % MOD * f[k] % MOD * (f[k] + 1) % MOD * inv2 % MOD; 
				if(i == j && j == k && k == l) ans += f[i] * (f[i] + 1) % MOD * (f[i] + 2) % MOD * (f[i] + 3) % MOD * inv24 % MOD;
				ans %= MOD;
			}
	return ans;
}

int main() {
    n = read(), m = (n - 1) / 2;
    f[0] = 1, f[1] = 1;
    for(int i = 2; i <= n / 2; ++i) {
        for(int j = 0; j <= n / 2; ++j)
            for(int k = j; k + j <= n / 2; ++k) { // j <= k <= l
                int l = i - 1 - j - k; if(l < k) break;
                if(j < k && k < l) f[i] += f[j] * f[k] % MOD * f[l] % MOD;
                else if(j != k && k == l) f[i] += f[j] * f[k] % MOD * (f[k] + 1) % MOD * inv2 % MOD;
                else if(j == k && k != l) f[i] += f[l] * f[k] % MOD * (f[k] + 1) % MOD * inv2 % MOD;
                else if(j == k && k == l) f[i] += (f[j] + 2) * (f[j] + 1) % MOD * f[j] % MOD * inv6 % MOD;
                f[i] %= MOD;
            }
    }
    LL ans = calc();
    if(!(n & 1)) ans += f[n / 2] * (f[n / 2] + 1) % MOD * inv2 % MOD;
    cout << ans % MOD << endl;
    return 0;
}
```

#### 烷基计数优化代码

时间复杂度 $O(n \log n)$，可处理 $n \leq 10^5$ 的数据。

```cpp

#pragma GCC optimize(2)
#include <iostream>
#include <cstring>
#include <cstdio>

#define Maxn 800005
#define rg register
#define LL long long

using namespace std;

const LL MOD = 998244353, G = 3;

inline int read() {
    int x = 0, f = 1;
    char c = getchar();
    while(c < '0' || c > '9') {
        if(c == '-') f = -1;
        c = getchar();
    }
    while('0' <= c && c <= '9') {
        x = x * 10 + c - '0';
        c = getchar();
    }
    return x * f;
}

inline LL Pow(LL a, LL b) {
    LL ans = 1, base = a;
    while(b) {
        if(b & 1) ans = ans * base % MOD;
        base = base * base % MOD;
        b >>= 1;
    }
    return ans;
}
const LL invG = Pow(G, MOD - 2);

LL m, trans[Maxn], f[Maxn], inv[Maxn];

struct Poly {
    LL pwI[Maxn], pwG[Maxn];
    void Ready() {
        for(int len = 2; len < Maxn; len <<= 1) {
            pwI[len] = Pow(invG, (MOD - 1) / len);
            pwG[len] = Pow(G, (MOD - 1) / len);
            inv[len] = Pow(len, MOD - 2);
        }
    }
    int lstn = 0;
    inline void NTT(LL *F, bool op, int n) {
        if(n != lstn) for(int i = 0; i < n; ++i) trans[i] = (trans[i >> 1] >> 1) | ((i & 1) ? (n >> 1) : 0);
        for(rg int i = 0; i < n; ++i) if(i < trans[i]) swap(F[i], F[trans[i]]);
        for(int len = 2; len <= n; len <<= 1) {
            LL base = op ? pwI[len] : pwG[len];
            for(int i = 0; i < n; i += len) {
                LL now = 1;
                for(int k = i; k < i + (len >> 1); ++k) {
                    LL temp = F[k + (len >> 1)] * now;
                    F[k + (len >> 1)] = (F[k] - temp) % MOD;
                    F[k] = (F[k] + temp) % MOD;
                    now = now * base % MOD;
                }
            }
        }
        lstn = n;
    }
    LL t1[Maxn], t2[Maxn];
    inline void Getinv(LL *F, LL *invF, int n) {
        // Get the inv of F, record it in invF
        for(int i = 0; i < n; ++i) invF[i] = t1[i] = t2[i] = 0;
        invF[0] = Pow(F[0], MOD - 2);
        for(int len = 2; len <= n; len <<= 1) {
            for(int i = 0; i < (len >> 1); ++i) t1[i] = 2 * invF[i];
            for(int i = 0; i < len; ++i) t2[i] = F[i];
            NTT(invF, 0, len << 1); NTT(t2, 0, len << 1);
            for(int i = 0; i < (len << 1); ++i) invF[i] = invF[i] * invF[i] % MOD * t2[i];
            NTT(invF, 1, len << 1);
            for(int i = 0; i < len; ++i) invF[i] = (t1[i] - invF[i] * inv[len << 1]) % MOD;
            for(int i = len; i < (len << 1); ++i) invF[i] = 0;
        }
    }
}poly;

LL a[Maxn], b[Maxn], c[Maxn], d[Maxn], e[Maxn], inve[Maxn];
inline void init(int len) {
    for(int i = 0; i < len; ++i) a[i] = b[i] = c[i] = d[i] = e[i] = 0;
}

void newton_method() {
    f[0] = 1;
    int n = 1; while(n <= m) n <<= 1;
    for(int len = 2; len <= n; len <<= 1) {
        init(len);
        for(int i = 0; i < (len >> 1); ++i) if(i * 2 < len) a[i * 2] = f[i]; else break;
        for(int i = 0; i < (len >> 1); ++i) if(i * 3 < len) b[i * 3] = f[i]; else break;
        for(int i = 0; i < (len >> 1); ++i) c[i] = f[i];
        poly.NTT(a, 0, len << 1); poly.NTT(b, 0, len << 1); poly.NTT(c, 0, len << 1);
        for(int i = 0; i < (len << 1); ++i) {
            d[i] = c[i] * c[i] % MOD * c[i] + 3 * a[i] * c[i] + 2 * b[i];
            e[i] = 3 * c[i] * c[i] + 3 * a[i];
        }
        poly.NTT(d, 1, len << 1); poly.NTT(e, 1, len << 1);
        for(int i = len; i < len * 1.5; ++i) d[i] = e[i] = 0;
        for(int i = 0; i < len; ++i) d[i] = d[i] * inv[len << 1], e[i] = e[i] * inv[len << 1];
        for(int i = len - 1; i; --i) d[i] = d[i - 1] - 6 * f[i], e[i] = e[i - 1];
        d[0] = -6 * f[0] + 6; e[0] = -6;
        poly.Getinv(e, inve, len);
        poly.NTT(d, 0, len << 1); poly.NTT(inve, 0, len << 1);
        for(int i = 0; i < (len << 1); ++i) d[i] = d[i] * inve[i] % MOD;
        poly.NTT(d, 1, len << 1);
        for(int i = 0; i < len; ++i) f[i] = f[i] - d[i] * inv[len << 1];
    }
}
int main() {
    m = read();
    poly.Ready(); newton_method();
    cout << (f[m] % MOD + MOD) % MOD << endl;
    return 0;
}
```
