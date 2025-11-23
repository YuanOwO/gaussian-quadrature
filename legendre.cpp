#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using ld = long double;
using Matrix = std::vector<std::vector<ld>>;

const ld PI = acosl(-1);

////////////////////////////////////////////////////////////////////////////////

/**
 * @brief 產生 n 次 Legendre 多項式的伴隨矩陣
 *
 * @param n 多項式次數
 * @return Matrix n * n 對稱三對角矩陣
 */
Matrix legendre_companion_matrix(int n) {
    Matrix J(n, std::vector<ld>(n, 0.0L));

    // 填充 beta_k = k / sqrt(4k^2 - 1)
    for (int k = 1; k < n; k++) {
        ld beta = k / std::sqrt(4.0L * k * k - 1.0L);
        J[k][k - 1] = J[k - 1][k] = beta;
    }

    return J;
}

/**
 * @brief 使用 Jacobi 方法計算對稱矩陣的特徵值
 *
 * Jacobi 方法原理：
 *     - 對稱矩陣 A 透過一系列旋轉（Givens rotation）逐步消掉非對角元素
 *     - 收斂後，矩陣變成對角矩陣，對角線即為特徵值
 *
 * 注意：
 *     - 只返回特徵值，不計算特徵向量
 *     - 適合小到中型對稱矩陣 (n < 1000)
 *
 * @param _A 輸入對稱矩陣
 * @return std::vector<ld> 特徵值
 */
std::vector<ld> jacobi_eigen(const Matrix& _A) {
    int n = _A.size();
    Matrix A = _A;                  // 對矩陣做副本，避免修改原矩陣
    const ld eps = LDBL_EPSILON;    // 收斂閾值
    const int max_iter = 1000 * n;  // 最大迭代次數

    for (int iter = 0; iter < max_iter; iter++) {
        // 1. 找最大非對角元素 a[p][q]
        int p = 0, q = 1;
        ld max_off = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabsl(A[i][j]) > max_off) {
                    max_off = fabsl(A[i][j]);
                    p = i;
                    q = j;
                }

        // 2. 若所有非對角元素都小了 => 結束
        if (max_off < eps) break;

        ld app = A[p][p];
        ld aqq = A[q][q];
        ld apq = A[p][q];

        // 3. 計算 Jacobi rotation 參數
        ld tau = (aqq - app) / (2 * apq);
        ld t = (tau >= 0 ? 1.0L / (tau + sqrtl(1 + tau * tau)) : 1.0L / (tau - sqrtl(1 + tau * tau)));
        ld c = 1 / sqrtl(1 + t * t);
        ld s = c * t;

        // 4. 更新對角元素 (p,p), (q,q)
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = A[q][p] = 0;

        // 5. 更新其他元素
        for (int k = 0; k < n; k++) {
            if (k != p && k != q) {
                ld akp = A[k][p];
                ld akq = A[k][q];
                A[k][p] = A[p][k] = c * akp - s * akq;
                A[k][q] = A[q][k] = c * akq + s * akp;
            }
        }
    }

    // 6. 對角線即為 eigenvalues
    std::vector<ld> eigen(n);
    for (int i = 0; i < n; i++)
        eigen[i] = A[i][i];
    return eigen;
}

/**
 * @brief 計算 Legendre 多項式 P_n(x)
 *
 * 使用遞迴公式：
 *     - P_0(x) = 1
 *     - P_1(x) = x
 *     - P_k(x) = ((2k - 1) x P_{k-1}(x) - (k - 1) P_{k-2}(x)) / k
 *
 * @param n 多項式次數
 * @param x 輸入值 (-1 <= x <= 1)
 * @return ld P_n(x)
 */
ld legendre_value(int n, ld x) {
    // P_0(x) = 1
    if (n == 0) return 1.0L;

    // P_1(x) = x
    if (n == 1) return x;

    ld Pn_2 = 1.0L;  // P_0(x)
    ld Pn_1 = x;     // P_1(x)
    ld Pn;

    // 使用遞迴公式從 P_2(x) 計算到 P_n(x)
    for (int k = 2; k <= n; k++) {
        // P_k = ((2k - 1) x P_{k-1} - (k - 1) P_{k-2}) / k
        Pn = ((2.0L * k - 1.0L) * x * Pn_1 - (k - 1.0L) * Pn_2) / k;
        Pn_2 = Pn_1;
        Pn_1 = Pn;
    }

    return Pn;
}

/**
 * @brief 計算 Legendre 多項式 P_n(x) 的導數 P_n'(x)
 *
 * 使用關係式：
 *     P_n'(x) = n / (x^2 - 1) * (x P_n(x) - P_{n-1}(x))
 *
 * @param n 多項式次數
 * @param x 輸入值 (-1 <= x <= 1)
 * @return ld P_n'(x)
 */
ld legendre_derivative(int n, ld x) {
    // P_0(x) = 1 -> P_0'(x) = 0
    if (n == 0) return 0.0L;

    // P_1(x) = x -> P_1'(x) = 1
    if (n == 1) return 1.0L;

    ld Pn = legendre_value(n, x);
    ld Pn_1 = legendre_value(n - 1, x);

    return n / (x * x - 1.0L) * (x * Pn - Pn_1);
}

/**
 * @brief 計算 n 次 Legendre 多項式的根和權重
 *
 * 使用 Golub–Welsch 方法：
 *   1. 先生成 n 次 Legendre 多項式的 Jacobi（伴隨）矩陣
 *   2. 計算其特徵值作為初始根
 *   3. 可選用一次 Newton 方法修正根（提高精度）
 *   4. 計算對應權重
 *
 * @param n 多項式次數
 * @return std::pair<std::vector<ld>, std::vector<ld>> pair<roots, weights>
 */
std::pair<std::vector<ld>, std::vector<ld>> legendre_roots_weights(int n) {
    // first approximation of roots. We use the fact that the companion
    // matrix is symmetric in this case in order to obtain better zeros.
    Matrix J = legendre_companion_matrix(n);
    std::vector<ld> roots = jacobi_eigen(J);
    std::sort(roots.begin(), roots.end());  // 將 roots 由小排到大

    std::vector<ld> dy(n), df(n), fm(n), weights(n);
    ld max_df = 0, max_fm = 0;

    // improve roots by one application of Newton's method
    for (int i = 0; i < n; i++) {
        dy[i] = legendre_value(n, roots[i]);
        df[i] = legendre_derivative(n, roots[i]);
        roots[i] -= dy[i] / df[i];
        max_df = std::max(max_df, std::fabs(df[i]));
    }

    // compute the weights. We scale the factor to avoid possible numerical
    // overflow.
    for (int i = 0; i < n; i++) {
        fm[i] = legendre_value(n - 1, roots[i]);
        max_fm = std::max(max_fm, std::fabs(fm[i]));
    }

    for (int i = 0; i < n; i++) {
        fm[i] /= max_fm, df[i] /= max_df;
    }

    for (int i = 0; i < n; i++) {
        weights[i] = 1.0L / (fm[i] * df[i]);
    }

    // for Legendre we can also symmetrize the roots and weights
    std::vector<ld> w_copy = weights, r_copy = roots;
    ld w_sum = 0.0L;
    for (int i = 0; i < n; i++) {
        weights[i] = (w_copy[i] + w_copy[n - 1 - i]) / 2.0L;
        roots[i] = (r_copy[i] - r_copy[n - 1 - i]) / 2.0L;
        w_sum += weights[i];
    }

    // scale weights to get the right value
    for (int i = 0; i < n; i++) {
        weights[i] *= 2.0L / w_sum;
    }

    return {roots, weights};
}

////////////////////////////////////////////////////////////////////////////////

long double f(long double x, long double y) {
    return 3 * sin(8 * PI * x) * cos(8 * PI * y) + x + y + 1;
}

/**
 * @brief 計算 cell 的 Gauss-Legendre 積分
 *
 * 使用 Gauss-Legendre 積分法對區域 [a,b] x [c,d] 進行二重積分。
 * 這個函式將高斯節點映射到目標區域，並根據 Legendre 多項式的根與權重計算積分結果。
 *
 * 由於 Legendre 多項式的根和權重會重複使用，故以參考傳遞方式傳入
 *
 * @param a 起始 x 值 (區間 [a, b] 的左端點)
 * @param b 結束 x 值 (區間 [a, b] 的右端點)
 * @param c 起始 y 值 (區間 [c, d] 的下端點)
 * @param d 結束 y 值 (區間 [c, d] 的上端點)
 * @param roots Legendre 多項式 P_n(x) 的根
 * @param weights Legendre 多項式 P_n(x) 的權重
 * @return ld 積分結果
 */
ld integrate_cell_gauss_legendre(ld a, ld b, ld c, ld d, std::vector<ld>& roots, std::vector<ld>& weights) {
    int n = roots.size();  // 根的數量
    ld integral = 0.0L;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // 將 Gauss 節點轉換為區域 [a,b] x [c,d] 上的座標
            ld x = ((b - a) / 2.0L) * roots[i] + (a + b) / 2.0L;
            ld y = ((d - c) / 2.0L) * roots[j] + (c + d) / 2.0L;

            // 累加對應的權重和函數值
            integral += weights[i] * weights[j] * f(x, y);
        }
    }

    // 乘上 Jacobian 行列式
    return integral * ((b - a) / 2.0L) * ((d - c) / 2.0L);
}

/**
 * @brief 使用 Gaussian Quadrature 計算雙重積分
 *
 * 利用 Gauss-Legendre Quadrature 對區域 [a,b] x [c,d] 上的函數進行雙重積分。
 * 首先將區域劃分為小網格，然後對每個小區域 (cell) 使用高斯積分進行計算。
 *
 * @param a 起始 x 值 (區間的左端點)
 * @param b 結束 x 值 (區間的右端點)
 * @param c 起始 y 值 (區間的下端點)
 * @param d 結束 y 值 (區間的上端點)
 * @param n Legendre 多項式的次數
 * @param mesh_size 網格大小
 * @return ld 計算得到的雙重積分結果
 */
ld integrate_gauss_legendre(ld a, ld b, ld c, ld d, int n, int mesh_size) {
    // 取得 Legendre 多項式 P_n(x) 的根與權重
    auto [roots, weights] = legendre_roots_weights(n);

    // 計算每個網格單元的大小
    ld hx = (b - a) / mesh_size;
    ld hy = (d - c) / mesh_size;

    int cells_x = (b - a) / hx;  // x 方向上分割的 cell 數量
    int cells_y = (d - c) / hy;  // y 方向上分割的 cell 數量

    ld integral = 0.0L;

    // 對每個 cell 進行積分並累加結果
    for (int i = 0; i < cells_x; i++) {
        for (int j = 0; j < cells_y; j++) {
            // 計算每個 cell 的邊界
            ld cell_a = a + i * hx, cell_b = cell_a + hx;
            ld cell_c = c + j * hy, cell_d = cell_c + hy;

            // 計算該 cell 的積分並累加
            integral += integrate_cell_gauss_legendre(cell_a, cell_b, cell_c, cell_d, roots, weights);
        }
    }

    return integral;
}

int main() {
    // std::cout << std::fixed << std::setprecision(15);

    std::vector<int> mash_sizes = {1, 2, 4, 8, 16, 32, 64};
    std::vector<int> n_values = {0, 1, 2, 3, 4, 5};
    for (auto&& mash_size : mash_sizes) {
        for (auto&& n : n_values) {
            ld result;
            auto start = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < 100; i++) {
                result = integrate_gauss_legendre(2.0L, 6.0L, 2.0L, 6.0L, n, mash_size);
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = end - start;
            ld err = result - 144.0L;  // 已知解析解為 144

            std::cout << "mesh_size = " << mash_size << ", n = " << n << ", result = " << result << ", err = " << err
                      << ", " << err / 144.0L << ", time = " << elapsed.count() << "s" << std::endl;

            // std::cout << elapsed.count() << "s" << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
    return 0;
}
