#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#define ld long double

const long double PI = acosl(-1);

long double f(long double x, long double y) {
    return 3 * sin(8 * PI * x) * cos(8 * PI * y) + x + y + 1;
}

/**
 * @brief 使用梯形法則計算 f(x, y) 在區域 [a, b] x [c, d] 上使用 n x m 分割的雙重積分
 *
 * @param a 區域起始 x 值
 * @param b 區域結束 x 值
 * @param c 區域起始 y 值
 * @param d 區域結束 y 值
 * @param n x 方向分割數
 * @param m y 方向分割數
 * @return long double 雙重積分結果
 */
long double trapezoid(long double a, long double b, long double c, long double d, int n, int m) {
    long double hx = (b - a) / n, hy = (d - c) / m;
    long double integral = 0.0L;

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            long double x = a + i * hx, y = c + j * hy;
            long double c = 4.0L;  // 預設為內部點

            if (i == 0 || i == n) c /= 2.0L;  // 在 x 邊界上
            if (j == 0 || j == m) c /= 2.0L;  // 在 y 邊界上

            integral += c * f(x, y);
        }
    }

    return integral * (hx * hy / 4.0L);
}

int main() {
    std::cout << "long double precision: " << LDBL_DIG << " digits" << std::endl;
    std::cout << std::fixed << std::setprecision(15);

    for (int n = 4; n <= 128; n *= 2) {
        long double result = trapezoid(2, 6, 2, 6, n, n);
        long double err = result - 144.0L;  // 已知解析解為 144
        std::cout << "n = " << n << ", result = " << result << ", err = " << err << ", " << err / 144.0L << std::endl;
    }

    return 0;
}
