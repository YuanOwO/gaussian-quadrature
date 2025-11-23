#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>

using ld = long double;

const ld PI = acosl(-1);

ld f(ld x, ld y) {
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
 * @return ld 雙重積分結果
 */
ld trapezoid(ld a, ld b, ld c, ld d, int n, int m) {
    ld hx = (b - a) / n, hy = (d - c) / m;
    ld integral = 0.0L;

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            ld x = a + i * hx, y = c + j * hy;
            ld c = 4.0L;  // 預設為內部點

            if (i == 0 || i == n) c /= 2.0L;  // 在 x 邊界上
            if (j == 0 || j == m) c /= 2.0L;  // 在 y 邊界上

            integral += c * f(x, y);
        }
    }

    return integral * (hx * hy / 4.0L);
}

int main() {
    std::cout << "long double precision: " << LDBL_DIG << " digits" << std::endl;

    for (int n = 4; n <= 128; n *= 2) {
        ld result = trapezoid(2, 6, 2, 6, n, n);
        ld err = result - 144.0L;  // 已知解析解為 144
        std::cout << std::fixed << std::setprecision(LDBL_DIG - 3) << n << " & $" << result << "$ & ";
        std::cout << std::defaultfloat << std::setprecision(6) << "$\\sci{" << err << "}$ & $\\sci{" << err / 144.0L
                  << "}$ \\\\" << std::endl;
    }

    return 0;
}
