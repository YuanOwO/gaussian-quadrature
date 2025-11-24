# HW3: 2D Gaussian Quadrature

使用 C++ 實作二維 Gauss-Legendre 積分法，並比較不同網格劃分與節點數量對積分結果的影響。

## 檔案說明

-   `README.md`：專案說明文件。
-   `legendre.cpp`：主程式碼，包含 Gauss-Legendre 積分的實作與測試。
-   `trapezoid.cpp`：使用梯形法則計算雙重積分的程式碼。
-   `contour.py`：用於繪製函數等高線
-   `draw.py`：用於繪製積分結果的圖表。
-   `run.py`：用於批次執行測試的腳本。
-   `test.py`：用於驗證 legendre 找根與權重的正確性。

## 使用說明

1. 編譯程式碼：

    ```bash
    g++ -O2 -o legendre legendre.cpp
    ```

2. 執行程式：

    ```bash
    ./legendre
    ```

3. 查看輸出結果，包含不同網格劃分與節點數量下的積分結果、誤差與計算時間。
