from subprocess import PIPE, Popen, run

from numpy.polynomial.legendre import leggauss

EPS = 1e-13


def test_precision(n: int):
    x, w = leggauss(n)
    proc = Popen(["./a.exe"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    input_data = f"{n}\n".encode()
    stdout, stderr = proc.communicate(input=input_data)

    roots_weights = stdout.decode().strip().split("\n")

    for i, line in enumerate(roots_weights):
        x_wave, w_wave = map(float, line.split())

        if abs(x[i]) < EPS:
            error = abs(x_wave)
        else:
            error = abs(x_wave - x[i]) / abs(x[i])

        assert error < EPS, f"Root mismatch at index {i}: {x_wave} vs {x[i]}"

        if abs(w[i]) < EPS:
            error = abs(w_wave)
        else:
            error = abs(w_wave - w[i]) / abs(w[i])

        assert error < EPS, f"Weight mismatch at index {i}: {w_wave} vs {w[i]}"

    print(f"Test passed for n={n}")


if __name__ == "__main__":
    run(["g++", "-O2", "legendre.cpp"])

    for n in range(1, 129):
        test_precision(n)
    # xx = []
    # ww = []
    # for n in range(1, 129):
    #     x, w = leggauss(n)
    #     xx.append(x.tolist())
    #     ww.append(w.tolist())

    # with open("legendre_roots_weights.json", "w") as f:
    #     import json

    #     json.dump({"roots": xx, "weights": ww}, f)
