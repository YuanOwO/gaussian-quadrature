import json
from subprocess import PIPE, Popen, run


def compile_code():
    run(["g++", "-O2", "legendre.cpp"], check=True)


def execute_program(t, n):
    proc = Popen(["./a.exe"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate(input=f"{t} {n}\n".encode())

    # errors = list(map(float, stdout.decode().strip().splitlines()))
    errors = stdout.decode().strip().splitlines()
    return errors


if __name__ == "__main__":
    compile_code()

    results = {}

    for n in [2, 3, 4]:
        print(f"Running for n={n}...")
        results[n] = execute_program("n", n)

    with open("results_n.json", "w") as f:
        json.dump(results, f, indent=2)

    for m in [4, 16, 64]:
        print(f"Running for m={m}...")
        results[m] = execute_program("m", m)

    with open("results_m.json", "w") as f:
        json.dump(results, f, indent=2)
