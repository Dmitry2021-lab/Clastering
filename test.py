import random
import ledclstr


def main():
    Clster = ledclstr.TLindeBuzoGray()
    em = ledclstr.TGaussianMixture()

    nm = 1000
    vsize = 11
    data = [[random.random() for _ in range(vsize)] for _ in range(nm)]

    print(data[0][2])
    em.EM_Alg(data, 16)
    Clster.RandomInitPoint(data, 4)
    result = Clster.LindeBuzoGray_N(data, 16)
    print(result)


if __name__ == "__main__":
    main()