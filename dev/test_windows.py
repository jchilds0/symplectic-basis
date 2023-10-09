from datetime import datetime
import snappy
from multiprocessing import Pool
import itertools

import symplectic_basis

from test_base import is_symplectic, testing_string
import random

start = 0
end = 1
scale = 1000
test = "random"

if len(snappy.HTLinkExteriors(crossings=15)) == 0:
    file_name = "small-db"
else:
    file_name = "large-db"

test_set = snappy.HTLinkExteriors(knots_vs_links="links")


def process_manifold(index: int, output: bool = True):
    if test == "random":
        index = random.randint(1, len(snappy.HTLinkExteriors(knots_vs_links="links")) - 1)

    M = snappy.HTLinkExteriors(knots_vs_links="links")[index]
    label = M.identify()[0] if len(M.identify()) > 0 else ""
    # print(label)

    if index == 0:
        return True

    basis = symplectic_basis.symplectic_basis(M)
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    if output:
        with open("logs/links-0.log", "a") as file:
            file.write(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}\n")

    print(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}")

    return result


def test_link_complements_pool(manifolds):
    with open("logs/total.log", "a") as file:
        if test == "database":
            length = len(manifolds)
        else:
            length = scale * (end - start)

        file.write(testing_string(length))
        print(testing_string(length))

    with Pool(maxtasksperchild=25) as pool:
        if test == "database":
            result = pool.imap(process_manifold, manifolds)
        elif test == "random":
            result = pool.imap(process_manifold, [random.randint(1, len(test_set)) for _ in range((end - start) * scale)])
        else:
            result = pool.imap(process_manifold, range(start * scale, end * scale))

        for _ in range(start, end):
            lst = list(itertools.islice(result, scale))

        # lst = list(result)
            time = datetime.now().strftime('%d-%m-%y %H:%M:%S')
            print(f"[{time}]    Passed: {sum(lst)} / {len(lst)}")

            with open("logs/total.log", "a") as file:
                file.write(f"[{time}]    Passed: {sum(lst)} / {len(lst)}\n")


if __name__ == "__main__":
    # with open(file_name, "r") as file:
    #     lst = file.readlines()

    # manifolds = list(set([ for x in lst]))
    test_link_complements_pool([])
    # test_link_complements()
    # generate_tests()
    # unittest.main()
