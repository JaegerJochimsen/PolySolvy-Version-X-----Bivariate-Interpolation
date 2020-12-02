import math

"""
array of numbers (repr method ouotputs formatted string)

- generate coeff via vandermond, put these in array
"""
"""
def funct(...):
    generate our function 
    return function obj
"""


def populate(num_points):
    """
    Takes a squaer number of points (i.e. 1,4,9,etc.) and populates two lists, one of y terms and the other of x terms
    with each term having exponents beginning at 0 and increasing to sqrt(num_points).

    Ex/ num_points = 4
    -> xterms = [x^0, x^1, x^2],
    -> yterms = [y^0, y^1, y^2]

    param: num_points: the square-number of points for which to find the interpolating polynomial for.

    return: List[xterms, yterms] -- a list of 2 lists, where xterms and yterms are as depicted in the above example.
    """
    superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    highest_deg = int(math.sqrt(num_points))
    x_terms = []
    y_terms = []
    for i in range(highest_deg):
        x_terms.append(f'(x{i})'.translate(superscript))
        y_terms.append(f'(y{i})'.translate(superscript))
    return [x_terms, y_terms]

def make_order_right(string):
    first = string[:4]
    second = string[4:]
    return f'{second}{first}'

def generate_generic_form(num_points):
    subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    x_terms = populate(num_points)[0]
    y_terms = populate(num_points)[1]
    combined_terms =[]
    for termx in x_terms:
        for termy in y_terms:
            combined_terms.append(f'{termx}{termy}')

    const = 0
    for term_index in range(len(combined_terms)):
        if combined_terms[term_index][2] == "⁰":            # if x term has 0 exponent
            if combined_terms[term_index][6] == "⁰":        # if y term ALSO has 0 exponent
                combined_terms[term_index] = "c0"        # replace with appropriate constant
            else:                                           # if just the x term is ^0
                combined_terms[term_index] = combined_terms[term_index].replace("(x⁰)", f'c{const}')

        elif combined_terms[term_index][6] == "⁰":          # if just the y term is ^0
            combined_terms[term_index] = combined_terms[term_index].replace(f"(y⁰)", f"c{const}")
            combined_terms[term_index] = make_order_right(combined_terms[term_index])

        else:
            combined_terms[term_index] = f'c{const}{combined_terms[term_index]}'
        const += 1
    return " + ".join(combined_terms).translate(subscript)


def main():
    print("General function using 4 nodes: ", generate_generic_form(4))
    print("General functio using 9 nodes: ", generate_generic_form(9))
    print("General functio using 16 nodes: ", generate_generic_form(16))

if __name__ == "__main__":
    main()


