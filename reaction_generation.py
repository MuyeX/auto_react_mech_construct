from collections import Counter

def get_unique_intermediates(count, products, reactant):
    # Generate a list of unique intermediate identifiers
    # Check what is the last used letter in input
    if isinstance(reactant, str):
        reactant = [reactant]

    if isinstance(products, Counter):
        products = list(set(products)) 

    input_list = products + reactant
    last_char = sorted(input_list)[-1]
    next_char_num = ord(last_char)+1

    # Generate intermediates starting from the next char
    return [chr(next_char_num + i) for i in range(count)]

def create_reactions(reactant, products, intermediates):
    reactions = []
    remaining_products = Counter(products)
    #TODO:currently it is only checking for 'B' but it could be a different char based on the input, fix it

    # Generate the initial reaction with the first intermediate
    if remaining_products['B'] > 0:
        reactions.append(f"{reactant} -> B + {intermediates[0]}")
        remaining_products['B'] -= 1

    # Generate reactions involving intermediates
    for i in range(len(intermediates)):
        # If this is the last intermediate, it should produce the final products
        if i == len(intermediates) - 1:
            product_str = ' + '.join([product for product in remaining_products.elements()])
            reactions.append(f"{intermediates[i]} -> {product_str}")
        else:
            # If 'B' is still needed and we're not at the last intermediate
            if remaining_products['B'] > 0:
                reactions.append(f"{intermediates[i]} -> B + {intermediates[i+1]}")
                remaining_products['B'] -= 1
            else:
                # If all 'B's are accounted for, just produce the next intermediate
                reactions.append(f"{intermediates[i]} -> {intermediates[i+1]}")

    return reactions
