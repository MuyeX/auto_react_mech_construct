from collections import Counter
from copy import deepcopy

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
    return [chr(next_char_num + i) for i in range(count - 1)]

def create_reactions(reactant, products, intermediates):
    reactions = []
    remaining_products = Counter(products)
    #TODO:currently it is only checking for 'B' but it could be a different 
    #char based on the input, fix it

    # Generate the initial reaction with the first intermediate
    if remaining_products['B'] > 0:
        reactions.append(f"{reactant} -> B + {intermediates[0]}")
        remaining_products['B'] -= 1

    # Generate reactions involving intermediates
    for i in range(len(intermediates)):
        # If this is the last intermediate, it should produce the final products
        if i == len(intermediates) - 1:
            product_str = ' + '.join(
                [product for product in remaining_products.elements()]
                )
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

def generate_alternative_chains(base_chain, reactant=['A'], products=['B', 'B', 'B', 'C']):

    if isinstance(products, Counter):
        products_list = []
        for key in products.keys():
            l = [key] * products[key]
            products_list += l
        # print("list", products_list)
        products = products_list

    intermediates = []
    for equation in base_chain:
        for char in equation:
            if char not in ['+', '-', '>', ' '] + reactant + products:
                intermediates.append(char)
    
    intermediates = sorted(list(set(intermediates)))
    products=Counter(['B', 'B', 'B', 'C'])
    # print(intermediates)

    min_req_steps = sum(products.values()) - 2
    num_free_steps = len(intermediates) - min_req_steps

    if num_free_steps == 0:
        return []
    num_all_steps = len(base_chain)
    # print(num_free_steps, num_all_steps)
    alternative_chains = []

    for i in range(num_all_steps - num_free_steps):
        chain = []
        if i == 0:
            chain.append((f"{reactant[0]} -> {intermediates[0]}"))
            remaining_products = deepcopy(products)
            intermediate_idx = 1
            while sum(remaining_products.values()) > 0:

                if sum(remaining_products.values()) == 2 and intermediate_idx == len(intermediates):
                    # print("HERE")
                    product1 = list(remaining_products.keys())[0]
                    product2 = list(remaining_products.keys())[1]
                    product_str = f'{product1} + {product2}'
                    chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                    # print(chain)
                    # print(remaining_products)
          
                    break
                
                if sum(remaining_products.values()) == 1 and intermediate_idx == len(intermediates):
                    # print("HERE")
                    product = list(remaining_products.keys())[0]
                    product_str = f'{product}'
                    chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                    # print(chain)
                    # print(remaining_products)
          
                    break

                product = list(remaining_products.keys())[0]
                product_str = f'{product} + {intermediates[intermediate_idx]}'
                chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                remaining_products[product] -= 1
                if remaining_products[product] == 0:
                    del remaining_products[product]
                intermediate_idx += 1
                # print(chain)
                # print(remaining_products)
            
                
        
        if i > 0:
            # print(f"products for {i}",products)
            remaining_products = deepcopy(products)
            product = list(remaining_products.keys())[0]
            chain.append((f"{reactant[0]} -> {product} + {intermediates[0]}"))
            remaining_products[product] -= 1

            intermediate_idx = 1
            while sum(remaining_products.values()) > 0:
                 
                if sum(remaining_products.values()) == 2 and intermediate_idx == len(intermediates):
                    # print("HERE")
                    product1 = list(remaining_products.keys())[0]
                    product2 = list(remaining_products.keys())[1]
                    product_str = f'{product1} + {product2}'
                    chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                    # print(chain)
                    # print(remaining_products)
          
                    break
                
                if sum(remaining_products.values()) == 1 and intermediate_idx == len(intermediates):
                    # print("HERE")
                    product = list(remaining_products.keys())[0]
                    product_str = f'{product}'
                    chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                    # print(chain)
                    # print(remaining_products)
          
                    break

                if intermediate_idx == i:
                    product_str = f'{intermediates[intermediate_idx]}'
                else:
                    product = list(remaining_products.keys())[0]
                    product_str = f'{product} + {intermediates[intermediate_idx]}'
                    remaining_products[product] -= 1
                    if remaining_products[product] == 0:
                        del remaining_products[product]
                
                chain.append(f"{intermediates[intermediate_idx-1]} -> {product_str}")
                intermediate_idx += 1
                # print(chain)
                # print(remaining_products)

        alternative_chains.append(chain)
    
    return alternative_chains

# # Base chain for the reaction A -> 3B + C with three intermediates
# base_chain = [
#     'A -> B + D',
#     'D -> B + E',
#     'E -> B + F',
#     'F -> C'
# ]

# # Generate alternative chains
# alternative_chains = generate_alternative_chains(base_chain)
# for chain in alternative_chains:
#     print(chain)
