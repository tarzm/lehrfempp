basis_1d = ["1", "y", "(1.5y^2 - 0.5)", "x", "x * y", "(1.5y^2 - 0.5) * x", "(1.5x^2 - 0.5)", "(1.5x^2 -0.5) * y", "(1.5x^2 - 0.5) * (1.5y^2 - 0.5)"]


# for j in range(0,len(basis_1d)):
#     for i in range(j, len(basis_1d)):
#         mij_string = f'{basis_1d[i]} * {basis_1d[j]}'
#         print (f'Entry ({i}, {j}) is {mij_string}')
    
#     print("")

exp_str = "e^(x*y)"

for i in range(0, len(basis_1d)):
    print(f'Entry {i} i is {exp_str} * {basis_1d[i]}')


