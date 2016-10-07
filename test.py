import soap

at = [[1],[7]]
c = [[[0.1,0.1,0.1]],[[0.2,1,0.6]]]

M = soap.kernelf(at,c,at,c,True)
print(M)
