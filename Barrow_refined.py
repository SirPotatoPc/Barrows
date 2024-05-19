from sympy import binomial, Symbol, Sum, Matrix


n = Symbol('n') #total amount of items interested in
k = Symbol('k') #variable used to sum over. the amount of items interested
#  in during intermediate calculations


#R=102, d=24. It was done this way to not have to do fractions in the matrix
R = Symbol('R') #1/rate of items
d = Symbol('d') #total amount of items

#Transition matrix. last row represents probability of getting at least one
# of 'n' interesting items
M = Matrix([
[1-1/R,0,0,0,0,0,0,0],
[(d-k)/(d*R),1-1/R,0,0,0,0,0,0],
[0,(d-1-k)/((d-1)*R),1-1/R,0,0,0,0,0],
[0,0,(d-2-k)/((d-2)*R),1-1/R,0,0,0,0],
[0,0,0,(d-3-k)/((d-3)*R),1-1/R,0,0,0],
[0,0,0,0,(d-4-k)/((d-4)*R),1-1/R,0,0],
[0,0,0,0,0,(d-5-k)/((d-5)*R),1-1/R +((d-6)-k)/((d-6)*R),0],
[k/(d*R),k/((d-1)*R),k/((d-2)*R),k/((d-3)*R),k/((d-4)*R),k/((d-5)*R),k/((d-6)*R),1]
]).subs({R:102,d:24})



chests = Symbol('N') #How many chests opened

# |A_k|
Prob = Sum((-1)**k*binomial(n,k)*((1-(M**7)[-1,0])**chests),(k,0,n))


amount_of_sets = Symbol('I')
# |S1 u S2 u ... u S6|
ProbAny = Sum((-1)**(amount_of_sets+1) * binomial(6, amount_of_sets)* Prob.subs(n,amount_of_sets*4), (amount_of_sets, 1, 6))



#Calculating probabilities for specific amounts of chests
amount_of_chests = [1,10,25,50,75,100,250,500,750,1000,2000]
print(f"    N | {'one':^8} | {'all':^8} | {'any':^8} |")
for N in amount_of_chests:
    One = float(Prob.evalf(8,subs={chests:N,n:4}))
    All = float(Prob.evalf(8,subs={chests:N,n:24}))
    Any = float(ProbAny.evalf(8,subs={chests:N}))
    print (f"{N:>5} | {One:<8.{'2E' if One<1e-6 else '6f'}} | {All:<8.{'2E' if All<1e-6 else '6f'}} | {Any:<8.{'2E' if Any<1e-6 else '6f'}} |" )