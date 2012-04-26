import sympy as sp
import numpy as np
import inspect
import copy

INF = np.inf

def lie_bracket(f1, f2, dx):
    
    n = len(dx)
    f1dx = sp.matrices.zeros([n,n])
    f2dx = sp.matrices.zeros([n,n])
    
    for i in range(n):
        for j in range(n):
            f1dx[i,j] = sp.diff(f1[i], dx[j])
            f2dx[i,j] = sp.diff(f2[i], dx[j])
        
    lb = (f2dx*f1 - f1dx*f2).expand()
    
    return lb
    
    
def lie_algebra(F, dx, n_iterations=INF, show=True):
    # F should either be a list of sp.Matrix elements corresponding to the vector fields, or a sp.Matrix where each column is a vector field
    # dx should be a list of the sp.symbols that correspond to the partial derivatives

    # Make a matrix of Vector Fields
    if type(F) is list:
        n = len(F)
        newF = F[0]
        for i in range(1, n):
            newF = newF.row_join(F[i])
        F = newF
    else:
        newF = F
        
    # number of original vector fields
    nF = F.shape[1]
        
    # original strings:
    F_strings = ['f'+str(i) for i in range(nF)]

    # initialize lie algebra dictionaries
    LA_zero = {} # lie brackets that equal zero
    LA_rejected = {} # lie brackets that are linearly dependent on other terms
    LA_complete = {} # the complete lie algebra
    LA_original = {} # the original vector fields
    tested_lie_brackets = [] # keep track of all the computed lie brackets in string form
        
    # match vector fields with their string names
    for i, s in enumerate(F_strings):
        LA_original.setdefault(F_strings[i], F[:,i])
    # initialize complete Lie Algebra with original
    LA_complete = copy.copy(LA_original)
    
    # temporary lie algebra, used for linear dependence calculations
    newF_tmp = newF
    
    # initialize    
    have_new_lie_brackets = True
    iterations = 0
    
    # run loop
    while have_new_lie_brackets is True and iterations < n_iterations:
        print 'ITERATION: ', iterations, 'out of: ', n_iterations
        
        # re-initialize
        have_new_lie_brackets = False
        iterations += 1
        n = newF.shape[1]
        
        original_vector_field_strings = LA_original.keys()
        lie_algebra_strings = LA_complete.keys()
        
        for i in range(len(original_vector_field_strings)):
            for j in range(len(lie_algebra_strings)):
            
                k1 = original_vector_field_strings[i]
                k2 = lie_algebra_strings[j]
            
                print
                #print i,j
                
                lie_bracket_string_name = '[' + k1 + ', ' + k2 + ']'
                print 'calculating: '
                print lie_bracket_string_name

                if lie_bracket_string_name in tested_lie_brackets:
                    print 'already computed'
                    continue
                else:
                    tested_lie_brackets.append(lie_bracket_string_name)
                
                # hack to make sure f0 term is on the left:
                if k2 == 'f0':
                    print 'skipping - will get covered another time'
                    continue
                
                lb = lie_bracket(LA_original[k1], LA_complete[k2], dx)
                if np.sum(np.abs(lb)) == 0:
                    print 'lie bracket equals 0'
                    LA_zero.setdefault(lie_bracket_string_name, lb)
                    continue
                    
                # temporarily add the new lie bracket to our lie algebra    
                newF_tmp = lb.row_join(newF)
                
                # check to see if there is any linear dependence of the first term on the other terms
                nullspace = newF_tmp.nullspace()
                for ii in range(len(nullspace)):
                    if nullspace[ii][0] != 0:
                        nullspace[ii] *= nullspace[ii][0]**-1
                nullspace_linear_components = [is_vec_linear(nullspace[ii], dx) for ii in range(len(nullspace))]
                
                
                # if the nullspace is linear, then we do not have a linearly independent lie bracket, else we do
                if np.sum(nullspace_linear_components) == 0:
                    newF = newF_tmp
                    have_new_lie_brackets = True
                    LA_complete.setdefault(lie_bracket_string_name, lb)
                    print '* Added to Lie Algebra'
                else:
                    LA_rejected.setdefault(lie_bracket_string_name, lb)
                    print '* Linearly Dependent on other terms'
                    
    if show:
        print '*'*25
        print 'showing results for ', n_iterations, ' iterations'
    
        print '*'*25
        print 'original vector fields: '
        print_LA_simple(LA_original)
                    
        print '*'*25
        print 'complete lie algebra: '
        print_LA_simple(LA_complete)
        
        print '*'*25
        print 'linearly dependent lie brackets: '
        print_LA_simple(LA_rejected)
        
        print '*'*25
        print 'lie brackets = zero: '
        print_LA_simple(LA_zero)
                    
    return LA_original, LA_complete, LA_rejected, LA_zero
    
def print_LA_simple(LA):
    # shows the values of the dictionary in a practical way

    keys = LA.keys()
    keys.sort(key = len) # sort shortest to longest
    
    for key in keys:
        print
        print key, ': '
        print LA[key]
        
    
def is_term_linear(m, var_list):
    
    if 0: # derivative method: too slow
        typsok = []
        '''
        for name, obj in inspect.getmembers(sp.core.numbers):
            if inspect.isclass(obj):
                typsok.append(obj)
        '''
        typsok = [sp.core.numbers.Real, sp.core.numbers.Zero, sp.core.numbers.One, sp.core.numbers.Pi, sp.core.numbers.Integer, sp.core.numbers.Rational]
        
        typs = []
        for var in var_list:
            typ = type(sp.diff(m, var, 1).evalf())
            if typ in typsok:
                typs.append(1)
            else:
                typs.append(0)
        if np.sum(typs) == len(typs):
            return True
        else:
            return False
            
    if 1: # poly method
            
        try:
            p = m.as_poly()
            
            if p is None:
                return False
            else:
                deg = p.degree
                if deg <= 1:
                    return True
                else:
                    return False 
            
        except (sp.SymbolsError, AttributeError):
            return True
            
        except:
            print 'help!'
            return False
            
            
        
        
def is_vec_linear(vec, var_list):
    typs = []
    for v in vec:
        typs.append(is_term_linear(v, var_list))
    if np.sum(typs) == len(typs):
        return True
    else:
        return False
    

    
if __name__ == "__main__":
    # example

    w1,w2,w3 = sp.symbols('w1,w2,w3')
    g1, g2, g3 = sp.symbols('g1,g2,g3')

    f0 = sp.Matrix([g1*w2*w3, g2*w1*w3, g3*w1*w2])
    f1 = sp.Matrix([1, 0, 0])
    f2 = sp.Matrix([0, 1, 0])
    dx = [w1, w2, w3]

    F = [f0, f1, f2]

    LA_original, LA_complete, LA_rejected, LA_zero = lbt.lie_algebra(F, dx, 2, show=True)



    
