import sympy as sp
import numpy as np
import inspect

INF = np.inf

def lie_bracket(f1, f2, dx):
    
    n = len(dx)
    f1dx = sp.matrices.zeros([2,2])
    f2dx = sp.matrices.zeros([2,2])
    
    for i in range(n):
        for j in range(n):
            f1dx[i,j] = sp.diff(f1[i], dx[j])
            f2dx[i,j] = sp.diff(f2[i], dx[j])
        
    lb = f2dx*f1 - f1dx*f2
    
    return lb
    
    
def lie_algebra(F, dx, n_iterations=INF):
    
    # Make a matrix of Vector Fields
    if type(F) is list:
        n = len(F)
        newF = F[0]
        for i in range(1, n):
            newF = newF.row_join(F[i])
    else:
        newF = F
    
    newF_tmp = newF
    newF_strings = ['f'+str(i) for i in range(n)]
    
    failedF = []
    failedF_strings = []
    
    have_new_lie_brackets = True
    iterations = 0
    while have_new_lie_brackets is True and iterations < n_iterations:
        have_new_lie_brackets = False
        iterations += 1
        n = newF.shape[1]
        for i in range(n):
            for j in range(i+1, n):
                print
                print i,j
                
                lie_bracket_string_name = '[' + newF_strings[i] + ', ' + newF_strings[j] + ']'
                if lie_bracket_string_name in newF_strings or lie_bracket_string_name in failedF_strings:
                    continue
                
                lb = lie_bracket(newF[:,i], newF[:,j], dx)
                if np.sum(lb) == 0:
                    continue
                    
                print lb
                
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
                    newF_strings.insert(0, lie_bracket_string_name)
                    have_new_lie_brackets = True
                    print 'yes'
                else:
                    failedF.insert(0,lb)
                    failedF_strings.insert(0,lie_bracket_string_name)
                    print 'no'
                    
    return newF, newF_strings, failedF, failedF_strings
    
    
    
    
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
    x,y,z = sp.symbols('x,y,z')
    dx = [x,z]
    f1 = sp.Matrix([0, x**2])
    f2 = sp.Matrix([z+x, 0])
    
    F = [f1, f2]
    newF, newF_strings, failedF, failedF_strings = lie_algebra(F, dx, 2)
    print 
    print
    print 'Lie Algebra:'
    print 
    print newF_strings
    print newF
    print
    print 'Failed Lie Brackets:'
    print failedF_strings
    print failedF
