import numpy as np

def terminal_branch(Br_El):#Br_El is el_array output from import_exelem_tree subroutine

    T_Br=np.ones((len(Br_El),1))#if it is 1, it is terminal br, if it is zero it is not, initialise all branches as terminal

    for i in range(0, len(Br_El)):
        for j in range(0,len(Br_El)):     
        
            if Br_El[j][1] == Br_El[i][2]:# if the first node of an elment is equal to 2nd node another element 
               T_Br[i] = 0 #the element is not terminal (zero means 'not terminal')
             
       
    Terminal_El_num=np.array(np.where(T_Br[:,0]==1))#looking for the element number which are terminal branch 
    
    return {'terminal_el': Terminal_El_num, 'total_terminal_el': len(Terminal_El_num[0])}

