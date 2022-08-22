import numpy as np
import pandas as pd
import string

def string_length(len):
    lower_case = string.ascii_lowercase
    upper_case = string.ascii_uppercase
    numbers = string.printable[0:10]
    special_chars = string.printable[62:94]
    total = list(upper_case + lower_case + numbers+special_chars)
    total.pop(75)
    if len <= 93:
        total = total[0:len]
        return total
    else:
        return total


def subtitution(partitions, save_text_bool,filename):
    protein_encoding = string_length(partitions)
    a = np.full((1, partitions), 1.0) #
    d = np.diag(a[0])                      # creating a matrix of -1s, with 1s as it's diagonal elements
    d[d==0] = -1                           #
    indexes_list = list(protein_encoding)
    df = pd.DataFrame(data=d, index = indexes_list, columns = indexes_list) #for ascii char indexes: example->  !  "  #  $  %  &  '  (  )  *  +  ...  s  t  u  v  w  x  y  z  {  |  }  ~
    print(df)
    if save_text_bool:
        base_filename = filename
        with open('data/'+base_filename,'w') as outfile:
            df.to_string(outfile)


partitions =80
file_name = 'substitution_matrices/synthetic_substitution_matrix_'+str(partitions)+'_chars.txt'
save_text_bool = 1
subtitution(partitions, save_text_bool, file_name)
