import sys
import numpy as np
import pandas as pd
import genepattern
import gp
from gp.data import _obtain_io, _apply_backwards_compatibility
from scipy import stats
from cuzcatlan import elemental
#import ccalnoir
#from ccalnoir import get_file_from_server


# t test for detecting differentially expressed genes
def expression_t_test(num_rows, info, data, p):

    # use a dictionary to store all uniquely expressed clusters
    # the key is index of probe, and the value is the index of cluster
    # a tuple (_which_cluster_, _p_value_) is stored in each entry
    diff_expression = {}

    # Find uniquely expressed/repressed genes of the provided clusters
    for i in range(len(info)):
        start = info[i][0]
        length = info[i][1]

        # for each gene, flag the current cluster as test, and others as
        # compared calculate the mean for each and do a two-sample t-test to
        # see if the cluster gene expression is uniquely expressed
        for j in range(num_rows):
            test = np.copy(data[j][start:start+length])
            compared = np.copy(data[j])
            for item in test:
                compared = np.delete(compared, item)
            
            # if the expression level is statistically significant (0.05),
            # both tails, then include this in the different expressed clusters
            p_value = stats.ttest_ind(test, compared)[1]
            if (p_value < p):

                # if there is already another cluster has the same differently
                # expressed gene, check the p_value and the smaller one wins
                if j in diff_expression:
                    old_p_value = diff_expression[j][1]
                    if p_value < old_p_value:
                        diff_expression[j] = (i, p_value)
                
                # else, we annotate this diff expressed gene to the
                # corresponding cluster
                else:
                    diff_expression[j] = (i, p_value)

    return diff_expression


# main function for sorting genes based on differentially expressed level
def kmSortFile(data, num_clusters, P_VALUE=0.005, report_only_diff_expr_genes=1):

    # original data with combined clusters
    #original = open(name+'.gct')

################

    # get the input filename and job number
    #jobNum = name.split("/")[-2]
    #input_file_Name = name.split("/")[-1]


    # get the GenePattern input job object and my username
    #lastJob = gp.GPJob(genepattern.get_session(0), jobNum)
    #myUserId = genepattern.get_session(0).username

    # Handle all the various initialization types and get an IO object
    #file_io = _obtain_io(lastJob.get_file(input_file_Name))
    #data = pd.read_table(file_io, skiprows=2)
    #print(data)

#################

    original = data[0]

    #clusters = []
    #for i in range(1,num_clusters+1):
        #with open(name+'-'+str(i)+'.gct') as myfile:
            #head = [next(myfile) for x in range(2)]
            #clusters.append(head)
        #clusters.append(open(name + '-' + str(i) + '.gct'))

    #temp = original.read().split('\n')
    num_rows = original.shape[0]
    num_columns = original.shape[1]-2

    #temp_data = []      # temporarily holds data for future numpy format

    # pull down the data to different chunks.
    # 1) data part -> matrix
    #for i in range(3, len(temp)):
    #    if temp[i] != '':
    #        processing = temp[i].split('\t')
    #        temp_data.append(list(processing[2:]))
    #        del processing
    #matrix = np.array(temp_data, dtype=float)
    matrix = original.values[:,2:]

    # get clusters info (starting location and num of items in each cluster)
    clusters_info = []
    starting_loc = 0
    for i in range(1, num_clusters+1):
        temp_num = data[i].shape[1]-2
        clusters_info.append((starting_loc, temp_num))
        starting_loc += temp_num
    
    #print(clusters_info)
    #return

    diff_expression = expression_t_test(num_rows, clusters_info, matrix, P_VALUE)

    # construct a int array containing length of num of rows, start writing in
    # rows that are differently expressed by each cluster one by one,
    # simultaneously remove them from this array. Finally, write the rest rows
    # IN ORDER to the file
    rest_rows = [x for x in range(num_rows)]

    df_to_process = original
    new_order = []

    # rearrange the rows to pass the matrix for an np ndarray, then construct
    # pandas DataFrame
    for key, value in sorted(diff_expression.items(), key=lambda item :(item[1][0],item[0]) ):
        rest_rows.remove(key)
        new_order.append(key)

    certainty = len(new_order)  # the ones that we know have certainty of diff expressed
    for item in rest_rows:
        new_order.append(item)

    # reorder the pandas table by the new order, then change the indexing back
    # to 'Name' column
    df_to_process = df_to_process.reindex(new_order)
    df_to_process.set_index('Name',inplace=True)

    # report whole genes or only the differentially expressed ones
    if report_only_diff_expr_genes:
        elemental.df2gct(df_to_process.head(certainty), 1, True, name+'-sorted.gct', False)
    else:
        elemental.df2gct(df_to_process, 1, True, name+'-sorted.gct', False)

# test line for individually calling the .py file with command line inputs
#kmSortFile(sys.argv[1], int(sys.argv[2]))#, float(sys.argv[3]), int(sys.argv[4]))

